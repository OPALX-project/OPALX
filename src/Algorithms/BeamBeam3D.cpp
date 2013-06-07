// ------------------------------------------------------------------------
// $RCSfile: BeamBeam3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// Definitions for class: BeamBeam3D
//   Defines a concrete beam-beam interaction.
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/BeamBeam3D.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/ElementImage.h"
#include "Algorithms/PartBunch.h"
#include "BeamlineGeometry/Matrix3D.h"
#include "BeamlineGeometry/Vector3D.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "Physics/Physics.h"
#include "Utilities/ComplexErrorFun.h"
#include "Utilities/Gauss.h"
#include "Utilities/InverseGauss.h"
#include "Utilities/TpsWerrf.h"
#include <cmath>
#include <complex>

using std::complex;
using std::max;
using namespace Physics;


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(BeamBeam3D::*get)() const;
        void (BeamBeam3D::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &BeamBeam3D::getElementLength,
            0
        },
        { 0, 0, 0 }
    };

    template <class T>
    inline T sqr(const T &x) { return x * x; }
}


// Class BeamBeam3D
// ------------------------------------------------------------------------

BeamBeam3D::BeamBeam3D():
    BeamBeam(), geometry(), errorFunction(Werrf), F(0.0), slices(1), nsli(1) {
    lf.betax = lf.betay = lf.emitx = lf.emity = 1;
}


BeamBeam3D::BeamBeam3D(const BeamBeam3D &rhs):
    BeamBeam(rhs), geometry(), errorFunction(Werrf), F(rhs.F),
    slices(rhs.slices), nsli(rhs.nsli) {
    lf.betax = rhs.lf.betax;
    lf.betay = rhs.lf.betay;
    lf.emitx = rhs.lf.emitx;
    lf.emity = rhs.lf.emity;
}


BeamBeam3D::BeamBeam3D(const string &name):
    BeamBeam(name), geometry(), errorFunction(Werrf), F(0.0),
    slices(1), nsli(1) {
    lf.betax = lf.betay = lf.emitx = lf.emity = 1;
}


BeamBeam3D::~BeamBeam3D()
{}


void BeamBeam3D::accept(BeamlineVisitor &visitor) const {
    visitor.visitComponent(*this);
}


ElementBase *BeamBeam3D::clone() const {
    return new BeamBeam3D(*this);
}


Channel *BeamBeam3D::getChannel(const string &aKey, bool create) {
    return ElementBase::getChannel(aKey, create);
}


NullField &BeamBeam3D::getField() {
    return field;
}


const NullField &BeamBeam3D::getField() const {
    return field;
}


NullGeometry &BeamBeam3D::getGeometry() {
    return geometry;
}

const NullGeometry &BeamBeam3D::getGeometry() const {
    return geometry;
}


ElementImage *BeamBeam3D::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


const string &BeamBeam3D::getType() const {
    static const string type("BeamBeam3D");
    return type;
}


double BeamBeam3D::getBunchCharge() const {
    return 0.0;
}


const Matrix3D &BeamBeam3D::getBunchMoment() const {
    static const Matrix3D dummy;
    return dummy;
}


const Vector3D &BeamBeam3D::getBunchDisplacement() const {
    static const Vector3D dummy;
    return dummy;
}


void BeamBeam3D::setErrorFunctionPointer
(complex<double> (*fun)(complex<double>)) {
    errorFunction = fun;
}


void BeamBeam3D::setCrossingAngle(double value) {
    phi = value;
    cphi = cos(phi);
    sphi = sin(phi);
    tphi = sphi / cphi;
    computeSlices();
}


void BeamBeam3D::setBeamBeamParameter(double value) {
    xiyn = value;
    computeF();
}


void BeamBeam3D::setBeamDescription(const Vector3D &disp, const Beta &beta) {
    displacement = disp;
    lf = beta;
    computeF();
    computeSlices();
}


void BeamBeam3D::setSlices(int slices) {
    nsli = slices;
    computeF();
    computeSlices();
}


void BeamBeam3D::trackBunch
(PartBunch &bunch, const PartData &, bool revBeam, bool revTrack) const {
    boost(bunch);
    synchroBeamCollision(bunch);
    boosti(bunch);
    double ax, ay;
    bunch.maximumAmplitudes(D, ax, ay);
    axmax = max(axmax, ax);
    aymax = max(aymax, ay);
}


void BeamBeam3D::trackMap
(FVps<double, 6> &map, const PartData &, bool revBeam, bool revTrack) const {
    boost(map);
    synchroBeamCollision(map);
    boosti(map);
}


void BeamBeam3D::boost(PartBunch &bunch) const {

    for(unsigned int i = 0; i < bunch.getLocalNum(); i++) {
        Particle part = bunch.get_part(i);

        double a = (sqr(part.px()) + sqr(part.py())) / sqr(1.0 + part.pt());
        double sqr1a = sqrt(1.0 - a);
        double h = (part.pt() + 1.0) * a / (1.0 + sqr1a);
        part.px() = (part.px() - tphi * h) / cphi;
        part.py() = part.py() / cphi;
        part.pt() = part.pt() - sphi * part.px();
        double a1 = (sqr(part.px()) + sqr(part.py())) / sqr(1.0 + part.pt());
        sqr1a = sqrt(1.0 - a1);
        double hd1 = (1.0 + part.pt()) * sqr1a;
        double h1x = part.px() / hd1;
        double h1y = part.py() / hd1;
        double h1z = a1 / (1.0 + sqr1a) / sqr1a;

        double x1 = tphi * part.t() + (1.0 + sphi * h1x) * part.x();
        part.y() = part.y() + sphi * h1y * part.x();
        part.t() = part.t() / cphi - sphi * h1z * part.x();
        part.x() = x1;
        bunch.set_part(part, i);
    }
}


void BeamBeam3D::boost(FVps<double, 6> &map) const {
    Series a = (sqr(map[1]) + sqr(map[3])) / sqr(1.0 + map[5]);
    Series sqr1a = sqrt(1.0 - a);
    Series h = (map[5] + 1.0) * a / (1.0 + sqr1a);
    map[1] = (map[1] - tphi * h) / cphi;
    map[3] = map[3] / cphi;
    map[5] = map[5] - sphi * map[1];
    Series a1 = (sqr(map[1]) + sqr(map[3])) / sqr(1.0 + map[5]);
    sqr1a = sqrt(1.0 - a1);
    Series hd1 = (1.0 + map[5]) * sqr1a;
    Series h1x = map[1] / hd1;
    Series h1y = map[3] / hd1;
    Series h1z = a1 / (1.0 + sqr1a) / sqr1a;

    Series x1 = tphi * map[4] + (1.0 + sphi * h1x) * map[0];
    map[2] = map[2] + sphi * h1y * map[0];
    map[4] = map[4] / cphi - sphi * h1z * map[0];
    map[0] = x1;
}


void BeamBeam3D::boosti(PartBunch &bunch) const {
    for(unsigned int i = 0; i < bunch.getLocalNum(); i++) {
        Particle part = bunch.get_part(i);

        double a1 = (sqr(part.px()) + sqr(part.py())) / sqr(1.0 + part.pt());
        double sqr1a = sqrt(1.0 - a1);
        double h1d = (1.0 + part.pt()) * sqr1a;
        double h1 = (part.pt() + 1.0) * a1 / (1.0 + sqr1a);
        double h1x = part.px() / h1d;
        double h1y = part.py() / h1d;
        double h1z = a1 / (1.0 + sqr1a) / sqr1a;
        double det = 1.0 + sphi * (h1x - sphi * h1z);
        part.x() = (part.x() - sphi * part.t()) / det;
        part.t() = cphi * (part.t() + h1z * sphi * part.x());
        part.y() = part.y() - h1y * sphi * part.x();
        part.pt() = part.pt() + sphi * part.px();
        part.px() = (part.px() + sphi * h1) * cphi;
        part.py() = part.py() * cphi;
        bunch.set_part(part, i);
    }
}


void BeamBeam3D::boosti(FVps<double, 6> &map) const {
    Series a1 = (sqr(map[1]) + sqr(map[3])) / sqr(1.0 + map[5]);
    Series sqr1a = sqrt(1.0 - a1);
    Series h1d = (1.0 + map[5]) * sqr1a;
    Series h1 = (map[5] + 1.0) * a1 / (1.0 + sqr1a);
    Series h1x = map[1] / h1d;
    Series h1y = map[3] / h1d;
    Series h1z = a1 / (1.0 + sqr1a) / sqr1a;
    Series det = 1.0 + sphi * (h1x - sphi * h1z);
    map[0] = (map[0] - sphi * map[4]) / det;
    map[4] = cphi * (map[4] + h1z * sphi * map[0]);
    map[2] = map[2] - h1y * sphi * map[0];
    map[5] = map[5] + sphi * map[1];
    map[1] = (map[1] + sphi * h1) * cphi;
    map[3] = map[3] * cphi;
}


void BeamBeam3D::computeF() {
    double sigxxn = lf.emitx * lf.betax;
    double sigyyn = lf.emity * lf.betay;
    F = (xiyn * two_pi * sqrt(sigyyn) * (sqrt(sigyyn) + sqrt(sigxxn))) /
        (lf.betay * double(nsli));
}


void BeamBeam3D::computeSlices() {
    static const double border = 8.0;
    double bord1 = - border;
    double bord2;
    double sigz = lf.sigt / cphi;

    for(int i = 1; i <= nsli; ++i) {
        if(i != nsli) {
            bord2 = InverseGauss(double(i) / double(nsli));
        } else {
            bord2 = + border;
        }

        BeamBeam3D::Slice slice;
        slice.zstar = (exp(-sqr(bord2) / 2.0) - exp(-sqr(bord1) / 2.0)) /
                      sqrt(two_pi) * double(nsli);
        bord1 = bord2;

        slice.zstar = displacement.getZ() + slice.zstar * sigz;
        slice.xstar = displacement.getX() + slice.zstar * sphi;
        slice.ystar = displacement.getY();
        slice.sigx  =  lf.emitx * lf.betax + sqr(lf.etax * lf.sige);
        slice.sigpx = (lf.emitx / lf.betax + sqr(lf.etapx * lf.sige)) / sqr(cphi);
        slice.sigy  =  lf.emity * lf.betay + sqr(lf.etay * lf.sige);
        slice.sigpy = (lf.emity / lf.betay + sqr(lf.etapy * lf.sige)) / sqr(cphi);
        slices.push_back(slice);
    }
}


void BeamBeam3D::synchroBeamCollision(PartBunch &bunch) const {
    std::vector<Slice>::const_iterator last_slice = slices.end();
    std::vector<Slice>::const_iterator slice = slices.begin();


    for(; slice != last_slice; ++slice) {

        for(unsigned int i = 0; i < bunch.getLocalNum(); i++) {
            Particle part = bunch.get_part(i);
            double s  = (part.t() - slice->zstar) / 2.0;
            double sx = slice->sigx + slice->sigpx * s * s;
            double sy = slice->sigy + slice->sigpy * s * s;

            double sepx = part.x() + part.px() * s - slice->xstar;
            double sepy = part.y() + part.py() * s - slice->ystar;
            double bbfx, bbfy, bbgx, bbgy;

            if(std::abs(sx - sy) < 1.0e-6 * std::abs(sx + sy)) {
                double x = sepx * sepx + sepy * sepy;
                if(x != 0.0) {
                    double factor = x / (sx + sy);
                    double expfac = exp(-factor);
                    bbfx = 2.0 * sepx * (1.0 - expfac) / x;
                    bbfy = 2.0 * sepy * (1.0 - expfac) / x;
                    double comfac = -sepx * bbfx + sepy * bbfy;
                    bbgx = (+comfac + 4.0 * sepx * sepx * factor * expfac / x) / (2.0 * x);
                    bbgy = (-comfac + 4.0 * sepy * sepy * factor * expfac / x) / (2.0 * x);
                }
            } else if(sx > sy) {
                bbf(sepx, sepy, sx, sy, bbfx, bbfy, bbgx, bbgy);
            } else {
                bbf(sepy, sepx, sy, sx, bbfy, bbfx, bbgy, bbgx);
            }

            bbfx *= F;
            bbfy *= F;
            bbgx *= F;
            bbgy *= F;

            part.x()  += s * bbfx;
            part.px() -= bbfx;
            part.y()  += s * bbfy;
            part.py() -= bbfy;
            part.pt() -= s * (slice->sigpx * bbgx + slice->sigpy * bbgy) +
                         (bbfx * (part.px() - bbfx / 2.0) + bbfy * (part.py() - bbfy / 2.0)) / 2.0;
            bunch.set_part(part, i);
        }
    }
}


void BeamBeam3D::synchroBeamCollision(FVps<double, 6> &map) const {
    std::vector<Slice>::const_iterator last_slice = slices.end();
    std::vector<Slice>::const_iterator slice = slices.begin();

    for(; slice != last_slice; ++slice) {
        Series s  = (map[4] - slice->zstar) / 2.0;
        Series sx = slice->sigx + slice->sigpx * s * s;
        Series sy = slice->sigy + slice->sigpy * s * s;

        Series sepx = map[0] + map[1] * s - slice->xstar;
        Series sepy = map[2] + map[3] * s - slice->ystar;
        Series bbfx, bbfy, bbgx, bbgy;

        if(std::abs(sx[0] - sy[0]) < 1.0e-6 * std::abs(sx[0] + sy[0])) {
            Series x = sepx * sepx + sepy * sepy;
            if(x[0] != 0.0) {
                Series factor = x / (sx + sy);
                Series expfac = exp(-factor);
                bbfx = 2.0 * sepx * (1.0 - expfac) / x;
                bbfy = 2.0 * sepy * (1.0 - expfac) / x;
                Series comfac = -sepx * bbfx + sepy * bbfy;
                bbgx = (+comfac + 4.0 * sepx * sepx * factor * expfac / x) / (2.0 * x);
                bbgy = (-comfac + 4.0 * sepy * sepy * factor * expfac / x) / (2.0 * x);
            }
        } else if(sx[0] > sy[0]) {
            bbf(sepx, sepy, sx, sy, bbfx, bbfy, bbgx, bbgy);
        } else {
            bbf(sepy, sepx, sy, sx, bbfy, bbfx, bbgy, bbgx);
        }

        bbfx *= F;
        bbfy *= F;
        bbgx *= F;
        bbgy *= F;

        map[0] += s * bbfx;
        map[1] -= bbfx;
        map[2] += s * bbfy;
        map[3] -= bbfy;
        map[5] -= s * (slice->sigpx * bbgx + slice->sigpy * bbgy) +
                  (bbfx * (map[1] - bbfx / 2.0) + bbfy * (map[3] - bbfy / 2.0)) / 2.0;
    }
}


void BeamBeam3D::bbf(double sepx,  double sepy,  double sigxx, double sigyy,
                     double &bbfx, double &bbfy, double &bbgx, double &bbgy)
const {
    static const double sqrpi2 = 2.0 * sqrt(pi);

    double x = sepx * sepx / sigxx + sepy * sepy / sigyy;
    double fac2 = 2.0 * std::abs(sigxx - sigyy);
    double fac  = sqrt(fac2);
    double sigxy = sqrt(sigxx / sigyy);
    double expfac = 0.0;
    complex<double> w1 =
        (*errorFunction)(complex<double>(std::abs(sepx), std::abs(sepy)) / fac);

    if(x < 100.0) {
        expfac = exp(- x * 0.5);
        w1 -= expfac * (*errorFunction)
              (complex<double>(std::abs(sepx) / sigxy, std::abs(sepy) * sigxy) / fac);
    }

    complex<double> bbf = (sqrpi2 / fac) * w1;
    bbfx = std::imag(bbf);
    bbfy = std::real(bbf);
    if(sepx < 0) bbfx = - bbfx;
    if(sepy < 0) bbfy = - bbfy;
    double comfac = sepx * bbfx + sepy * bbfy;
    bbgx = - (comfac + 2.0 * (expfac / sigxy - 1.0)) / fac2;
    bbgy = (comfac + 2.0 * (expfac * sigxy - 1.0)) / fac2;
}


void BeamBeam3D::bbf(const Series &sepx, const Series &sepy,
                     const Series &sigxx, const Series &sigyy,
                     Series &bbfx, Series &bbfy, Series &bbgx, Series &bbgy)
const {
    static const double sqrpi2 = 2.0 * sqrt(pi);

    Series x = sepx * sepx / sigxx + sepy * sepy / sigyy;
    Series fac2 = 2.0 * ((sigxx[0] > sigyy[0]) ? sigxx - sigyy : sigyy - sigxx);
    Series fac  = sqrt(fac2);
    Series sigxy = sqrt(sigxx / sigyy);
    Series expfac = 0.0;
    Series argx = (sepx[0] < 0.0) ? -sepx : sepx;
    Series argy = (sepy[0] < 0.0) ? -sepy : sepy;
    Series wx;
    Series wy;
    TpsWerrf(argx / fac, argy / fac, wx, wy);

    if(x[0] < 100.0) {
        expfac = exp(- x * 0.5);
        Series wx2;
        Series wy2;
        TpsWerrf(argx / (sigxy * fac), (argy * sigxy) / fac, wx2, wy2);
        wx -= expfac * wx2;
        wy -= expfac * wy2;
    }

    bbfx = (sqrpi2 / fac) * wy;
    bbfy = (sqrpi2 / fac) * wx;

    if(sepx[0] < 0) bbfx = - bbfx;
    if(sepy[0] < 0) bbfy = - bbfy;
    Series comfac = sepx * bbfx + sepy * bbfy;
    bbgx = - (comfac + 2.0 * (expfac / sigxy - 1.0)) / fac2;
    bbgy = (comfac + 2.0 * (expfac * sigxy - 1.0)) / fac2;
}
