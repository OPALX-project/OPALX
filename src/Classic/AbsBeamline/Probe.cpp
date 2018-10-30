// ------------------------------------------------------------------------
// $RCSfile: Probe.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Probe
//   Defines the abstract interface for a septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 09:32:32 $
// $Author: Bi, Yang $
// 2012/03/01: fix bugs and change the algorithm in the checkProbe()
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Probe.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Structure/PeakFinder.h"
#include "Utilities/Options.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class Probe
// ------------------------------------------------------------------------

Probe::Probe():
    Component(),
    filename_m(""),
    position_m(0.0),
    step_m(0.0){
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0);
}


Probe::Probe(const Probe &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    step_m(right.step_m){
    setDimensions(right.xstart_m, right.xend_m, right.ystart_m, right.yend_m, right.width_m);
}


Probe::Probe(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    step_m(0.0){
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0);
}

Probe::~Probe() {}

void Probe::accept(BeamlineVisitor &visitor) const {
    visitor.visitProbe(*this);
}

void Probe::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    position_m = startField;
    startField -= 0.005;
    endField = position_m + 0.005;
}

void Probe::initialise(PartBunchBase<double, 3> *bunch) {
    *gmsg << "* Initialize probe" << endl;

    std::string name;
    if (filename_m == std::string(""))
        name = getName();
    else
        name = filename_m.substr(0, filename_m.rfind("."));

    bool singlemode = (bunch->getTotalNum() == 1) ? true : false;
    peakfinder_m = std::unique_ptr<PeakFinder>  (new PeakFinder  (name, rstart_m, rend_m, step_m, singlemode));
    lossDs_m     = std::unique_ptr<LossDataSink>(new LossDataSink(name, !Options::asciidump));
}

void Probe::finalise() {
    *gmsg << "* Finalize probe " << getName() << endl;
    peakfinder_m->save();
    lossDs_m->save();
}

bool Probe::bends() const {
    return false;
}

void Probe::goOffline() {
    *gmsg << "* Probe goes offline " << getName() << endl;
    online_m = false;
    peakfinder_m->save();
    lossDs_m->save();
}

void Probe::setDimensions(double xstart, double xend, double ystart, double yend, double width) {
    xstart_m = xstart;
    ystart_m = ystart;
    xend_m   = xend;
    yend_m   = yend;
    width_m  = width;
    rstart_m = std::sqrt(xstart*xstart + ystart * ystart);
    rend_m   = std::sqrt(xend * xend   + yend * yend);
    // start position is the one with lowest radius
    if (rstart_m > rend_m) {
        std::swap(xstart_m, xend_m);
        std::swap(ystart_m, yend_m);
        std::swap(rstart_m, rend_m);
    }

    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = std::sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}

void  Probe::setStep(double step) {
    step_m = step;
}

double  Probe::getXstart() const {
    return xstart_m;
}

double  Probe::getXend() const {
    return xend_m;
}

double  Probe::getYstart() const {
    return ystart_m;
}

double  Probe::getYend() const {
    return yend_m;
}

double  Probe::getWidth() const {
    return width_m;
}

double  Probe::getStep() const {
    return step_m;
}

void Probe::setGeom(const double dist) {

    double slope;
    if (xend_m == xstart_m)
      slope = 1.0e12;
    else
      slope = (yend_m - ystart_m) / (xend_m - xstart_m);

    double coeff2 = sqrt(1 + slope * slope);
    double coeff1 = slope / coeff2;
    double halfdist = dist / 2.0;
    geom_m[0].x = xstart_m - halfdist * coeff1;
    geom_m[0].y = ystart_m + halfdist / coeff2;

    geom_m[1].x = xstart_m + halfdist * coeff1;
    geom_m[1].y = ystart_m - halfdist / coeff2;

    geom_m[2].x = xend_m + halfdist * coeff1;
    geom_m[2].y = yend_m - halfdist  / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;
}

bool  Probe::checkProbe(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    bool flagprobed = false;
    Vector_t rmin, rmax, probepoint;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute minimum and maximum
    double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_min = std::hypot(xmin, ymin);
    double rbunch_max = std::hypot(xmax, ymax);

    if( rbunch_max > rstart_m - 10.0 && rbunch_min < rend_m + 10.0 ) {
        size_t tempnum = bunch->getLocalNum();
        int pflag = 0;

        Vector_t meanP(0.0, 0.0, 0.0);
        for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
            for(int d = 0; d < 3; ++d) {
                meanP(d) += bunch->P[i](d);
            }
        }
        reduce(meanP, meanP, OpAddAssign());
        meanP = meanP / Vector_t(bunch->getTotalNum());

        double sk1, sk2, stangle = 0.0;
        if ( B_m == 0.0 ){
            sk1 = meanP(1)/meanP(0);
            if (sk1 == 0.0)
                stangle = 1.0e12;
            else
                stangle = std::abs(1/sk1);
        } else if (meanP(0) == 0.0 ) {
            sk2 = - A_m/B_m;
            if ( sk2 == 0.0 )
                stangle = 1.0e12;
            else
                stangle = std::abs(1/sk2);
        } else {
            sk1 = meanP(1)/meanP(0);
            sk2 = - A_m/B_m;
            stangle = std::abs(( sk1-sk2 )/(1 + sk1*sk2));
        }
        // change probe width depending on step size on angle of particle
        double lstep = (sqrt(1.0-1.0/(1.0+dot(meanP, meanP))) * Physics::c) * tstep*1.0e-6; // [mm]
        double Swidth = lstep / sqrt( 1 + 1/stangle/stangle );
        setGeom(Swidth);

        for(unsigned int i = 0; i < tempnum; ++i) {
            pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
            if(pflag != 0) {
                // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
                double dist1 = (A_m*bunch->R[i](0)+B_m*bunch->R[i](1)+C_m)/R_m/1000.0;
                double k1, k2, tangle = 0.0;
                if ( B_m == 0.0 ){
                    k1 = bunch->P[i](1)/bunch->P[i](0);
                    if (k1 == 0.0)
                        tangle = 1.0e12;
                    else
                        tangle = std::abs(1/k1);
                } else if ( bunch->P[i](0) == 0.0 ) {
                        k2 = -A_m/B_m;
                        if (k2 == 0.0)
                            tangle = 1.0e12;
                        else
                            tangle = std::abs(1/k2);
                } else {
                    k1 = bunch->P[i](1)/bunch->P[i](0);
                    k2 = -A_m/B_m;
                    tangle = std::abs(( k1-k2 )/(1 + k1*k2));
                }
                double dist2 = dist1 * sqrt( 1+1/tangle/tangle );
                double dt = dist2/(sqrt(1.0-1.0/(1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c)*1.0e9;

                probepoint(0) = (B_m*B_m*bunch->R[i](0) - A_m*B_m*bunch->R[i](1)-A_m*C_m)/(R_m*R_m);
                probepoint(1) = (A_m*A_m*bunch->R[i](1) - A_m*B_m*bunch->R[i](0)-B_m*C_m)/(R_m*R_m);
                probepoint(2) = bunch->R[i](2);
                lossDs_m->addParticle(probepoint, bunch->P[i], bunch->ID[i], t+dt, turnnumber);
                flagprobed = true;
            
                peakfinder_m->addParticle(bunch->R[i]);
            }
        }
    }

    reduce(&flagprobed, &flagprobed + 1, &flagprobed, OpBitwiseOrAssign());
    return flagprobed;
}

int Probe::checkPoint(const double &x, const double &y) {
    int    cn = 0;

    for(int i = 0; i < 4; i++) {
        if(((geom_m[i].y <= y) && (geom_m[i+1].y > y))
           || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if(x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)

}

void Probe::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd   = position_m + 0.005;
}

ElementBase::ElementType Probe::getType() const {
    return PROBE;
}