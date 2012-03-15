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
// $Author: bi $

// ------------------------------------------------------------------------

#include "AbsBeamline/Probe.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include <iostream>
#include <fstream>
using Physics::pi;

extern Inform *gmsg;

// Class Probe
// ------------------------------------------------------------------------

Probe::Probe():
    Component(),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    step_m(0)
{}


Probe::Probe(const Probe &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    width_m(right.width_m),
    step_m(right.step_m)
{}


Probe::Probe(const string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    step_m(0)
{}


Probe::~Probe() {
    idrec_m.clear();
}


void Probe::accept(BeamlineVisitor &visitor) const {
    visitor.visitProbe(*this);
}

bool Probe::apply(const int &i, const double &t, double E[], double B[]) {
    *gmsg << "septum1" << endl;
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Probe::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    *gmsg << "septum2" << endl;
    return false;
}

bool Probe::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    *gmsg << "septum3" << endl;
    return false;
}

void Probe::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    position_m = startField;
    startField -= 0.005;
    endField = position_m + 0.005;
}

void Probe::initialise(PartBunch *bunch, const double &scaleFactor) {
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(1000000, doH5);
    lossDs_m->openH5(getName());
}

void Probe::finalise() {
    *gmsg << "Finalize probe" << endl;
    if(lossDs_m)
        delete lossDs_m;
}

bool Probe::bends() const {
    return false;
}
void Probe::goOnline() {
    online_m = true;
}

void Probe::goOffline() {
    online_m = false;
}

void  Probe::setXstart(double xstart) {
    xstart_m = xstart;
}

void  Probe::setXend(double xend) {
    xend_m = xend;
}

void  Probe::setYstart(double ystart) {
    ystart_m = ystart;
}


void  Probe::setYend(double yend) {
    yend_m = yend;
}
void  Probe::setWidth(double width) {
    width_m = width;
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


bool  Probe::checkProbe(PartBunch &bunch, int turnnumber, double deltaTheta) {

    bool flag = false;
    Vector_t rmin;
    Vector_t rmax;
    bunch.get_bounds(rmin, rmax);
    double slope = (yend_m - ystart_m) / (xend_m - xstart_m);
    double intcept = ystart_m - slope * xstart_m;
    double r_ref = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);

    double line_ref = slope * rmax(0) + intcept;
    double rang = deltaTheta * sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1)) * 5;
    double line_ref1 = slope * rmax(0) + intcept - rang;
    double line_ref2 = slope * rmax(0) + intcept + rang;
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    double  slope_cen, tangle;

    if(r1 > r_ref - 100) {
        if(line_ref1 < rmax(1) && line_ref2 > rmax(1)) {
            Vector_t meanR = Vector_t(0.0, 0.0, 0.0);
            Vector_t meanP = Vector_t(0.0, 0.0, 0.0);

            for(int ii = 0; ii < (bunch.getLocalNum()); ii++) {
                for(int j = 0; j < 3; j++) {
                    meanR(j) += bunch.R[ii](j);
                    meanP(j) += bunch.P[ii](j);
                }
            }
            reduce(meanR, meanR, OpAddAssign());
            meanR /= Vector_t((bunch.getTotalNum()));
            reduce(meanP, meanP, OpAddAssign());
            meanP /= Vector_t((bunch.getTotalNum()));

            slope_cen = meanP(1) / meanP(0);
            tangle = abs(slope_cen + 1.0 / slope) / abs(1 - 1.0 / slope * slope_cen);
            double stepsize = deltaTheta * sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1)) * 2;
            width_m = stepsize / sqrt(1 + tangle * tangle);
        }
        double intcept1 = intcept - width_m / 2.0 * sqrt(slope * slope + 1);
        double intcept2 = intcept + width_m / 2.0 * sqrt(slope * slope + 1);

        for(int i = 0; i < bunch.getLocalNum(); ++i) {
            Vector_t R = bunch.R[i];
            Vector_t P = bunch.P[i];

            double line1 = slope * R(0) + intcept1;
            double line2 = slope * R(0) + intcept2;

            if(line1 <= R(1) && line2 >= R(1) && R(0) > xstart_m) {
                // trace back to the position of the middle line of the probe,
                // and get the beam parameters on this line.
                double slope_par = P(1) / P(0);
                double intcept_par = R(1) - slope_par * R(0);
                double x_prob = (intcept_par - intcept) / (slope - slope_par);
                double y_prob = slope * x_prob + intcept;
                Vector_t prob = Vector_t(x_prob, y_prob, R(2));
                lossDs_m->addParticle(prob, turnnumber, bunch.ID[i]);

            }
        }
    }
    lossDs_m->save(getName());
    return flag;
}


// angle range [0~2PI) degree
double Probe::calculateAngle(double x, double y) {
    double thetaXY;

    if(x < 0)                   thetaXY = pi + atan(y / x);
    else if((x > 0) && (y >= 0))  thetaXY = atan(y / x);
    else if((x > 0) && (y < 0))   thetaXY = 2.0 * pi + atan(y / x);
    else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    else if((x == 0) && (y < 0)) thetaXY = 3.0 / 2.0 * pi;

    return thetaXY;

}
void Probe::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}


const string &Probe::getType() const {
    static const string type("Probe");
    return type;
}

