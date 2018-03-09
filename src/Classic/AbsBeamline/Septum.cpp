// ------------------------------------------------------------------------
// $RCSfile: Septum.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Septum
//   Defines the abstract interface for a septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Septum.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h" // OPAL file
#include <iostream>
#include <fstream>
using Physics::pi;

extern Inform *gmsg;

// Class Septum
// ------------------------------------------------------------------------

Septum::Septum():
    Component(),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0)
{}


Septum::Septum(const Septum &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    width_m(right.width_m)
{}


Septum::Septum(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0)
{}


Septum::~Septum() {
}


void Septum::accept(BeamlineVisitor &visitor) const {
    visitor.visitSeptum(*this);
}

void Septum::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    position_m = startField;
    startField -= 0.005;
    endField = position_m + 0.005;


}

void Septum::initialise(PartBunchBase<double, 3> *bunch) {
    *gmsg << "Septum initialise" << endl;
}

void Septum::finalise()
{}

bool Septum::bends() const {
    return false;
}

void Septum::goOffline() {
    online_m = false;
}

void  Septum::setXstart(double xstart) {
    xstart_m = xstart;
}

void  Septum::setXend(double xend) {
    xend_m = xend;
}

void  Septum::setYstart(double ystart) {
    ystart_m = ystart;
}


void  Septum::setYend(double yend) {
    yend_m = yend;
}
void  Septum::setWidth(double width) {
    width_m = width;
}


double  Septum::getXstart() const {
    return xstart_m;

}

double  Septum::getXend() const {
    return xend_m;

}

double  Septum::getYstart() const {
    return ystart_m;

}

double  Septum::getYend() const {
    return yend_m;

}
double  Septum::getWidth() const {
    return width_m;

}

bool  Septum::checkSeptum(PartBunchBase<double, 3> *bunch) {

    bool flag = false;
    Vector_t rmin;
    Vector_t rmax;
    bunch->get_bounds(rmin, rmax);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    if(r1 > sqrt(xstart_m * xstart_m + ystart_m * ystart_m) - 100)  {
        for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
            Vector_t R = bunch->R[i];
            double slope = (yend_m - ystart_m) / (xend_m - xstart_m);
            double intcept = ystart_m - slope * xstart_m;
            double intcept1 = intcept - width_m / 2.0 * sqrt(slope * slope + 1);
            double intcept2 = intcept + width_m / 2.0 * sqrt(slope * slope + 1);

            double line1 = fabs(slope * R(0) + intcept1);
            double line2 = fabs(slope * R(0) + intcept2);



            if(fabs(R(1)) > line2 && fabs(R(1)) < line1 && R(0) > xstart_m && R(0) < xend_m && R(1) > ystart_m && R(1) < yend_m) {

                bunch->lossDs_m->addParticle(R, bunch->P[i], bunch->ID[i]);
                bunch->Bin[i] = -1;
                flag = true;

            }

        }
    }
    return flag;
}



// angle range [0~2PI) degree
double Septum::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    // if(x < 0)                   thetaXY = pi + atan(y / x);
    // else if((x > 0) && (y >= 0))  thetaXY = atan(y / x);
    // else if((x > 0) && (y < 0))   thetaXY = 2.0 * pi + atan(y / x);
    // else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    // else if((x == 0) && (y < 0)) thetaXY = 3.0 / 2.0 * pi;

    return thetaXY >= 0 ? thetaXY: thetaXY + Physics::two_pi;

}
void Septum::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}


ElementBase::ElementType Septum::getType() const {
    return SEPTUM;
}