// ------------------------------------------------------------------------
// $RCSfile: Stripper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Stripper
//   Defines the abstract interface for a stripper
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2011/07/08 11:16:04 $
// $Author: Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include <iostream>
#include <fstream>


using Physics::pi;
using Physics::q_e;

extern Inform *gmsg;

using namespace std;

// Class Stripper
// ------------------------------------------------------------------------

Stripper::Stripper():
    Component(),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    opcharge_m(0.0),
    opmass_m(0.0),
    opyield_m(1.0),
    stop_m(true),
    step_m(0) {
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}


Stripper::Stripper(const Stripper &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    width_m(right.width_m),
    opcharge_m(right.opcharge_m),
    opmass_m(right.opmass_m),
    opyield_m(right.opyield_m),
    stop_m(right.stop_m),
    step_m(right.step_m) {
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}


Stripper::Stripper(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    opcharge_m(0.0),
    opmass_m(0.0),
    opyield_m(1.0),
    stop_m(true),
    step_m(0){
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}

void Stripper::setGeom(const double dist) {

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


Stripper::~Stripper() {
    idrec_m.clear();
}


void Stripper::accept(BeamlineVisitor &visitor) const {
    visitor.visitStripper(*this);
}

void Stripper::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}

void Stripper::initialise(PartBunchBase<double, 3> *bunch) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}

void Stripper::finalise() {
    *gmsg << "* Finalize Stripper " << getName() << endl;
    lossDs_m->save();
}

bool Stripper::bends() const {
    return false;
}

void Stripper::goOffline() {
    online_m = false;
    lossDs_m->save();
}

void  Stripper::setXstart(double xstart) {
    xstart_m = xstart;
}

void  Stripper::setXend(double xend) {
    xend_m = xend;
}

void  Stripper::setYstart(double ystart) {
    ystart_m = ystart;
}

void  Stripper::setYend(double yend) {
    yend_m = yend;
}
void  Stripper::setWidth(double width) {
    width_m = width;
}

void  Stripper::setOPCharge(double charge) {
    opcharge_m = charge;
}

void  Stripper::setOPMass(double mass) {
    opmass_m = mass;
}

void  Stripper::setOPYield(double yield) {
    opyield_m = yield;
}

void  Stripper::setStop(bool stopflag) {
    stop_m = stopflag;

}

double  Stripper::getXstart() const {
    return xstart_m;
}

double  Stripper::getXend() const {
    return xend_m;
}

double  Stripper::getYstart() const {
    return ystart_m;
}

double  Stripper::getYend() const {
    return yend_m;
}
double  Stripper::getWidth() const {
    return width_m;
}

double  Stripper::getOPCharge() const {
    return opcharge_m;
}

double  Stripper::getOPMass() const {
    return opmass_m;
}

double  Stripper::getOPYield() const {
    return opyield_m;
}

bool  Stripper::getStop () const {
    return stop_m;
}


//change the stripped particles to outcome particles
bool  Stripper::checkStripper(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    bool flagresetMQ = false;
    Vector_t rmin, rmax, strippoint;
    bunch->get_bounds(rmin, rmax);
    double r_ref = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));

    if(r1 > r_ref - 10.0 ){

        size_t count = 0;
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
            if(sk1 == 0.0)
                stangle =1.0e12;
            else
                stangle = abs(1/sk1);
        }else if (meanP(0) == 0.0 ){
            sk2 = - A_m/B_m;
            if(sk2 == 0.0)
              stangle =1.0e12;
            else
              stangle = abs(1/sk2);
        }else {
            sk1 = meanP(1)/meanP(0);
            sk2 = - A_m/B_m;
            stangle = abs(( sk1-sk2 )/(1 + sk1*sk2));
        }
        double lstep = (sqrt(1.0-1.0/(1.0+dot(meanP, meanP))) * Physics::c) * tstep*1.0e-6; // [mm]
        double Swidth = lstep /  sqrt( 1+1/stangle/stangle ) * 1.2;
        setGeom(Swidth);

        for(unsigned int i = 0; i < tempnum; ++i) {
            if(bunch->PType[i] == ParticleType::REGULAR || bunch->PType[i] == ParticleType::STRIPPED) {
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
                            tangle = abs(1/k1);
                    }else if (bunch->P[i](0) == 0.0 ){
                        k2 = -A_m/B_m;
                        if (k2 == 0.0)
                            tangle = 1.0e12;
                        else
                            tangle = abs(1/k2);
                    }else {
                        k1 = bunch->P[i](1)/bunch->P[i](0);
                        k2 = -A_m/B_m;
                        tangle = abs(( k1-k2 )/(1 + k1*k2));
                    }
                    double dist2 = dist1 * sqrt( 1+1/tangle/tangle );
                    double dt = dist2/(sqrt(1.0-1.0/(1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c)*1.0e9;
                    strippoint(0) = (B_m*B_m*bunch->R[i](0) - A_m*B_m*bunch->R[i](1)-A_m*C_m)/(R_m*R_m);
                    strippoint(1) = (A_m*A_m*bunch->R[i](1) - A_m*B_m*bunch->R[i](0)-B_m*C_m)/(R_m*R_m);
                    strippoint(2) = bunch->R[i](2);
                    lossDs_m->addParticle(strippoint, bunch->P[i], bunch->ID[i], t+dt, turnnumber);

                    if (stop_m) {
                        bunch->Bin[i] = -1;
                        flagNeedUpdate = true;
                        *gmsg << "* Total number of particles lost in stripper = " << bunch->getTotalNum() << endl;
                    }
                    else {
                    	bunch->updateNumTotal();
                    	*gmsg << "* Total number of particles after stripper = " << bunch->getTotalNum() << endl;

                        flagNeedUpdate = true;
                        // change charge and mass of PartData when the reference particle hits the stripper.
                        if(bunch->ID[i] == 0)
                          flagresetMQ = true;

                        // change the mass and charge
                        bunch->M[i] = opmass_m;
                        bunch->Q[i] = opcharge_m * q_e;
                        bunch->PType[i] = ParticleType::STRIPPED;

                        int j = 1;
                        //create new particles
                        while (j < opyield_m){
                          bunch->create(1);
                          bunch->R[tempnum+count] = bunch->R[i];
                          bunch->P[tempnum+count] = bunch->P[i];
                          bunch->Q[tempnum+count] = bunch->Q[i];
                          bunch->M[tempnum+count] = bunch->M[i];
                          // once the particle is stripped, change PType from 0 to 1 as a flag so as to avoid repetitive stripping.
                          bunch->PType[tempnum+count] = ParticleType::STRIPPED;
                          count++;
                          j++;
                        }

                        if(bunch->weHaveBins())
                            bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];

                    }
                }
            }
        }
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

    if(!stop_m){
        reduce(&flagresetMQ, &flagresetMQ + 1, &flagresetMQ, OpBitwiseOrAssign());
        if(flagresetMQ){
            bunch->resetM(opmass_m * 1.0e9); // GeV -> eV
            bunch->resetQ(opcharge_m);     // elementary charge
        }
    }

    return flagNeedUpdate;
}


void Stripper::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}

ElementBase::ElementType Stripper::getType() const {
    return STRIPPER;
}


int Stripper::checkPoint(const double &x, const double &y) {
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
