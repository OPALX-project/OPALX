// ------------------------------------------------------------------------
// $RCSfile: BeamStripping.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamStripping
//   Defines the abstract interface for a beam collimator for cyclotrons.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/BeamStripping.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Solvers/BeamStrippingPhysics.hh"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Structure/LossDataSink.h"
#include "Utilities/LogicalError.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <memory>
#include <fstream>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Ippl.h"

using Physics::m_p;
using Physics::m_e;
using Physics::r_e;
using Physics::z_p;
using Physics::kB;
using Physics::q_e;

extern Inform *gmsg;

using namespace std;

// Class BeamStripping
// ------------------------------------------------------------------------

BeamStripping::BeamStripping():
    Component(),
	filename_m(""),
    informed_m(false),
	pressure_m(0.0),
	temperature_m(0.0),
	sigma_m(0.0),
	energycs_m(0.0),
//	Eng(0.0),
//	CS(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    parmatintbst_m(NULL)
{}

BeamStripping::BeamStripping(const BeamStripping &right):
    Component(right),
    filename_m(right.filename_m),
    informed_m(right.informed_m),
	pressure_m(right.pressure_m),
	temperature_m(right.temperature_m),
	sigma_m(right.sigma_m),
	energycs_m(right.energycs_m),
//	Eng(right.Eng),
//	CS(right.CS),
    losses_m(0),
    lossDs_m(nullptr),
    parmatintbst_m(NULL)
{
    setGeom();
}

BeamStripping::BeamStripping(const std::string &name):
    Component(name),
    filename_m(""),
    informed_m(false),
    pressure_m(0.0),
	temperature_m(0.0),
	sigma_m(0.0),
	energycs_m(0.0),
//	Eng(0.0),
//	CS(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    parmatintbst_m(NULL)
{}


BeamStripping::~BeamStripping() {
    if (online_m)
        goOffline();
}


void BeamStripping::accept(BeamlineVisitor &visitor) const {
    visitor.visitBeamStripping(*this);
}



void BeamStripping::setCrossSection(vector<double> sigma) {
    sigma_m = sigma;
}
vector<double> BeamStripping::getCrossSection() const {
  return sigma_m;
}

void BeamStripping::setEnergyCS(vector<double> energycs) {
    energycs_m = energycs;
}
vector<double> BeamStripping::getEnergyCS() const {
  return energycs_m;
}

//double BeamStripping::CrossSection(){
//
//	double CS;
//	double Eng;
//	double m;
//	double n;
//	int a = energycs_m.size();
//	Eng=60000;
//	*gmsg << "* Energy = " << Eng << endl;
//
//	if(Eng < energycs_m[0]) {
//		m = (sigma_m[1]-sigma_m[0]) / (energycs_m[1]-energycs_m[0]);
//		n = sigma_m[0] - m * energycs_m[0];
//		CS = m * Eng + n;
//	}
//	else if(Eng > energycs_m[a-1]) {
//		m = (sigma_m[a-1]-sigma_m[a-2]) / (energycs_m[a-1]-energycs_m[a-2]);
//		n = sigma_m[a-1] - m * energycs_m[a-1];
//		CS = m * Eng + n;
//	}
//	else if(Eng > energycs_m[0] && Eng < energycs_m[a-1]){
//		for (int i=0; i<a; i++){
//			if(Eng == energycs_m[i]) {
//				CS = sigma_m[i];
//			}
//			else if(Eng > energycs_m[i] && Eng < energycs_m[i+1]) {
//				m = (sigma_m[i+1]-sigma_m[i]) / (energycs_m[i+1]-energycs_m[i]);
//				n = sigma_m[i] - m * energycs_m[i];
//				CS = m * Eng + n;
//			}
//		}
//	}
//	else {
//		*gmsg << "**Cross section not calculated**" << endl;
//	}
//	*gmsg << "* Cross section = " << CS << endl;
//	return CS;
//}



bool BeamStripping::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool BeamStripping::applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool BeamStripping::checkBeamStripping(Vector_t r, Vector_t rmin, Vector_t rmax) {

//    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
//    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
//    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    bool isDead = false;
//    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
//        if ( r1 > r_start - 10.0 && r1 < r_end + 10.0 ){
//            if (r(2) < zend_m && r(2) > zstart_m ) {
//                int pflag = checkPoint(r(0), r(1));
//                isDead = (pflag != 0);
//            }
//        }
//    }
    return isDead;
}


// rectangle collimators in cyclotron cyclindral coordinates
// without particlematterinteraction, the particle hitting collimator is deleted directly
bool BeamStripping::checkBeamStripping(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = true;
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);
//    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
//    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
//    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    std::pair<Vector_t, double> boundingSphere;
    boundingSphere.first = 0.5 * (rmax + rmin);
    boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);

//    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
//        // if ( r1 > r_start - 10.0 && r1 < r_end + 10.0 ){
//        if ( r1 > r_start - 100.0 && r1 < r_end + 100.0 ){
//            size_t tempnum = bunch->getLocalNum();
//            int pflag = 0;
//            for (unsigned int i = 0; i < tempnum; ++i) {
//                if (bunch->PType[i] == ParticleType::REGULAR && bunch->R[i](2) < zend_m && bunch->R[i](2) > zstart_m ) {
//                    pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
//		    /// bunch->Bin[i] != -1 makes sure the partcile is not stored in more than one collimator
//                    if ((pflag != 0) && (bunch->Bin[i] != -1))  {
//		      if (!parmatintbst_m)
//			lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i]);
//		      bunch->Bin[i] = -1;
//		      flagNeedUpdate = true;
//                    }
//                }
//            }
//        }
//    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if (flagNeedUpdate && parmatintbst_m) {
        parmatintbst_m->apply(bunch, boundingSphere);
    }
    return flagNeedUpdate;
}

void BeamStripping::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    initialise(bunch);
}

void BeamStripping::initialise(PartBunchBase<double, 3> *bunch) {
    RefPartBunch_m = bunch;

    parmatintbst_m = getParticleMatterInteraction();

    // if (!parmatintbst_m) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
    // }
    goOnline(-1e6);
}

void BeamStripping::finalise() {
    if (online_m)
        goOffline();
    *gmsg << "* Finalize Beam Stripping" << endl;
}

void BeamStripping::goOnline(const double &) {
//    CrossSection();
    print();
//    // PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // PosZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    // time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
////    id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//    online_m = true;
}

void BeamStripping::print() {
    if (RefPartBunch_m == NULL) {
        if (!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            ERRORMSG(errormsg << endl);
            if (Ippl::myNode() == 0) {
                std::ofstream omsg("errormsg.txt", std::ios_base::app);
                omsg << errormsg << std::endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }
//    int a = sigma_m.size();
//    for (int i=0; i<a; i++){
//    	*gmsg << "Energy - cross section  " << energycs_m[i] << "  -  " << sigma_m[i] << endl;
//    }
}

void BeamStripping::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    lossDs_m.reset(0);
    online_m = false;
}

bool BeamStripping::bends() const {
    return false;
}

void BeamStripping::setOutputFN(std::string fn) {
    filename_m = fn;
}

std::string BeamStripping::getOutputFN() {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

void BeamStripping::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementBase::ElementType BeamStripping::getType() const {
    return BEAMSTRIPPING;
}

std::string BeamStripping::getBeamStrippingShape() {
    return "BeamStripping";
}

void BeamStripping::setGeom() {
//    double slope;
//    if (xend_m == xstart_m)
//        slope = 1.0e12;
//    else
//        slope = (yend_m - ystart_m) / (xend_m - xstart_m);
//
//    double coeff2 = sqrt(1 + slope * slope);
//    double coeff1 = slope / coeff2;
//    double halfdist = width_m / 2.0;
//    geom_m[0].x = xstart_m - halfdist * coeff1;
//    geom_m[0].y = ystart_m + halfdist / coeff2;
//
//    geom_m[1].x = xstart_m + halfdist * coeff1;
//    geom_m[1].y = ystart_m - halfdist / coeff2;
//
//    geom_m[2].x = xend_m + halfdist * coeff1;
//    geom_m[2].y = yend_m - halfdist / coeff2;
//
//    geom_m[3].x = xend_m - halfdist * coeff1;
//    geom_m[3].y = yend_m + halfdist / coeff2;
//
//    geom_m[4].x = geom_m[0].x;
//    geom_m[4].y = geom_m[0].y;
//
//    if (zstart_m > zend_m){
//        double tempz = zstart_m;
//        zstart_m = zend_m;
//        zend_m = tempz;
//    }
}

int BeamStripping::checkPoint() {
    int cn = 1;
//    for (int i = 0; i < 4; i++) {
//        if (((geom_m[i].y <= y) && (geom_m[i+1].y > y))
//            || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {
//
//            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
//            if (x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
//                ++cn;
//        }
//    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}
