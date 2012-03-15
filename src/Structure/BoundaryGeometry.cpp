// ------------------------------------------------------------------------
// $RCSfile: Geometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//   The class for the OPAL GEOMETRY command.
//
// $Date: 2003/08/11 22:09:00 $
// $Author: A. Adelmann $
//
// ------------------------------------------------------------------------

//FIXME: cleanup
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Elements/OpalBeamline.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "hdf5.h"
#include "Distribution/Distribution.h"
#include <sys/stat.h>
#include "Ippl.h"
using namespace Expressions;
using namespace Physics;
using Physics::c;
extern Inform *gmsg;//debug

// Class BoundaryGeometry
// -----------------------------------------------------------------------

// The attributes of class Geometry.
namespace {
    enum {
        FGEOM,    // file holding the geometry
        LENGHT,   // length of elliptic tube or boxcorner
        S,        // start of the geometry
	L1,       // in case of BOXCORNER first part of geometry with hight B
	L2,       // in case of BOXCORNER second part of geometry with hight B-C
        A,        // major semi-axis of elliptic tube
        B,        // minor semi-axis of ellitpic tube
	C,        // in case of BOXCORNER hight of corner
        TOPO,     // BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written
        DISTR,    // Add distribution to generate physics model on the surface
        DISTRS,   // Add distribution array to generate physics model on the surface
	ZSHIFT,   // Shift in z direction
        SIZE
    };

}

/*namespace BGphysics {
    enum TPHYACTION {
        Nop = 0x01,// tringle is transparent to particle like beam window
        Absorption = 0x02,// triangle has no field emission and secondary emission
        FNEmission = 0x04,// triangle has field emission
        SecondaryEmission = 0x08// trangle has secondary emission 
    };
    }*/


BoundaryGeometry::BoundaryGeometry() : Definition(SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry.") {
    itsAttr[FGEOM] = Attributes::makeString
	("FGEOM", "Specifies the geometry file [h5fed]", "");

    itsAttr[TOPO] = Attributes::makeString
	("TOPO", "BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written ", "ELLIPTIC");

    itsAttr[LENGHT] = Attributes::makeReal
	("LENGHT", "Specifies the length of a tube shaped elliptic beam pipe [m]", 1.0);

    itsAttr[S] = Attributes::makeReal
	("S", "Specifies the start of a tube shaped elliptic beam pipe [m]", 0.0);

    itsAttr[A] = Attributes::makeReal
	("A", "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);

    itsAttr[B] = Attributes::makeReal
	("B", "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);

    itsAttr[L1] = Attributes::makeReal
	("L1", "In case of BOXCORNER Specifies first part with hight == B [m]", 0.5);

    itsAttr[L2] = Attributes::makeReal
        ("L2", "In case of BOXCORNER Specifies first second with hight == B-C [m]", 0.2);

    itsAttr[C] = Attributes::makeReal
        ("C", "In case of BOXCORNER Specifies hight of corner C [m]", 0.01);


    itsAttr[DISTR] = Attributes::makeString
	("DISTR", "Distribution to be generated on the surface", "");
    itsAttr[DISTRS] = Attributes::makeStringArray
	("DISTRS", "Distribution array to be generated on the surface");
    
    //itsAttr[ZFLAG] = Attributes::makeReal
    //("ZFLAG", "Flag control the changing of x/y axis with z axis", 0);

    itsAttr[ZSHIFT] = Attributes::makeReal
	("ZSHIFT", "Shift in z direction", 0.0);

    BoundaryGeometry *defGeometry = clone("UNNAMED_GEOMETRY");
    defGeometry->builtin = true;

    TPInside_m = IpplTimings::getTimer("Particle Inside");
    TPreProc_m = IpplTimings::getTimer("Pre Processing");
    TRayTrace_m = IpplTimings::getTimer("Ray tracing");
    h5FileName_m = Attributes::getString(itsAttr[FGEOM]);
    
    try {
        defGeometry->update();
        OPAL.define(defGeometry);
    } catch(...) {
        delete defGeometry;
    }
    if(!h5FileName_m.empty())
        initialize();
}


BoundaryGeometry::BoundaryGeometry(const string &name, BoundaryGeometry *parent):
    Definition(name, parent),
    h5FileName_m(parent->h5FileName_m) {
    h5FileName_m = Attributes::getString(itsAttr[FGEOM]);

   
    if(!h5FileName_m.empty())
        initialize();
    TPInside_m = IpplTimings::getTimer("Particle Inside");
    TPreProc_m = IpplTimings::getTimer("zshiftPre Processing");
    TRayTrace_m = IpplTimings::getTimer("Ray tracing");

}


BoundaryGeometry::~BoundaryGeometry() {
    if(allbfaces_m)
        delete allbfaces_m;
    if(Tribarycent_m)
        delete Tribarycent_m;
    if(TriPrPartlossZ_m)
        delete TriPrPartlossZ_m;
    if(TriSePartlossZ_m)
        delete TriSePartlossZ_m;
    if(TriFEPartlossZ_m)
        delete TriFEPartlossZ_m;
    if(TriPrPartloss_m )
	delete TriPrPartloss_m ;
    if(TriFEPartloss_m)
	delete TriFEPartloss_m;
    if(TriSePartloss_m)
	delete TriSePartloss_m;
}


bool BoundaryGeometry::canReplaceBy(Object *object) {
    // Can replace only by another GEOMETRY.
    return dynamic_cast<Geometry *>(object) != 0;
}


BoundaryGeometry *BoundaryGeometry::clone(const string &name) {
    return new BoundaryGeometry(name, this);
}


void BoundaryGeometry::execute() {
    update();
    TPInside_m = IpplTimings::getTimer("Particle Inside");
    TPreProc_m = IpplTimings::getTimer("Pre Processing");
    TRayTrace_m = IpplTimings::getTimer("Ray tracing");
}


BoundaryGeometry *BoundaryGeometry::find(const string &name) {
    BoundaryGeometry *geom = dynamic_cast<BoundaryGeometry *>(OPAL.find(name));

    if(geom == 0)
        throw OpalException("BoundaryGeometry::find()", "Geometry \"" + name + "\" not found.");
    return geom;
}


string BoundaryGeometry::getFilename() {
    return (string)Attributes::getString(itsAttr[FGEOM]);
}

string BoundaryGeometry::getTopology() {
    return (string)Attributes::getString(itsAttr[TOPO]);
}

double BoundaryGeometry::getA() {
    return (double)Attributes::getReal(itsAttr[A]);
}

double BoundaryGeometry::getB() {
    return (double)Attributes::getReal(itsAttr[B]);
}

double BoundaryGeometry::getC() {
    return (double)Attributes::getReal(itsAttr[C]);
}

double BoundaryGeometry::getS() {
    return (double)Attributes::getReal(itsAttr[S]);
}

double BoundaryGeometry::getLenght() {
    return (double)Attributes::getReal(itsAttr[LENGHT]);
}

double BoundaryGeometry::getL1() {
    return (double)Attributes::getReal(itsAttr[L1]);
}

double BoundaryGeometry::getL2() {
    return (double)Attributes::getReal(itsAttr[L2]);
}
string BoundaryGeometry::getDistribution() {

    return (string)Attributes::getString(itsAttr[DISTR]);
}

vector<string> BoundaryGeometry::getDistributionArray() {

    return Attributes::getStringArray(itsAttr[DISTRS]);
}

/*int BoundaryGeometry::getZflag() {

return (int)(Attributes::getReal(itsAttr[ZFLAG]));
}*/

double BoundaryGeometry::getZshift() {

    return (double)(Attributes::getReal(itsAttr[ZSHIFT]));
}

size_t BoundaryGeometry::getN() { 
	
    return partsr_m.size();
}
void BoundaryGeometry::update() {
    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_GEOMETRY");
}
/**
 * We define some tags in namespace BGphysics for each surface triangle to identify
 * the physical reactions for each triangle when amplitude of electrostatic field
 * exceeds some threshold or particles incident the surface.
 * Caution: we temporarily use this function to identify the particle exiting area
 * of ctf3 gun or cyclotron cavity now. This should be reconsidered later.
 */


void BoundaryGeometry::setBGphysicstag() {
    if( Options::ppdebug ) {
	for(int i = 0; i < numbfaces_global_m; i++) {
	    short temp = 0x00;
	    temp = ((BGphysics::Absorption) | (BGphysics::FNEmission) | (BGphysics::SecondaryEmission));
	    TriBGphysicstag_m.push_back(temp);
	}

    } else {// temp for ciae RF Cavity
	for(int i = 0; i < numbfaces_global_m; i++) {
	    short temp = 0x00;
	    /*if(((Tribarycent_m[i](0) > maxcoords_m[0] - 0.00001) && (Tribarycent_m[i](0) < maxcoords_m[0] + 0.00001))
	      || ((Tribarycent_m[i](1) > maxcoords_m[1] - 0.00001) && (Tribarycent_m[i](1) < maxcoords_m[1] + 0.00001))
	      || ((Tribarycent_m[i](0) > mincoords_m[0] - 0.00001) && (Tribarycent_m[i](0) < mincoords_m[0] + 0.00001))
	      || ((Tribarycent_m[i](1) > mincoords_m[1] - 0.00001) && (Tribarycent_m[i](1) < mincoords_m[1] + 0.00001))) */
	    /*	Tribarycent_m[i] = (geo3Dcoords_m[allbfaces_m[4*i+1]] + geo3Dcoords_m[allbfaces_m[4*i+2]] + geo3Dcoords_m[allbfaces_m[4*i+3]]) / 3.0 ;*/
	    double z_upper = maxcoords_m[2]/2.0 + 0.02;
	    double z_lower = maxcoords_m[2]/2.0 - 0.02;
	    double k_p = (0.733901-0.136929)/(1.84435-0.0286504);
	    double k_n = -1.0*k_p;
	    double y_ptest = k_p*(Tribarycent_m[i](0)-0.0286504)+0.136929;
	    double y_ntest = k_n*(Tribarycent_m[i](0)-0.0286504)-0.136929;
	    if(((Tribarycent_m[i](2) > z_lower - 0.00001) && (Tribarycent_m[i](2) < z_upper + 0.00001))
	       && (((Tribarycent_m[i](1) > y_ptest - 0.00001) && (Tribarycent_m[i](1) < y_ptest + 0.00001))
		   || ((Tribarycent_m[i](1) > y_ntest - 0.00001) && (Tribarycent_m[i](1) < y_ntest + 0.00001)))) {
		temp = BGphysics::Absorption;
	    } else {
		temp = ((BGphysics::Absorption) | (BGphysics::FNEmission) | (BGphysics::SecondaryEmission)); //Other positions can have obsorption, FNEmission and SecondaryEmission.
	    }
	    TriBGphysicstag_m.push_back(temp);
	}
    }
}

/**
 * Determine physical behaviour when particle hits the boundary triangle, non secondary emission version.
 */
int BoundaryGeometry::doBGphysics(const Vector_t &intecoords, const int &triId) {
    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if((BGtag & (BGphysics::Nop)) == (BGphysics::Nop)) {
        ret = -1;
    } else if ( ((BGtag & (BGphysics::Absorption)) == (BGphysics::Absorption)) && ((BGtag & (BGphysics::FNEmission)) != (BGphysics::FNEmission))) {

	ret = 0;

    } else {
	ret = 1;
    }
    return ret;
}



/**
 * Determine physical behaviour when particle hits the boundary triangle, call Furman-Pivi's secondary emission model.
 */


int BoundaryGeometry::doBGphysics(const Vector_t &intecoords, const int &triId, const double &incEnergy, const double &incQ, const Vector_t &incMomentum, PartBunch *itsBunch, double &seyNum) {
    Inform msg("BGphyDebug");
    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if((BGtag & (BGphysics::Nop)) == (BGphysics::Nop)) {
        ret = -1;
    } else if ( ((BGtag & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((BGtag & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((BGtag & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))   ) {

	ret = 0;

    }else {
        // Secondary Emission;
	ret = 1;
        int se_Num = 0;
        int seType = 0;
        if((BGtag & (BGphysics::SecondaryEmission)) == (BGphysics::SecondaryEmission)) {

            double cosTheta = -dot(incMomentum, TriNormal_m[triId]) / sqrt(dot(incMomentum, incMomentum)); //cosTheta should be positive
            if(cosTheta < 0) {
                cout << "cosTheta = " << cosTheta <<" intecoords " << intecoords(0)<<" "<<intecoords(1)<<" "<<intecoords(2)<<" "<< endl;
                cout<<"incident momentum=(" << incMomentum(0) << ","<<incMomentum(1)<<","<<incMomentum(2)<<") triNormal=(" << TriNormal_m[triId](0) <<","<< TriNormal_m[triId](1)<<","<< TriNormal_m[triId](2)<< ") "<< endl;
            }
            //assert(cosTheta>=0);
            if((intecoords[0] != geo3Dcoords_m[allbfaces_m[4*triId+1]](0))
               || (intecoords[1] != geo3Dcoords_m[allbfaces_m[4*triId+1]](1))
               || (intecoords[2] != geo3Dcoords_m[allbfaces_m[4*triId+1]](2))) {// intersection is not the 1st vertex
                Vector_t local_x = geo3Dcoords_m[allbfaces_m[4*triId+1]];
                nSec(incEnergy, cosTheta, seBoundaryMatType_m, se_Num, seType, incQ, TriNormal_m[triId], intecoords, local_x, itsBunch, seyNum, ppVw_m, vVThermal_m, nEmissionMode_m);

            } else {// intersection is the 1st vertex

                Vector_t local_x = geo3Dcoords_m[allbfaces_m[4*triId+2]];
                nSec(incEnergy, cosTheta, seBoundaryMatType_m, se_Num, seType, incQ, TriNormal_m[triId], intecoords, local_x, itsBunch, seyNum, ppVw_m, vVThermal_m, nEmissionMode_m);


            }
        }

    }
    return ret;
}



/**
 * Determine physical behaviour when particle hits the boundary triangle, call Vaughan's secondary emission model.
 */


int BoundaryGeometry::doBGphysics(const Vector_t &intecoords, const int &triId, const double &incEnergy, const double &incQ, const Vector_t &incMomentum, PartBunch *itsBunch, double &seyNum, const int &para_null) {
    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if((BGtag & (BGphysics::Nop)) == (BGphysics::Nop)) {
        ret = -1;
    } else if ( ((BGtag & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((BGtag & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((BGtag & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))   ) {

	ret = 0;

    } else {
        // Secondary Emission;

        int se_Num = 0;
        int seType = 0;
        if((BGtag & (BGphysics::SecondaryEmission)) == (BGphysics::SecondaryEmission)) {

            double cosTheta = -dot(incMomentum, TriNormal_m[triId]) / sqrt(dot(incMomentum, incMomentum)); //cosTheta should be positive
            if(cosTheta < 0) {
                cout << "cosTheta = " << cosTheta << endl;
                INFOMSG("incident momentum=" << incMomentum << " triNormal=" << TriNormal_m[triId] << " dot=" << dot(incMomentum, TriNormal_m[triId]) << endl);
            }
            //assert(cosTheta>=0);
            if((intecoords[0] != geo3Dcoords_m[allbfaces_m[4*triId+1]](0))
               || (intecoords[1] != geo3Dcoords_m[allbfaces_m[4*triId+1]](1))
               || (intecoords[2] != geo3Dcoords_m[allbfaces_m[4*triId+1]](2))) {// intersection is not the 1st vertex
                Vector_t local_x = geo3Dcoords_m[allbfaces_m[4*triId+1]];
                nSec(incEnergy, cosTheta, se_Num, seType, incQ, TriNormal_m[triId], intecoords, local_x, itsBunch, seyNum, ppVw_m, vSeyZero_m, vEzero_m, vSeyMax_m, vEmax_m, vKenergy_m, vKtheta_m, vVThermal_m, nEmissionMode_m);

            } else {// intersection is the 1st vertex

                Vector_t local_x = geo3Dcoords_m[allbfaces_m[4*triId+2]];
                nSec(incEnergy, cosTheta, se_Num, seType, incQ, TriNormal_m[triId], intecoords, local_x, itsBunch, seyNum, ppVw_m, vSeyZero_m, vEzero_m, vSeyMax_m, vEmax_m, vKenergy_m, vKtheta_m, vVThermal_m, nEmissionMode_m);


            }
        }

    }
    return ret;
}

/**
 * Here we call field emission model.
 */


void BoundaryGeometry::doFNemission(OpalBeamline &itsOpalBeamline, PartBunch *itsBunch, const double t) {//Self-field is not considered at moment. Only 1D Child-Langmuir law is implemented for space charge limited current density.
    const double fa = parameterFNA_m / workFunction_m * fieldEnhancement_m * fieldEnhancement_m;
    /*  int node_num = Ippl::getNodes();

    size_t *count = new size_t [node_num];
    // itsBunch->getLocalNum();
    for(int i = 0; i < node_num; i++) {

    count[i] = 0;

    }*/
    size_t Nstp = 0;
    for(int i = 0; i < numbfaces_global_m; i++) {

        if((TriBGphysicstag_m[i] & (BGphysics::FNEmission)) == (BGphysics::FNEmission)) {

            Vector_t E = (0.0, 0.0, 0.0), B = (0.0, 0.0, 0.0);
            Vector_t centroid = (0.0, 0.0, 0.0);
            unsigned long rtvtmp = itsOpalBeamline.getFieldAt(Tribarycent_m[i], centroid,  t, E, B);
            double Enormal = dot(TriNormal_m[i], E);
            //
            // Enormal should be negative as E field direction should be opposite to inward normal of surface
            //
            if(Enormal < fieldFNthreshold_m) {
                vector<Vector_t> vertex;
                vertex.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]]);
                vertex.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]]);
                vertex.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]]);
                Fieldemission(itsBunch, fa, Enormal, parameterFNB_m, workFunction_m, parameterFNVYZe_m, parameterFNVYSe_m, parameterFNY_m, fieldEnhancement_m, maxFNemission_m, Triarea_m[i], vertex, TriNormal_m[i], Nstp);
            }

        }

    }
    *gmsg<<"* Newly field emitted particle number= " << Nstp << endl;
}

void BoundaryGeometry::setNEmissionMode(bool nEmissionMode) {
    nEmissionMode_m = nEmissionMode;
}

void BoundaryGeometry::setWorkFunction(double workFunction) {
    workFunction_m = workFunction;
}

void BoundaryGeometry::setFieldEnhancement(double fieldEnhancement) {
    fieldEnhancement_m = fieldEnhancement;
}

void BoundaryGeometry::setMaxFN(size_t maxFNemission) {
    maxFNemission_m = maxFNemission;
}

void BoundaryGeometry::setFNTreshold(double fieldFNthreshold) {
    fieldFNthreshold_m = -1.0e6 * fieldFNthreshold;
}

void BoundaryGeometry::setFNParameterA(double parameterFNA) {
    parameterFNA_m = parameterFNA;
}

void BoundaryGeometry::setFNParameterB(double parameterFNB) {
    parameterFNB_m = parameterFNB;
}

void BoundaryGeometry::setFNParameterY(double parameterFNY) {
    parameterFNY_m = parameterFNY;
}

void BoundaryGeometry::setFNParameterVYZe(double parameterFNVYZe) {
    parameterFNVYZe_m = parameterFNVYZe;
}

void BoundaryGeometry::setFNParameterVYSe(double parameterFNVYSe) {
    parameterFNVYSe_m = parameterFNVYSe;
}

void BoundaryGeometry::setBoundaryMatType(int BoundaryMatType) {
    seBoundaryMatType_m = BoundaryMatType;
}

void BoundaryGeometry::setEInitThreshold(double einitthreshold) {eInitThreshold_m = 1.0e6 * einitthreshold;}
void BoundaryGeometry::setvSeyZero(double vSeyZero) {vSeyZero_m = vSeyZero;}
void BoundaryGeometry::setvEZero(double vEZero) {vEzero_m = vEZero;}
void BoundaryGeometry::setvSeyMax(double vSeyMax) {vSeyMax_m = vSeyMax;}
void BoundaryGeometry::setvEmax(double vEmax){vEmax_m = vEmax;}
void BoundaryGeometry::setvKenergy(double vKenergy){vKenergy_m = vKenergy;}
void BoundaryGeometry::setvKtheta(double vKtheta){vKtheta_m = vKtheta;}
void BoundaryGeometry::setvVThermal(double vVThermal){vVThermal_m = vVThermal;}
void BoundaryGeometry::setVw(double ppVw){ppVw_m = ppVw;}
/**
 * Determine if a line segment  has intersects with a triangle.
 * Algorithm :http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
 * @param x0 and @param x1 define the line segment.
 * @param i is the id of the triangle.
 */

Vector_t BoundaryGeometry::LineInsTri(Vector_t x0, Vector_t x1, size_t i) {
    //    INFOMSG(" LineInsTri called "<< endl);
    Vector_t lseg = x1 - x0;
    Vector_t t0, t1, t2;
    t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
    t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
    t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
    Vector_t u = t1 - t0;
    Vector_t v = t2 - t0;
    Vector_t lt = t0 - x0;
    Vector_t n;
    n[0] = u[1] * v[2] - v[1] * u[2];
    n[1] = u[2] * v[0] - v[2] * u[0];
    n[2] = u[0] * v[1] - v[0] * u[1];

    if(fabs(dot(n,lseg)) < 1.0e-10) {
        // cout<<"Triangle parallell to line segment, return x1 "<< endl;
        return x1;
    } else {

        double rI = dot(n,lt) / dot(n, lseg);
        Vector_t ItSec = x0 + rI * lseg;
        Vector_t w = ItSec - t0;
        if((rI < 0) || (rI > 1)) {
            //   cout<<"Intesect is on the extended line, return x1 "<< endl;
            return x1;
        } else {
            double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
            tmp1 = dot(u, v);
            tmp2 = dot(w, v);
            tmp3 = dot(u, w);
            tmp4 = dot(u, u);
            tmp5 = dot(v, v);
	    tmp6 = (tmp1 * tmp1 - tmp4 * tmp5);
            double sI = (tmp1 * tmp2 - tmp5 * tmp3) / tmp6;
            double tI = (tmp1 * tmp3 - tmp4 * tmp2) / tmp6;
            if((sI >= 0) && (tI >= 0) && ((sI + tI) <= 1)) {
                return ItSec;
            } else {
                //     cout<<"Intesect is on the extended plane, return x1 "<< endl;
                return x1;
            }
        }
    }
}



/**
 * Determine if a point is outside (return false), inside or just on (return true) the boundary.
 * @param x stands for coordinates of the test point.
 * The basic idea is if a line segment starting from the test point has odd intersects with a closed boundary,
 * then the test point is inside the geometry; if the intersects have even number, then the test points is outside the geometry;
 * if the test point is amoung the intersects, then the test point is just on the boundary. Makesure the end point of the line
 * segment is outside the geometry boundary.
 */

bool  BoundaryGeometry::isInside(Vector_t x) {
    Vector_t x0 = x;
    Vector_t x1;
    x1[0] = x0[0];
    //x1[1] = x0[1];
    RANLIB_class *rGen = new RANLIB_class(265314159, 4);
    x1[1] = maxcoords_m[1] * (1.1 + rGen->uniform(0.0, 1.0));
    x1[2] = maxcoords_m[2] * (1.1 + rGen->uniform(0.0, 1.0));
    //x1[2] = x0[2];
    delete rGen;
    //*gmsg<<"x1: "<<x1<<endl;

    /*=================================================debug=====================================================================
    std::ofstream of;
    of.open(string("vtk/testLine.vtk").c_str());
    assert(of.is_open());
    of.precision(6);
    of << "# vtk DataFile Version 2.0" << endl;                      // first line of a vtk file
    of << "generated using BoundaryGeometry::makeBoundaryIndexSet" << endl;   // header
    of << "ASCII" << endl << endl;                                   // file format
    of << "DATASET UNSTRUCTURED_GRID" << endl;                 // dataset structrue
    of << "POINTS " << 2 << " float" << endl;  // data num and type
   
    of << x0[0] << " " << x0[1]  << " " << x0[2] << endl;
    
    of << x1[0] << " " << x1[1]  << " " << x1[2] << endl;
       
    

    of << "CELLS " << 1 << " " << 3 << endl;
   
    of << "2 " << 0 << " " << 1 << endl;
    of << "CELL_TYPES " << 1 << endl;
    of << "3" << endl;
    of << endl;

    ====================================================================================================*/
    /*
     * Random number could avoid some specific situation,
     * like line parallel to boundary......
     * x1 could be any point outside the boundary ;
     */

    vector<Vector_t> IntesecNum = PartBoundaryInteNum(x0, x1);
    /* vector<Vector_t>::iterator InIt;
    for(InIt = IntesecNum.begin(); InIt != IntesecNum.end(); InIt++) {
	*gmsg<<"IntesecNum " << *InIt <<" IntesecNum.size() "<<IntesecNum.size()<<endl;
	}*/
    if(IntesecNum[0] == x0) {
	return true;// x0 is just on the boundary;
    } else {
        if(((IntesecNum.size() % 2) == 0) || (*IntesecNum.begin() == x1)) {
            return false;// x0 is  outside the boundary;
        } else
            return true;// x0 is inside the boundary;
    }
}

vector<Vector_t>  BoundaryGeometry::GridIntersection(Vector_t x0, Vector_t x1) {

    using namespace std;

    vector<Vector_t> SegDiscrete;
    vector<Vector_t> Isp;
    vector<Vector_t>::iterator TriIscIt;

    vector<int> TriId;
    vector<int>::iterator TriIdIt;
    double hr_tmp=hr_m[0];
    if(hr_m[1]<hr_tmp)
        hr_tmp=hr_m[1];
    if(hr_m[2]<hr_tmp)
        hr_tmp=hr_m[2];//interval should be smaller than the size of box;
    //hr_tmp*=0.001;
    int Seglen = (((int)floor(sqrt(dot(x0 - x1, x0 - x1)) / hr_tmp)) + 1);
    int count = 0;
    for(int i = 0; i < Seglen ; i++) {
        SegDiscrete.push_back(x0 + hr_tmp*i * (x1 - x0) / sqrt(dot(x0 - x1, x0 - x1)));
    }
    SegDiscrete.push_back(x1);


    for(vector<Vector_t>::iterator myit = SegDiscrete.begin() ; myit != SegDiscrete.end() ; myit++) {

        size_t id = f(*myit);

        vector<size_t> ::iterator idIt;
        map< size_t, vector<size_t> >::iterator It = CubicLookupTable_m.find(id);
        if(It != CubicLookupTable_m.end()) {
            /*
             * If cubic box is the boundary bounding box then do the following tests.
             */
            for(idIt = (*It).second.begin(); idIt != (*It).second.end(); idIt++) {
                Vector_t tmp =  LineInsTri(x0, x1, *idIt);
                TriIdIt = std::find(TriId.begin(), TriId.end(), *idIt);
                if(((tmp != x1) && ((TriIdIt == TriId.end()) || TriId.size() == 0))) {
                    TriIscIt = std::find(Isp.begin(), Isp.end(), tmp);
                    if((TriIscIt == Isp.end()) || (Isp.size() == 0)) {
                        TriId.push_back(*idIt);
                        Isp.push_back(tmp);
                        ++count;
                    }
                }
            }
        }

    }
    if(count == 0)
        Isp.push_back(x1);
    return Isp;
}



void  BoundaryGeometry::initialize() {
    *gmsg<<"* Iniitializing Boundary Geometry ... ..."<<endl;
    IpplTimings::startTimer(TPreProc_m);
    string dir("vtk");
    mkdir(dir.c_str(), O_CREAT);
    chmod(dir.c_str(), S_IRWXU);
   
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS); //Property list identifier
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    // Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fopen(h5FileName_m.c_str(), H5F_ACC_RDONLY, plist_id);
    assert(file_id >= 0);
    H5Pclose(plist_id);
   
    /////////////////////////////////////////////
    //   Read dataset "surface" from .h5 file
    ////////////////////////////////////////////

    hsize_t dimsf[2];//dataset dimensions
    herr_t  status;

    hid_t dset_id = H5Dopen(file_id, "/surface");
    assert(dset_id >= 0);
    // Create the dataspace for the dataset.
    hid_t x = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(x, dimsf, NULL);
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);

    numbfaces_global_m = dimsf[0];
    allbfaces_m = new int [numbfaces_global_m * 4];

   
    Tribarycent_m = new Vector_t[numbfaces_global_m];
    // Create property list for collective dataset write.
    hid_t plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id2, H5FD_MPIO_COLLECTIVE);
    status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace, plist_id2, allbfaces_m);
    assert(status >= 0);
    // Store local boundary faces, discard others
    nof_sym_m = 0; // Count number of symmetry planes.
    int const nogeosym_flag = 0; // heronion outputs a index, case 0 indicates no symmetry plane, 1 for (x,y) sym, 2 for (y,z), 3 for (z,x),none 0 means the current triangle is not on a real surface of geometry but on a symmetry plane.

   
    for(int i = 0; i < numbfaces_global_m; i ++) {
	if(allbfaces_m[4 * i] > nogeosym_flag) {
            nof_sym_m += 1;
            if(i < numbfaces_global_m - 1) {
                for(int j = 0; j < numbfaces_global_m - i; j ++) {
                    allbfaces_m[4*(i+j)] = allbfaces_m[4*(i+j+1)];
                    allbfaces_m[4*(i+j)+1] = allbfaces_m[4*(i+j+1)+1];
                    allbfaces_m[4*(i+j)+2] = allbfaces_m[4*(i+j+1)+2];
                    allbfaces_m[4*(i+j)+3] = allbfaces_m[4*(i+j+1)+3];
                }
            } else
                numbfaces_global_m = numbfaces_global_m - 1;
        }
    }
    H5Dclose(dset_id);
    H5Sclose(filespace);
    *gmsg<<"*  Surface index read in done."<<endl;
    ////////////////////////////////
    //  Also read dataset "coords"
    ////////////////////////////////
    hsize_t dimsf_c[2];
    herr_t status_c;
    hid_t dset_id_c = H5Dopen(file_id, "/coords");
    assert(dset_id_c >= 0);

    // Create the dataspace for the dataset.
    hid_t cp_space = H5Dget_space(dset_id_c);
    H5Sget_simple_extent_dims(cp_space, dimsf_c, NULL);
    hid_t filespace_c = H5Screate_simple(2, dimsf_c, NULL);

    numpoints_global_m = dimsf_c[0];
    double *point_coords = new double[3 * numpoints_global_m];
    hid_t plist_id3 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id3, H5FD_MPIO_COLLECTIVE);
    status_c = H5Dread(dset_id_c, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_c, plist_id3, point_coords);
    assert(status_c >= 0);
    Vector_t geo3d_tmp;

    //int zflg = getZflag();//flag control the changing of axis.
    double zshift = getZshift();
    //if (zflg == 0) {

    for(int i = 0; i < numpoints_global_m; i++) {
	geo3d_tmp[0] = point_coords[3*i];
	geo3d_tmp[1] = point_coords[3*i+1];
	geo3d_tmp[2] = point_coords[3*i+2] + zshift;
	geo3Dcoords_m.push_back(geo3d_tmp);
    }
    *gmsg<<"*  Vertex built done."<<endl;
    mincoords_m= getMinExtend();
    TriPrPartloss_m = new double [numbfaces_global_m];
    TriFEPartloss_m = new double [numbfaces_global_m];
    TriSePartloss_m = new double [numbfaces_global_m];
    for(int i = 0; i < numbfaces_global_m; i ++) {
        Tribarycent_m[i] = (geo3Dcoords_m[allbfaces_m[4*i+1]] + geo3Dcoords_m[allbfaces_m[4*i+2]] + geo3Dcoords_m[allbfaces_m[4*i+3]]) / 3.0 ;
        Triarea_m.push_back(TriangleArea(i));

        TriPrPartloss_m[i] = 0.0;
        TriFEPartloss_m[i] = 0.0;
	TriSePartloss_m[i] = 0.0;
    }
    if(point_coords)
        delete point_coords;
    // Close HDF5 stuff
    H5Dclose(dset_id_c);
    H5Sclose(filespace_c);
    *gmsg<<"*  Triangle barycent built done."<<endl;  
    maxcoords_m = getMaxExtend();
  
    len_m = maxcoords_m - mincoords_m;  
    getMaxDimenssion();// get maximum and minimum of triangle side.

    // In principal, the number of discretization nr_m is maxmum lenth in each dimension divided by maxmum of triangle length. But if the hot spot, i.e., the multipacting/field emission zone is too small that normal bounding box covers the whole hot spot, the expensive triangle-line intersection tests will be frequently called. In these cases, we need to use smaller bounding box size to speed up simulation by setting a larger nr_m, for example use nr_m(0) /= 0.2 in the following. 
    // Todo: The relation between bounding box size and simulation time step & geometry shape maybe need to be summarized and modeled in a more flexible manner and could be adjusted in input file.
    nr_m(0) = (int)floor(len_m(0) / triangle_max_m / 0.2);//2.0 
    nr_m(1) = (int)floor(len_m(1) / triangle_max_m / 0.2);
    nr_m(2) = (int)floor(len_m(2) / triangle_max_m / 0.1);//2.0

    hr_m = len_m / nr_m;
    out_m = maxcoords_m + hr_m;
    TriPrPartlossZ_m = new double [nr_m(2)];
    TriSePartlossZ_m = new double [nr_m(2)];
    TriFEPartlossZ_m = new double [nr_m(2)];
     *gmsg<<"*  Geometry interval built done."<<endl; 
    makeBoundaryIndexSet();
     *gmsg<<"*  BoundaryIndexSet built done."<<endl; 
    makeTriNormal();
     *gmsg<<"*  Triangle Normal built done."<<endl;
    setBGphysicstag();
    *gmsg<<*this<<endl;
    IpplTimings::stopTimer(TPreProc_m);
}


/**
 * Initialize some darkcurrent particles near the surface with inward momenta.
 */
void  BoundaryGeometry::createParticlesOnSurface(size_t n, double darkinward,  OpalBeamline &itsOpalBeamline, PartBunch &itsBunch) {
    int tag = 1002;
    int Parent=0;
    if (Ippl::myNode()==0) {
	for(int i = 0; i < n; i ++) {
	    short BGtag = BGphysics::Absorption;
	    int k = 0;
	    Vector_t E = (0.0, 0.0), B = (0.0, 0.0);
	    while ( (((BGtag & (BGphysics::Absorption)) == (BGphysics::Absorption)) 
		     && ((BGtag & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) 
		     && ((BGtag & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))
		    || ( (fabs(E(0)) < eInitThreshold_m) 
			 && (fabs(E(1)) < eInitThreshold_m) 
			 && (fabs(E(2)) < eInitThreshold_m) ) ) {
		
		
		
		E = (0.0, 0.0, 0.0), B = (0.0, 0.0, 0.0);
		int tmp = (int)(IpplRandom() * numbfaces_global_m) ;
		BGtag = TriBGphysicstag_m[tmp];
		k = tmp;
		Vector_t centroid = (0.0, 0.0);
		const unsigned long sindex = 0;
		unsigned long temp = itsOpalBeamline.getFieldAt(Tribarycent_m[k]+darkinward * TriNormal_m[k], centroid, itsBunch.getdT(), E, B);
		
	    }
	    partsr_m.push_back(Tribarycent_m[k] + darkinward * TriNormal_m[k]);
	    
	}
	Message *mess = new Message();
	putMessage(*mess, partsr_m.size());
	for (vector<Vector_t>::iterator myIt = partsr_m.begin(); myIt != partsr_m.end(); myIt++) {
	    putMessage(*mess,*myIt);
	    
	}
	Ippl::Comm->broadcast_all(mess,tag);
	
	
    } else {
	// receive particle position message
	size_t nData = 0;
	Message *mess = Ippl::Comm->receive_block(Parent,tag);
	getMessage(*mess,nData);
	for(size_t i=0; i<nData; i++) {
	    Vector_t tmp = Vector_t(0.0);
	    getMessage(*mess,tmp);
	    partsr_m.push_back(tmp);
	}
	
    }
    
}
    



/**
 * Initialize primary particles near the surface with inward momenta.
 */
void  BoundaryGeometry::createPriPart(size_t n, double darkinward,  OpalBeamline &itsOpalBeamline, PartBunch *itsBunch) {
  
    
    if( Options::ppdebug ) {
	int tag = 1001;
	int Parent=0;
	if (Ippl::myNode()==0) {
	      
	    double x_low = mincoords_m(0)+0.5*len_m(0)-0.49*len_m(0);// limit the initial particle in the center of the lower parallel plate. There is a distance of 0.01*length in x direction as margin.
	    double x_up = mincoords_m(0)+0.5*len_m(0)+0.49*len_m(0);// limit the initial particle in the center of the upper parallel plate. There is a distance of 0.01*length in x direction as margin.
	    double y_low = mincoords_m(1)+0.5*len_m(1)-0.49*len_m(1);// limit the initial particle in the center of the lower parallel plate. There is a distance of 0.01*length in y direction as margin.
	    double y_up = mincoords_m(1)+0.5*len_m(1)+0.49*len_m(1);// limit the initial particle in the center of the upper parallel plate. There is a distance of 0.01*length in y direction as margin.
	    for(int i = 0; i < n/2; i ++) {
            
		double zCoord = maxcoords_m(2);
		double xCoord = maxcoords_m(0);
		double yCoord = maxcoords_m(1);
		while ((zCoord>0.000001) 
		       || (zCoord<-0.000001)
		       || (xCoord>x_up)
		       || (xCoord<x_low)
		       || (yCoord>y_up)
		       || (yCoord<y_low)) {
			
		    int k = (int)(IpplRandom() * numbfaces_global_m) ;
		    zCoord = Tribarycent_m[k](2);
		    xCoord = Tribarycent_m[k](0);
		    yCoord = Tribarycent_m[k](1);
		    if ((Tribarycent_m[k](2)<0.000001)&&(Tribarycent_m[k](2)>-0.000001)&&(Tribarycent_m[k](0)<x_up)&&(Tribarycent_m[k](0)>x_low)&&(Tribarycent_m[k](1)<y_up)&&(Tribarycent_m[k](1)>y_low)) {//
			partsr_m.push_back(Tribarycent_m[k] + darkinward * TriNormal_m[k]);
			partsp_m.push_back(TriNormal_m[k]);
		    }
		    
		}
		    
	    }
	    for(int i = 0; i < n/2; i ++) {
		
		double zCoord = maxcoords_m(2);
		double xCoord = maxcoords_m(0);
		double yCoord = maxcoords_m(1);
		while ((zCoord>(maxcoords_m(2)+0.000001)) 
		       || (zCoord<(maxcoords_m(2)-0.000001))
		       || (xCoord>x_up)
		       || (xCoord<x_low)
		       || (yCoord>y_up)
		       || (yCoord<y_low)) {
		    
		    int k = (int)(IpplRandom() * numbfaces_global_m) ;
		    zCoord = Tribarycent_m[k](2);
		    xCoord = Tribarycent_m[k](0);
		    yCoord = Tribarycent_m[k](1);
		    if ((Tribarycent_m[k](2)<maxcoords_m(2)+0.000001)&&(Tribarycent_m[k](2)>maxcoords_m(2)-0.000001)&&(Tribarycent_m[k](0)<x_up)&&(Tribarycent_m[k](0)>x_low)&&(Tribarycent_m[k](1)<y_up)&&(Tribarycent_m[k](1)>y_low)) {//
			partsr_m.push_back(Tribarycent_m[k] + darkinward * TriNormal_m[k]);
			partsp_m.push_back(TriNormal_m[k]);
		    }
		    
		}
		    
	    }
	   
	    Message *mess = new Message();
	    putMessage(*mess, partsr_m.size());
	    for (vector<Vector_t>::iterator myIt = partsr_m.begin(), myItp = partsp_m.begin(); myIt != partsr_m.end(); myIt++,myItp++) {
		putMessage(*mess,*myIt);
		putMessage(*mess,*myItp);
	    }
	    Ippl::Comm->broadcast_all(mess,tag);
	} else {
	    // receive particle position message
	    size_t nData = 0;
	    Message *mess = Ippl::Comm->receive_block(Parent,tag);
	    getMessage(*mess,nData);
	    for(size_t i=0; i<nData; i++) {
		Vector_t tmpr = Vector_t(0.0);
		Vector_t tmpp = Vector_t(0.0);
		getMessage(*mess,tmpr);
		getMessage(*mess,tmpp);
		partsr_m.push_back(tmpr);
		partsp_m.push_back(tmpp);
	    }
	}
        
    } else {
	int tag = 1001;
	int Parent=0;
	if (Ippl::myNode()==0) {
	    for(int i = 0; i < n; i ++) {
		short BGtag = BGphysics::Absorption;
		int k = 0;
		Vector_t E = (0.0, 0.0), B = (0.0, 0.0);
		while ( (((BGtag & (BGphysics::Absorption)) == (BGphysics::Absorption)) 
			&& ((BGtag & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) 
			 && ((BGtag & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))
			|| ( (fabs(E(0)) < eInitThreshold_m) 
			     && (fabs(E(1)) < eInitThreshold_m) 
			     && (fabs(E(2)) < eInitThreshold_m) ) ) {
		   

		    E = (0.0, 0.0, 0.0), B = (0.0, 0.0, 0.0);
		    int tmp = (int)(IpplRandom() * numbfaces_global_m) ;
		    BGtag = TriBGphysicstag_m[tmp];
		    k = tmp;
		    Vector_t centroid = (0.0, 0.0);
		    const unsigned long sindex = 0;
		    unsigned long temp = itsOpalBeamline.getFieldAt(Tribarycent_m[k]+darkinward * TriNormal_m[k], centroid, itsBunch->getdT(), E, B);
		   
		    
		}
		partsr_m.push_back(Tribarycent_m[k] + darkinward * TriNormal_m[k]);
	
	    }
	    Message *mess = new Message();
	    putMessage(*mess, partsr_m.size());
	    for (vector<Vector_t>::iterator myIt = partsr_m.begin(); myIt != partsr_m.end(); myIt++) {
		putMessage(*mess,*myIt);

	    }
	    Ippl::Comm->broadcast_all(mess,tag);
	  
	  
	} else {
	    // receive particle position message
	    size_t nData = 0;
	    Message *mess = Ippl::Comm->receive_block(Parent,tag);
	    getMessage(*mess,nData);
	    for(size_t i=0; i<nData; i++) {
		Vector_t tmp = Vector_t(0.0);
		getMessage(*mess,tmp);
		partsr_m.push_back(tmp);
	    }
	  
	}
      

    }
}




/**
 * Make the boundary set by using triangle vertex, bounding box vertex ,as well as points in triangle central lines.
 */


void BoundaryGeometry::makeBoundaryIndexSet() {
    set<size_t>::iterator bbIt;
    vector<Vector_t> BBox;
    for(int i = 0; i < numbfaces_global_m; i++) {
        vector<Vector_t> discreteTri;
        vector<Vector_t>::iterator disIt;
        for(int j = 0; j < 200; j++) { // Discrete three central lines and three triangle sides to 200 segments to get a more complete boundary index set.
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]] + 0.005 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+2]] + geo3Dcoords_m[allbfaces_m[4*i+3]]) - geo3Dcoords_m[allbfaces_m[4*i+1]]));
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]] + 0.005 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+3]] + geo3Dcoords_m[allbfaces_m[4*i+1]]) - geo3Dcoords_m[allbfaces_m[4*i+2]]));
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]] + 0.005 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+1]] + geo3Dcoords_m[allbfaces_m[4*i+2]]) - geo3Dcoords_m[allbfaces_m[4*i+3]]));// Discrete three central lines.
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]] + 0.005 * (j + 1) * (geo3Dcoords_m[allbfaces_m[4*i+2]] - geo3Dcoords_m[allbfaces_m[4*i+1]]));
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]] + 0.005 * (j + 1) * (geo3Dcoords_m[allbfaces_m[4*i+3]] - geo3Dcoords_m[allbfaces_m[4*i+2]]));
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]] + 0.005 * (j + 1) * (geo3Dcoords_m[allbfaces_m[4*i+1]] - geo3Dcoords_m[allbfaces_m[4*i+3]]));// Discrete three  triangle sides.
        }
	discreteTri.push_back(Tribarycent_m[i]);
        for(disIt = discreteTri.begin(); disIt != discreteTri.end(); disIt++) {
            int id = f(*disIt);

            //  assert(id > 0);
            bbIt = boundary_ids_m.find(id);
	    if((bbIt == boundary_ids_m.end())) {
                Vector_t temp;
                temp[0] = floor(((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];
                temp[1] = floor(((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor(((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back(temp);
                boundary_ids_m.insert(id);
            }


            vector<size_t> tmp;
            map< size_t, vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find(id);
            if(It == CubicLookupTable_m.end()) {
                tmp.push_back(i);
                CubicLookupTable_m.insert(pair <size_t, vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back(i);
        }
	vector<Vector_t> ret, cubic_coords;
        ret = SetMinMaxBound(i);
        cubic_coords.push_back(ret[0]);
        cubic_coords.push_back((ret[0](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back((ret[0](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back((ret[0](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back((ret[1](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back((ret[1](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back((ret[1](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back((ret[1](0), ret[0](1), ret[0](2)));

	cubic_coords.push_back(ret[2]);
        cubic_coords.push_back((ret[2](0), ret[2](1), ret[3](2)));
        cubic_coords.push_back((ret[2](0), ret[3](1), ret[3](2)));
        cubic_coords.push_back((ret[2](0), ret[3](1), ret[2](2)));
        cubic_coords.push_back((ret[3](0), ret[3](1), ret[2](2)));
        cubic_coords.push_back((ret[3](0), ret[3](1), ret[3](2)));
        cubic_coords.push_back((ret[3](0), ret[2](1), ret[3](2)));
        cubic_coords.push_back((ret[3](0), ret[2](1), ret[2](2)));
        // vector<Vector_t>::iterator disIt;
        for(disIt = cubic_coords.begin(); disIt != cubic_coords.end(); disIt++) {
            int id = f(*disIt);
            bbIt = boundary_ids_m.find(id);

            //   assert(id > 0);


            if((bbIt == boundary_ids_m.end())) {
                Vector_t temp;
                temp[0] = floor(((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];// this may cause the increct display.
                temp[1] = floor(((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor(((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back(temp);
                boundary_ids_m.insert(id);
            }

            vector<size_t> tmp;
            map< size_t, vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find(id);
            if(It == CubicLookupTable_m.end()) {
                tmp.push_back(i);
                CubicLookupTable_m.insert(pair <size_t, vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back(i);
	}

    }
  
    
    /*-------------------------------------------------------------------------*/
    size_t numpoints = 8 * BBox.size();
    vector<Vector_t>::iterator bIt;
    std::ofstream of;
    of.open(string("vtk/testBBox.vtk").c_str());
    assert(of.is_open());
    of.precision(6);
    of << "# vtk DataFile Version 2.0" << endl;                      // first line of a vtk file
    of << "generated using BoundaryGeometry::makeBoundaryIndexSet" << endl;   // header
    of << "ASCII" << endl << endl;                                   // file format
    of << "DATASET UNSTRUCTURED_GRID" << endl;                 // dataset structrue
    of << "POINTS " << numpoints << " float" << endl;  // data num and type
    for(bIt = BBox.begin(); bIt != BBox.end() ; bIt ++) {
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1] << " " << (*bIt)[2] + hr_m[2] << endl;
    }
    of << endl;

    of << "CELLS " << BBox.size() << " " << 9 * BBox.size() << endl;
    for(int i = 0; i < BBox.size() ; i ++)
        of << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3 << " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << endl;
    of << "CELL_TYPES " << BBox.size() << endl;
    for(int i = 0; i <  BBox.size() ; i ++)
        of << "11" << endl;
    of << "CELL_DATA " << BBox.size() << endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << endl;
    of << "LOOKUP_TABLE " << "default" << endl ;
    for(int i = 0; i <  BBox.size() ; i ++)
        of << (float)(i) << endl;
    of << endl;
    of << "COLOR_SCALARS " << "BBoxColor " << 4 << endl ;
    for(int i = 0; i < BBox.size() ; i ++) {

        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << endl;

    }

    of << endl;
}

/**
 * Determine whether a particle with position @param r, momenta @param v, and time step @param dt will hit the boundary.
 *
 * Basic algorithms are as follows:
 * 1) Determine if the particle is near the boundary by checking r is in boundary bounding box index set.
 *    if r is in, then do the following checking step, else return -1 to indicate that particle is far from boundary and to be integrated.
 * 2) Traversal all the triangles in the bounding cubic box which the particle is in, as well as triangles in the adjacent 26 bounding cubic boxes
 *    if the momenta has oppsite direction with those triangles' normals, then check if the particle has intersection with those triangles. If
 *    intersection exsists, then return 0.
 */

int BoundaryGeometry::PartInside(const Vector_t r, const Vector_t v, const double dt, int Parttype, const double Qloss, Vector_t &intecoords, int &triId, double &Energy) {
    
    int ret;
    const double p_sq = dot(v, v);
    const double betaP = 1.0 / sqrt(1.0 + p_sq);
   
    const Vector_t temp1 = r;//particle position in tstep n;
    const Vector_t temp = r + ( c * betaP * v * dt); //particle position in tstep n+1;
    double rI = 0.0;
    Vector_t Isc = out_m;
    
    IpplTimings::startTimer(TPInside_m);

    //test if particle position in tstep n is inside the cubic bounding box. If true, do the following tests
    if(isInGeometry(temp1)) {
        int flg = 0;
        int id = f(temp1);
        
        // Build an array contains the ids of 27(3*3*3) cubic boxes. The id of the cubic box which contains the particle,
        // lies in the center of the 27 cubic boxes. Just for the situation that even the line segment cross more than one cubic box.
        int idc[27] = {id, id + 1, id - 1, id + nr_m(0), id - nr_m(0), id + nr_m(0) + 1, id + nr_m(0) - 1, id - nr_m(0) - 1, id - nr_m(0) + 1,
                       id + nr_m(0) *nr_m(1), id + nr_m(0) *nr_m(1) + 1, id + nr_m(0) *nr_m(1) - 1, id + nr_m(0) *nr_m(1) + nr_m(0),
                       id + nr_m(0) *nr_m(1) - nr_m(0), id + nr_m(0) *nr_m(1) + nr_m(0) + 1, id + nr_m(0) *nr_m(1) + nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1), id - nr_m(0) *nr_m(1) + 1, id - nr_m(0) *nr_m(1) - 1, id - nr_m(0) *nr_m(1) + nr_m(0), id - nr_m(0) *nr_m(1) - nr_m(0), id - nr_m(0) *nr_m(1) + nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1) + nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) + 1
	};

        // Test all the 27 cubic boxes to find if the line segment has intersection with the triangles in those cubic boxes.
        for(int k = 0; k < 27; k++) {


            map< size_t, vector<size_t> >::iterator It;
            vector<size_t> ::iterator myIt;
            It = CubicLookupTable_m.find(idc[k]);
           
            // Lookup table (stl map) contains the id of cubic box as key and triangle ids inside the id th cubic box as value
            if(It != CubicLookupTable_m.end()) { // If cubic box among the 27 cubic boxes is not the boundary bounding box, do not need to be tested at all. If cubic box is the boundary bounding box then do the following tests.
                for(myIt = (*It).second.begin(); myIt != (*It).second.end(); myIt++) {
                    // For each particle inside the cubic box, test if this particle have intersect with each triangle inside the same cubic box
                    if((v != 0) && (dot(v, TriNormal_m[*myIt]) <= 0.0)) {  // If the particle have none zero momenta and momenta has opposite direction with triangle normal, do the following tests.


                        IpplTimings::startTimer(TRayTrace_m);
                        FindIntersection(temp1, temp, *myIt, rI, Isc); // Ray tiangle test, parameters are particle position in tstep n, particle position in tstep n+1 and tiangle id, rI and Isc will return the ratio and intersection points.
                        IpplTimings::stopTimer(TRayTrace_m);

                        if((Isc != out_m)) {
                            // Test if the intersection is between the particle position in tstep n and particle position in tstep n+1 or is in the extension of line
                            // segment when particle position in tstep n is already outside the geometry( this may be not accurate and may be the source of problem.)
                            if(((rI >= -0.00001) && (rI <= 1.00001))  || ((rI < 0) && (dot(temp1 - Tribarycent_m[*myIt], TriNormal_m[*myIt]) <= 0.0))) {
				flg++;
                                intecoords = Isc;
                                triId = (*myIt);
                                assert((dot(TriNormal_m[*myIt], v) < 0)||(Isc==temp1));
                                Energy = Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9; //in eV
                                if(Parttype == 0) {
                                    TriPrPartloss_m[*myIt] += Qloss;
                                } else if(Parttype == 1) {
                                    TriFEPartloss_m[*myIt] += Qloss;
                                } else {
                                    TriSePartloss_m[*myIt] += Qloss;
                                }
				if (TriPrPartloss_m[*myIt]>0 || TriSePartloss_m[*myIt]>0 || TriFEPartloss_m[*myIt]>0) {

				    cout<<"* Loss Data"<<*myIt<<" : "<<TriPrPartloss_m[*myIt]<<" "<<TriSePartloss_m[*myIt]<<" "<<TriFEPartloss_m[*myIt]<<endl;

				}
                                break;
                            }
                        }
                    }
                }
                if(flg != 0) {
                    ret = 0;
                    break;
                }
                
            }
        }
        if(flg == 0) {
            ret = -1;
        }
    } else if(isInGeometry(temp)) {
        int flg = 0;
        int id = f(temp);
        
        // Build an array contains the ids of 27(3*3*3) cubic boxes. The id of the cubic box which contains the particle,
        // lies in the center of the 27 cubic boxes. Just for the situation that even the line segment cross more than one cubic box.
        int idc[27] = {id, id + 1, id - 1, id + nr_m(0), id - nr_m(0), id + nr_m(0) + 1, id + nr_m(0) - 1, id - nr_m(0) - 1, id - nr_m(0) + 1,
                       id + nr_m(0) *nr_m(1), id + nr_m(0) *nr_m(1) + 1, id + nr_m(0) *nr_m(1) - 1, id + nr_m(0) *nr_m(1) + nr_m(0),
                       id + nr_m(0) *nr_m(1) - nr_m(0), id + nr_m(0) *nr_m(1) + nr_m(0) + 1, id + nr_m(0) *nr_m(1) + nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1), id - nr_m(0) *nr_m(1) + 1, id - nr_m(0) *nr_m(1) - 1, id - nr_m(0) *nr_m(1) + nr_m(0), id - nr_m(0) *nr_m(1) - nr_m(0), id - nr_m(0) *nr_m(1) + nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1) + nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) + 1
	};

        // Test all the 27 cubic boxes to find if the line segment has intersection with the triangles in those cubic boxes.
        for(int k = 0; k < 27; k++) {


            map< size_t, vector<size_t> >::iterator It;
            vector<size_t> ::iterator myIt;
            It = CubicLookupTable_m.find(idc[k]);
           
            // Lookup table (stl map) contains the id of cubic box as key and triangle ids inside the id th cubic box as value
            if(It != CubicLookupTable_m.end()) { // If cubic box among the 27 cubic boxes is not the boundary bounding box, do not need to be tested at all. If cubic box is the boundary bounding box then do the following tests.
                for(myIt = (*It).second.begin(); myIt != (*It).second.end(); myIt++) {
                    // For each particle inside the cubic box, test if this particle have intersect with each triangle inside the same cubic box
                    if((v != 0) && (dot(v, TriNormal_m[*myIt]) <= 0.0)) {
                        // If the particle have none zero momenta and momenta has opposite direction with triangle normal, do the following tests.
                        IpplTimings::startTimer(TRayTrace_m);
                        //  Intersection tmp = FindIntersection(temp1, temp, *myIt); //Ray tiangle test, parameters are particle position in tstep n, particle position in tstep n+1 and tiangle id.
                        FindIntersection(temp1, temp, *myIt, rI, Isc);
                        IpplTimings::stopTimer(TRayTrace_m);
                        if((Isc != out_m)) {
                            // Test if the intersection is between the particle position in tstep n and particle position in tstep n+1 or is in the extension of line
                            // segment when particle position in tstep n is already outside the geometry( this may be not accurate and may be the source of problem.)
                            if(((rI >= -0.00001) && (rI <= 1.00001))  || ((rI < 0) && (dot(temp1 - Tribarycent_m[*myIt], TriNormal_m[*myIt]) <= 0.0))) {
				flg++;
                                intecoords = Isc;
                                triId = (*myIt);
                                assert(dot(TriNormal_m[*myIt], v) < 0);
                                Energy = Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9; //in eV
                                if(Parttype == 0)
                                    TriPrPartloss_m[*myIt] += Qloss;
                                else if(Parttype == 1)
                                    TriFEPartloss_m[*myIt] += Qloss;
				else 
                                    TriSePartloss_m[*myIt] += Qloss;
				if (TriPrPartloss_m[*myIt]>0 || TriSePartloss_m[*myIt]>0 || TriFEPartloss_m[*myIt]>0) {

				    cout<<"* Loss Data"<<*myIt<<" : "<<TriPrPartloss_m[*myIt]<<" "<<TriSePartloss_m[*myIt]<<" "<<TriFEPartloss_m[*myIt]<<" qloss: "<<Qloss<<endl;

				}
                                break;
                            }
                        }
                    }
                }

                if(flg != 0) {
                    ret = 0;
                    break;
                }
            }

        }
        if(flg == 0) {
            ret = -1;
        }
    } else {
        ret = -1;
    }
    IpplTimings::stopTimer(TPInside_m);
    return ret;
}



void BoundaryGeometry::updateElement(ElementBase *element) {

}



Inform &BoundaryGeometry::print(Inform &os) const {
    os << "* *************Boundary Geometry Info*********************************************** " << endl;
    os << "* GEOMETRY                   " << getOpalName() << '\n'
       << "* FGEOM                      " << Attributes::getString(itsAttr[FGEOM]) << '\n'
       << "* TOPO                       " << Attributes::getString(itsAttr[TOPO]) << '\n'
       << "* LENGHT                     " << Attributes::getReal(itsAttr[LENGHT]) << '\n'
       << "* S                          " << Attributes::getReal(itsAttr[S]) << '\n'
       << "* A                          " << Attributes::getReal(itsAttr[A]) << '\n'
       << "* B                          " << Attributes::getReal(itsAttr[B]) << '\n';
    if (getTopology()==string("BOXCORNER")) {
	  os << "* C                          " << Attributes::getReal(itsAttr[C]) << '\n'
	     << "* L1                         " << Attributes::getReal(itsAttr[L1]) << '\n'
	     << "* L1                         " << Attributes::getReal(itsAttr[L2]) << '\n';
	}
    os << "* Total triangle num         " << numbfaces_global_m <<'\n'
       << "* Total points num           " << numpoints_global_m <<'\n'
       << "* Triangle side(m)   Max=    " << triangle_max_m << '\n'
       << "*                    Min=    " << triangle_min_m <<'\n'
       << "* Geometry bounds(m) Max=    " << maxcoords_m << '\n'
       << "*                    Min=    " <<mincoords_m<<'\n'
       << "* Geometry length(m)         " <<len_m<<'\n'
       << "* Boundary box grid num      " <<nr_m<<'\n'
       << "* Boundary box size(m)       " <<hr_m<<'\n'
       << "* Size of boundary index set " << boundary_ids_m.size() <<'\n'
       << "* Number of all boxes        " << nr_m(0)*nr_m(1)*nr_m(2) << '\n'
       << "* Aligned Triangle Number    " << alignedT_m.size() <<endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}
