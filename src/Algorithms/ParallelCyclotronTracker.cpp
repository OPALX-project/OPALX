// ------------------------------------------------------------------------
// $RCSfile: ParallelCyclotronTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1 $initialLocalNum_m
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelCyclotronTracker
//   The class for tracking particles with 3D space charge in Cyclotrons and FFAG's
//
// ------------------------------------------------------------------------
//
// $Date: 2007/10/17 04:00:08 $
// $Author: adelmann, yang $
//
// ------------------------------------------------------------------------
#include <cfloat>
#include <iostream>
#include <fstream>
#include <vector>
#include "Algorithms/ParallelCyclotronTracker.h"

#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/Stripper.h"

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"

#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"

#include "Physics/Physics.h"

#include "Utilities/NumToStr.h"
#include "Utilities/OpalException.h"


#include "Ctunes.h"
#include "Ctunes.cc"
#include <cassert>


#include <hdf5.h>
#include "H5hut.h"

class Beamline;
class PartData;
using Physics::c;
using Physics::m_p; // GeV
using Physics::PMASS;
using Physics::PCHARGE;
using Physics::pi;
using Physics::q_e;

const double c_mmtns = c * 1.0e-6; // m/s --> mm/ns
const double mass_coeff = 1.0e18 * q_e / c / c; // from GeV/c^2 to basic unit: GV*C*s^2/m^2

#define PSdim 6

extern Inform *gmsg;

// typedef FVector<double, PSdim> Vector;

/**
 * Constructor ParallelCyclotronTracker
 *
 * @param beamline
 * @param reference
 * @param revBeam
 * @param revTrack
 */
ParallelCyclotronTracker::ParallelCyclotronTracker(const Beamline &beamline,
        const PartData &reference,
        bool revBeam, bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    sphys(NULL),
    myNode_m(Ippl::myNode()),
    initialLocalNum_m(0),
    initialTotalNum_m(0) {
    itsBeamline = dynamic_cast<Beamline *>(beamline.clone());
}

/**
 * Constructor ParallelCyclotronTracker
 *
 * @param beamline
 * @param bunch
 * @param ds
 * @param reference
 * @param revBeam
 * @param revTrack
 * @param maxSTEPS
 * @param timeIntegrator
 */
ParallelCyclotronTracker::ParallelCyclotronTracker(const Beamline &beamline,
                                                   PartBunch &bunch,
                                                   DataSink &ds,
                                                   const PartData &reference,
                                                   bool revBeam, bool revTrack,
                                                   int maxSTEPS, int timeIntegrator):
    Tracker(beamline, reference, revBeam, revTrack),
    sphys(NULL),
    maxSteps_m(maxSTEPS),
    timeIntegrator_m(timeIntegrator),
    myNode_m(Ippl::myNode()),
    initialLocalNum_m(bunch.getLocalNum()),
    initialTotalNum_m(bunch.getTotalNum()) {
    itsBeamline = dynamic_cast<Beamline *>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    //  scaleFactor_m = itsBunch->getdT() * c;
    scaleFactor_m = 1;
    multiBunchMode_m = 0;

    IntegrationTimer_m = IpplTimings::getTimer("Integration");
    TransformTimer_m   = IpplTimings::getTimer("Frametransform");
    DumpTimer_m        = IpplTimings::getTimer("Dump");
    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");
}

/**
 * Destructor ParallelCyclotronTracker
 *
 */
ParallelCyclotronTracker::~ParallelCyclotronTracker() {
    for(list<Component *>::iterator compindex = myElements.begin(); compindex != myElements.end(); compindex++) {
        delete(*compindex);
    }
    for(beamline_list::iterator fdindex = FieldDimensions.begin(); fdindex != FieldDimensions.end(); fdindex++) {
        delete(*fdindex);
    }
    delete itsBeamline;
}

/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::openFiles(string SfileName) {

    string  SfileName2 = SfileName + string("-Angle0.dat");

    outfTheta0_m.precision(8);
    outfTheta0_m.setf(ios::scientific, ios::floatfield);

    outfTheta0_m.open(SfileName2.c_str());
    outfTheta0_m << "#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
                 << endl;

    SfileName2 = SfileName + string("-Angle1.dat");

    outfTheta1_m.precision(8);
    outfTheta1_m.setf(ios::scientific, ios::floatfield);

    outfTheta1_m.open(SfileName2.c_str());
    outfTheta1_m << "#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
                 << endl;

    SfileName2 = SfileName + string("-Angle2.dat");

    outfTheta2_m.precision(8);
    outfTheta2_m.setf(ios::scientific, ios::floatfield);

    outfTheta2_m.open(SfileName2.c_str());
    outfTheta2_m << "#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
                 << endl;

    // for single Particle Mode, output after each turn, to define matched initial phase ellipse.

    SfileName2 = SfileName + string("-afterEachTurn.dat");

    outfThetaEachTurn_m.precision(8);
    outfThetaEachTurn_m.setf(ios::scientific, ios::floatfield);

    outfThetaEachTurn_m.open(SfileName2.c_str());
    outfThetaEachTurn_m << "# r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]" << endl;

}

/**
 * Close all files related to
 * special output in the Cyclotron
 * mode.
 */
void ParallelCyclotronTracker::closeFiles() {

    outfTheta0_m.close();
    outfTheta1_m.close();
    outfTheta2_m.close();
    outfThetaEachTurn_m.close();
}



/**
 *
 *
 * @param cycl
 */
void ParallelCyclotronTracker::visitCyclotron(const Cyclotron &cycl) {

    *gmsg << "* --------- Cyclotron ------------------------------" << endl;

    Cyclotron *elptr = dynamic_cast<Cyclotron *>(cycl.clone());
    myElements.push_back(elptr);

    double ri = elptr->getRinit();
    *gmsg << "* RINIT= " << ri << " [mm]" << endl;
    referenceR = ri;

    double pri = elptr->getPRinit();
    //msg << "PRINIT= " << pri << " [CU]" << endl;
    referencePr = pri;

    double phii = elptr->getPHIinit();
    *gmsg << "* PHIINIT= " << phii << " [deg]" << endl;
    referenceTheta = phii;
    if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron", "PHIINIT is out of [-180, 180)!");
    }

    referencePz = 0.0;
    referencePtot =  itsReference.getGamma() * itsReference.getBeta();
    referencePt = sqrt(referencePtot * referencePtot - referencePr * referencePr);
    if(referencePtot < 0.0) referencePt *= -1.0;

    sinRefTheta_m = sin(phii / 180.0 * pi);
    cosRefTheta_m = cos(phii / 180.0 * pi);

    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;

    *gmsg << "* Initial beta = " << itsReference.getBeta() << endl;

    *gmsg << "* Total reference momentum   = " << referencePtot * 1000.0 << " [MCU]" << endl;

    *gmsg << "* Reference azimuthal momentum  = " << referencePt * 1000.0 << " [MCU]" << endl;

    *gmsg << "* Reference radial momentum     = " << referencePr * 1000.0 << " [MCU]" << endl;

    double sym = elptr->getSymmetry();
    *gmsg << "* " << sym << " fold field symmerty " << endl;

    // ckr: this just returned the default value as defined in Component.h
    // double rff = elptr->getRfFrequ();
    // *gmsg << "* Rf frequency= " << rff << " [MHz]" << endl;

    string fmfn = elptr->getFieldMapFN();
    *gmsg << "* Field map file name= " << fmfn << " " << endl;

    string type = elptr->getType();
    *gmsg << "* Type of cyclotron= " << type << " " << endl;

    bool Sflag = elptr->getSuperpose();
    *gmsg << "* Electric field map are superpoesed ?  " << Sflag << " " << endl;


    double h = elptr->getCyclHarm();
    *gmsg << "* Harmonic number h= " << h << " " << endl;

    /**
     * To ease the initialise() function, set a integral parameter fieldflag internally.
     * Its value is  by the option "TYPE" of the element  "CYCLOTRON"
     * fieldflag = 1, readin PSI format measured field file (default)
     * fieldflag = 2, readin carbon cyclotron field file created by Jianjun Yang, TYPE=CARBONCYCL
     * fieldflag = 3, readin ANSYS format file for CYCIAE-100 created by Jianjun Yang, TYPE=CYCIAE
     * fieldflag = 4, readin AVFEQ format file for Riken cyclotrons
     * fieldflag = 5, readin FFAG format file for MSU/FNAL FFAG
     * fieldflag = 6, readin both median plane B field map and 3D E field map of RF cavity for compact cyclotron
     */
    int  fieldflag;
    if(type == string("CARBONCYCL")) {
        fieldflag = 2;
    } else if(type == string("CYCIAE")) {
        fieldflag = 3;
    } else if(type == string("AVFEQ")) {
        fieldflag = 4;
    } else if(type == string("FFAG")) {
        fieldflag = 5;
    } else if(type == string("BANDRF")) {
        fieldflag = 6;
    } else
        fieldflag = 1;

    // read field map on the  middle plane of cyclotron.
    // currently scalefactor is set to 1.0
    elptr->initialise(itsBunch, fieldflag, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++) BcParameter[i] = 0.0;
    string ElementType = "CYCLOTRON";
    BcParameter[0] = elptr->getRmin();
    BcParameter[1] = elptr->getRmax();

    // store inner radius and outer radius of cyclotron field map in the list
    buildupFieldList(BcParameter, ElementType, elptr);

}

/**
 * Not implemented and most probable never used
 *
 */
void ParallelCyclotronTracker::visitBeamBeam(const BeamBeam &) {
    *gmsg << "In BeamBeam tracker is missing " << endl;
}

/**
 *
 *
 * @param coll
 */
void ParallelCyclotronTracker::visitCollimator(const Collimator &coll) {

    Inform msg("visitCCollimator ");

    Collimator* elptr = dynamic_cast<Collimator *>(coll.clone());
    myElements.push_back(elptr);

    double angstart = elptr->getAngStart();
    msg << "AngStart= " << angstart << " [rad]" << endl;

    double angend = elptr->getAngEnd();
    msg << "AngEnd= " << angend << " [rad]" << endl;

    double rstart = elptr->getRStart();
    msg << "RStart= " << rstart << " [mm]" << endl;

    double rend = elptr->getREnd();
    msg << "REnd= " << rend << " [mm]" << endl;

    double width = elptr->getWidth();
    msg << "Width= " << width << " [mm]" << endl;

    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    string ElementType = "CCOLLIMATOR";
    BcParameter[0] = angstart ;
    BcParameter[1] = angend;
    BcParameter[2] = rstart ;
    BcParameter[3] = rend;
    BcParameter[4] = width ;
    buildupFieldList(BcParameter, ElementType, elptr);
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitCorrector(const Corrector &corr) {
    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Corrector *>(corr.clone()));
}

/**
 *
 *
 * @param diag
 */
void ParallelCyclotronTracker::visitDiagnostic(const Diagnostic &diag) {
    *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Diagnostic *>(diag.clone()));
}

/**
 *
 *
 * @param drift
 */
void ParallelCyclotronTracker::visitDrift(const Drift &drift) {
    *gmsg << "In drift L= " << drift.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Drift *>(drift.clone()));
}

/**
 *
 *
 * @param lamb
 */
void ParallelCyclotronTracker::visitLambertson(const Lambertson &lamb) {
    *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Lambertson *>(lamb.clone()));
}

/**
 *
 *
 * @param marker
 */
void ParallelCyclotronTracker::visitMarker(const Marker &marker) {
    //   *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Marker *>(marker.clone()));
    // Do nothing.
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitMonitor(const Monitor &corr) {
    //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Monitor *>(corr.clone()));
    //   applyDrift(flip_s * corr.getElementLength());
}

/**
 *
 *
 * @param mult
 */
void ParallelCyclotronTracker::visitMultipole(const Multipole &mult) {
    *gmsg << "In Multipole; L= " << mult.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Multipole *>(mult.clone()));
}

/**
 *
 *
 * @param prob
 */
void ParallelCyclotronTracker::visitProbe(const Probe &prob) {
    Inform msg("visitProbe ");

    Probe *elptr = dynamic_cast<Probe *>(prob.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXstart();
    msg << "XStart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXend();
    msg << "XEnd= " << xend << " [mm]" << endl;

    double ystart = elptr->getYstart();
    msg << "YStart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYend();
    msg << "YEnd= " << yend << " [mm]" << endl;

    double width = elptr->getWidth();
    msg << "Width= " << width << " [mm]" << endl;


    // initialise, do nothing
    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    string ElementType = "PROBE";
    BcParameter[0] = xstart ;
    BcParameter[1] = xend;
    BcParameter[2] = ystart ;
    BcParameter[3] = yend;
    BcParameter[4] = width ;

    // store probe parameters in the list
    buildupFieldList(BcParameter, ElementType, elptr);
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitRBend(const RBend &bend) {
    *gmsg << "In RBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RBend *>(bend.clone()));
}

/**
 *
 *
 * @param as
 */
void ParallelCyclotronTracker::visitRFCavity(const RFCavity &as) {

    *gmsg << "* --------- RFCavity ------------------------------" << endl;
    
    RFCavity *elptr = dynamic_cast<RFCavity *>(as.clone());
    myElements.push_back(elptr);

    if((elptr->getComponentType() != "SINGLEGAP") && (elptr->getComponentType() != "DOUBLEGAP")) {
        *gmsg << (elptr->getComponentType()) << endl;
        throw OpalException("ParallelCyclotronTracker::visitRFCavity",
                            "The ParallelCyclotronTracker can only play with cyclotron type RF system currently ...");
    }

    double rmin = elptr->getRmin();
    *gmsg << "* Minimal radius of cavity= " << rmin << " [mm]" << endl;

    double rmax = elptr->getRmax();
    *gmsg << "* Maximal radius of cavity= " << rmax << " [mm]" << endl;

    double rff = elptr->getCycFrequency();
    *gmsg << "* RF frequency (2*pi*f)= " << rff << " [rad/s]" << endl;

    string fmfn = elptr->getFieldMapFN();
    *gmsg << "* RF Field map file name= " << fmfn << endl;

    double angle = elptr->getAzimuth();
    *gmsg << "* Cavity azimuth position= " << angle << " [deg] " << endl;

    double gap = elptr->getGapWidth();
    *gmsg << "* Cavity gap width= " << gap << " [mm] " << endl;

    double pdis = elptr->getPerpenDistance();
    *gmsg << "* Cavity Shift distance= " << pdis << " [mm] " << endl;


    double phi0 = elptr->getPhi0();
    *gmsg << "* Initial RF phase (t=0)= " << phi0 << " [deg] " << endl;

    // read cavity voltage profile data from file.
    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    string ElementType = "CAVITY";
    BcParameter[0] = rmin;
    BcParameter[1] = rmax;
    BcParameter[2] = pdis;
    BcParameter[3] = angle;

    buildupFieldList(BcParameter, ElementType, elptr);
}

/**
 *
 *
 * @param rfq
 */
void ParallelCyclotronTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    *gmsg << "In RFQuadrupole; L= " << rfq.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RFQuadrupole *>(rfq.clone()));
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitSBend(const SBend &bend) {
    *gmsg << "In SBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<SBend *>(bend.clone()));
}

/**
 *
 *
 * @param sep
 */
void ParallelCyclotronTracker::visitSeparator(const Separator &sep) {
    *gmsg << "In Seapator L= " << sep.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Separator *>(sep.clone()));
}

/**
 *
 *
 * @param sept
 */
void ParallelCyclotronTracker::visitSeptum(const Septum &sept) {
    Inform msg("visitSeptum ");

    Septum *elptr = dynamic_cast<Septum *>(sept.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXstart();
    msg << "XStart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXend();
    msg << "XEnd= " << xend << " [mm]" << endl;

    double ystart = elptr->getYstart();
    msg << "YStart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYend();
    msg << "YEnd= " << yend << " [mm]" << endl;

    double width = elptr->getWidth();
    msg << "Width= " << width << " [mm]" << endl;


    // initialise, do nothing
    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    string ElementType = "SEPTUM";
    BcParameter[0] = xstart ;
    BcParameter[1] = xend;
    BcParameter[2] = ystart ;
    BcParameter[3] = yend;
    BcParameter[4] = width ;

    // store septum parameters in the list
    buildupFieldList(BcParameter, ElementType, elptr);
}

/**
 *
 *
 * @param solenoid
 */
void ParallelCyclotronTracker::visitSolenoid(const Solenoid &solenoid) {
    myElements.push_back(dynamic_cast<Solenoid *>(solenoid.clone()));
    Component *elptr = *(--myElements.end());
    if(!elptr->hasAttribute("ELEMEDGE")) {
        *gmsg << "Solenoid: no position of the element given!" << endl;
        return;
    }
}

/**
 *
 *
 * @param pplate
 */
void ParallelCyclotronTracker::visitParallelPlate(const ParallelPlate &pplate) {//do nothing

    //*gmsg << "ParallelPlate: not in use in ParallelCyclotronTracker!" << endl;

    //buildupFieldList(startField, endField, elptr);

}

/**
 *
 *
 * @param cv
 */
void ParallelCyclotronTracker::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}
/**
 * not used
 *
 * @param angle
 * @param curve
 * @param field
 * @param scale
 */
void ParallelCyclotronTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {

}

/**
 *
 *
 * @param stripper
 */

void ParallelCyclotronTracker::visitStripper(const Stripper &stripper) {

    *gmsg << "* ---------Stripper------------------------------" << endl;

    Stripper *elptr = dynamic_cast<Stripper *>(stripper.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXstart();
    *gmsg << "XStart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXend();
    *gmsg << "XEnd= " << xend << " [mm]" << endl;

    double ystart = elptr->getYstart();
    *gmsg << "YStart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYend();
    *gmsg << "YEnd= " << yend << " [mm]" << endl;

    double width = elptr->getWidth();
    *gmsg << "Width= " << width << " [mm]" << endl;

    double opcharge = elptr->getOPCharge();
    *gmsg << "Charge of outcome particle = +e * " << opcharge << endl;

    double opmass = elptr->getOPMass();
    *gmsg << "Mass of the outcome particle = " << opmass << " [GeV/c^2]" << endl;

    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    string ElementType = "STRIPPER";
    BcParameter[0] = xstart ;
    BcParameter[1] = xend;
    BcParameter[2] = ystart ;
    BcParameter[3] = yend;
    BcParameter[4] = width ;
    BcParameter[5] = opcharge;
    BcParameter[6] = opmass;

    buildupFieldList(BcParameter, ElementType, elptr);
}


void ParallelCyclotronTracker::applyExitFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {

}


/**
 *
 *
 * @param BcParameter
 * @param ElementType
 * @param elptr
 */
void ParallelCyclotronTracker::buildupFieldList(double BcParameter[], string ElementType, Component *elptr) {
    beamline_list::iterator sindex;

    type_pair *localpair = new type_pair();
    localpair->first = ElementType;

    for(int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;

    // always put cyclotron as the first element in the list.
    if(ElementType == "CYCLOTRON") {
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);

    // FixMe: so why not free up the memory?
    // free memory
    //  if(localpair)
    //    free(localpair);
}

/**
 *
 *
 * @param bl
 */
void ParallelCyclotronTracker::visitBeamline(const Beamline &bl) {
    itsBeamline->iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}


/**
 *
 *
 * @param s
 * @param nlp
 */
void checkNumPart(string s, int nlp) {
    int minnlp;
    int maxnlp;
    minnlp = 0;
    maxnlp = 111111;
    reduce(nlp, minnlp, OpMinAssign());
    reduce(nlp, maxnlp, OpMaxAssign());
    *gmsg << s << " min local particle number " << minnlp << " max local particle number: " << maxnlp << endl;
}

/**
 *
 *
 */
void ParallelCyclotronTracker::execute() {

    /*
      Initialize common variables and structures
      for the integrators
    */

    step_m = 0;
    restartStep0_m = 0;

    // record how many bunches has already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();

    // For the time being, we set bin number equal to bunch number. FixMe: not used
    BinCount_m = BunchCount_m;

    itsBeamline->accept(*this);

    // display the selected elements
    *gmsg << "-----------------------------" << endl;
    *gmsg << "The selected Beam line elements are :" << endl;
    for(beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++)
        *gmsg << ((*sindex)->first) << endl;
    *gmsg << "-----------------------------" << endl;

    // external field arrays for dumping
    for(int k = 0; k < 2; k++)
        FDext_m[k] = Vector_t(0.0, 0.0, 0.0);
    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    if(timeIntegrator_m == 0) {
        *gmsg << "* 4th order Runge-Kutta integrator" << endl;
        Tracker_RK4();
    } else if(timeIntegrator_m == 1) {
        *gmsg << "* 2nd order Leap-Frog integrator" << endl;
        Tracker_LF();
    } else if(timeIntegrator_m == 2) {
        *gmsg << "* Multiple time stepping (MTS) integrator" << endl;
        Tracker_MTS();
    } else {
        *gmsg << "ERROR: Invalid name of TIMEINTEGRATOR in Track command" << endl;
        exit(1);
    }

    *gmsg << "-----------------------------" << endl;
    *gmsg << "Finalize i.e. write data and close files :" << endl;
    for(beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++) {
        (((*sindex)->second).second)->finalise();
    }
    *gmsg << "-----------------------------" << endl;
}

/**
   In general the two tracker have much code in common.
   This is a great source of errors.
   Need to avoid this

*/



/**
 *
 *
 */
void ParallelCyclotronTracker::Tracker_LF() {

    Inform *gmsgAll;
    gmsgAll = new  Inform("CycTracker LF", INFORM_ALL_NODES);

    BorisPusher pusher;

    // time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();

    beamline_list::iterator sindex = FieldDimensions.begin();
    Cyclotron *elptr = dynamic_cast<Cyclotron *>(((*sindex)->second).second);
    if (elptr == NULL)         
        throw OpalException("ParallelCyclotronTracker::Tracker_LF()",
                            "The first item in the FieldDimensions list does not seem to be a cyclotron element");

    const double harm = elptr-> getCyclHarm();

    // load time
    const double dt = itsBunch->getdT() * 1.0e9 * harm; //[s]-->[ns]

    // find the injection time interval
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << stepsPerTurn *dt << "[ns]" << endl;
    }

    initTrackOrbitFile();

    // parameter for reset bin in multi-bunch run
    const  double eta = 0.01;

    int SteptoLastInj = itsBunch->getSteptoLastInj();

    // get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun() && numBunch_m > 1) {

        restartStep0_m = itsBunch->getTrackStep();
        step_m = restartStep0_m;
        itsBunch->resetPartBinID2(eta);
        *gmsg << "Restart at integration step " << restartStep0_m << endl;
    }

    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        itsBunch->setT(0.0);
    }

    *gmsg << "* Beginning of this run is at t= " << itsBunch->getT() * 1e9 << " [ns]" << endl;
    *gmsg << "* The time step is set to dt= " << dt << " [ns]" << endl;

    // for single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1)
        openFiles(inputFileNameWithoutExtension());

    double const initialReferenceTheta = referenceTheta / 180.0 * pi;
    if(!OpalData::getInstance()->inRestartRun()) {
        PathLength_m = 0.0;
        // Force the initial phase space values of the particle with ID=0 to zero, to set it as a reference particle.
        if(initialTotalNum_m > 2) {  // only for mulit-particle mode
            for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
                if(itsBunch->ID[ii] == 0) {
                    itsBunch->R[ii] = Vector_t(0.0);
                    itsBunch->P[ii] = Vector_t(0.0);
                }
            }
        }

        if(Options::psDumpLocalFrame) {
            // dump the initial distribution
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
            itsDataSink->writeStatData(*itsBunch, FDext_m, 0, 0, 0);
            *gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the local frame " << endl;
        }

        itsBunch->R *= Vector_t(1000.0); // m --> mm
        
        // Initialize global R
        Vector_t const initMeanR = Vector_t(referenceR * cosRefTheta_m, referenceR * sinRefTheta_m, 0.0);
        localToGlobal(itsBunch->R, initialReferenceTheta, initMeanR);
        
        // Initialize global P
        for(unsigned int ii = 0; ii < initialLocalNum_m; ++ii) {
            itsBunch->P[ii](0) += referencePr;
            itsBunch->P[ii](1) += referencePt;
        }
        localToGlobal(itsBunch->P, initialReferenceTheta);
        
        // Initialize the bin number of the first bunch to 0
        for(unsigned int ii = 0; ii < initialLocalNum_m; ++ii) {
            itsBunch->Bin[ii] = 0;
        }

        if(!(Options::psDumpLocalFrame)) {

            // dump the initial distribution
            itsBunch->R /= Vector_t(1000.0); // mm --> m
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);

            itsBunch->R *= Vector_t(1000.0); // m --> mm
            *gmsg << "meanR=( " << initMeanR(0) << " " << initMeanR(1) << " " << initMeanR(2) << " ) [mm] " << endl;
            *gmsg << "phase space dumping in glabol frame at initial position " << endl;
            *gmsg << "Dump step = " << lastDumpedStep_m << endl;
        }

        // AUTO mode
        if(multiBunchMode_m == 2) {
            Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m->mm
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;

            *gmsg << "initial Radial position = " << RThisTurn_m << " [mm]" << endl;
        }

        if((initialTotalNum_m > 2) && (numBunch_m > 1) && (multiBunchMode_m == 1)) {

            // backup initial distribution
	    initialR_m = new Vector_t[initialLocalNum_m];
            initialP_m = new Vector_t[initialLocalNum_m];

            for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
                initialR_m[ii] = itsBunch->R[ii];
                initialP_m[ii] = itsBunch->P[ii];
            }
        }

    } else {
        PathLength_m = itsBunch->getLPath();

        // AUTO mode
        if(multiBunchMode_m == 2) {
            itsBunch->calcBeamParameters_cycl();

            Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m->mm
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;

            *gmsg << "Radial position at restart position = " << RThisTurn_m << " [mm]" << endl;
        }

        itsBunch->R *= Vector_t(1000.0); // m --> mm

    }

    //  read in some control parameters
    const int SinglePartDumpFreq = Options::sptDumpFreq;
    const int resetBinFreq = Options::rebinFreq;
    const int scSolveFreq = Options::scSolveFreq;
    const bool doDumpAfterEachTurn = Options::psDumpEachTurn;


    int boundpDestroyFreq = 10; // todo: is it better treat as a control parameter

    // prepare for dump after each turn
    double oldReferenceTheta = initialReferenceTheta;

    *gmsg << "single particle trajectory dump frequency is set to " << SinglePartDumpFreq << endl;
    *gmsg << "particles repartition frequency is set to " << Options::repartFreq << endl;
    if(numBunch_m > 1)
        *gmsg << "particles energy bin ID reset frequency is set to " << resetBinFreq << endl;

    // if initialTotalNum_m = 2, trigger SEO mode
    // prepare for transverse tuning calculation
    vector<double> Ttime, Tdeltr, Tdeltz;
    // prepare for transverse tuning calculation
    vector<int> TturnNumber;
    int turnnumber = 1;


    // flag to determine when to transit from single-bunch to multi-bunches mode
    bool flagTransition = false;
    // step point determining the next time point of check for transition
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();

    const  double deltaTheta = pi / (stepsPerTurn);
    // record at which angle the space charge are solved
    double angleSpaceChargeSolve = 0.0;
    // local particle number
    int nlp;

    if(initialTotalNum_m == 1) {
        *gmsg << "* *---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," << endl
              << "* The initial distribution file must be specified which should contain only one line for the single particle " << endl
              << "* *------------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
        *gmsg << "* *------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." << endl
              << "* This mode does NOT include any RF cavities. The initial distribution file must be specified" << endl
              << "* In the file the first line is for reference particle and the second line is for offcenter particle." << endl
              << "* The tunes are calculated by FFT routines based on these two particles. " << endl
              << "* *------------NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    } else {

        itsBunch->boundp();
    }

    // do repartation for restart run
    if(OpalData::getInstance()->inRestartRun() && numBunch_m > 1) {
        Vector_t const meanR = calcMeanR();
        Vector_t const meanP = calcMeanP();
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
        globalToLocal(itsBunch->R, phi, meanR);
        
        //scale coordinates
        itsBunch->R /= Vector_t(1000.0); // mm --> m

        itsBunch->boundp();
        *gmsg << "Do binaryRepartation for restart run" << endl;
        nlp = itsBunch->getLocalNum();
        checkNumPart(string("Before repartation: "), nlp);
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        IpplTimings::stopTimer(BinRepartTimer_m);
        Ippl::Comm->barrier();
        nlp = itsBunch->getLocalNum();
        checkNumPart(string("After repartation: "), nlp);
        *gmsg << "Done binaryRepartation for restart run" << endl;

        //scale coordinates back
        itsBunch->R *= Vector_t(1000.0); // m --> mm

        localToGlobal(itsBunch->R, phi, meanR);
    }

    // *****************II***************
    // main integration loop
    // *****************II***************
    *gmsg << "---------------------------- Start tracking ----------------------------" << endl;
    for(; step_m < maxSteps_m; step_m++) {
        bool dumpEachTurn = false;
        if(step_m % SinglePartDumpFreq == 0) {
            singleParticleDump();
        }
        Ippl::Comm->barrier();

        // Push for first half step
        itsBunch->R *= Vector_t(0.001);
        push(0.5 * dt * 1e-9);
        itsBunch->R *= Vector_t(1000.0);

        // bunch injection
        if(numBunch_m > 1) {

            if((BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)) {
                if(step_m == stepsNextCheck) {
                    // under 3 conditions, following code will be execute
                    // to check the distance between two neighborring bunches
                    // 1.multi-bunch mode, AUTO sub-mode
                    // 2.After each revolution
                    // 3.only one bunch exists

                    *gmsg << "checking for automatically injecting new bunch ..." << endl;

                    itsBunch->R /= Vector_t(1000.0); // mm --> m
                    itsBunch->calcBeamParameters_cycl();
                    itsBunch->R *= Vector_t(1000.0); // m --> mm

                    Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm

                    RThisTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));

                    Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

                    double XYrms =  sqrt(pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0));


                    // if the distance between two neighbour bunch is less than CoeffDBunches_m times of its 2D rms size
                    // start multi-bunch simulation, fill current phase space to initialR and initialP arrays

                    if((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                        // since next turn, start multi-bunches
                        saveOneBunch();
                        flagTransition = true;

                        *gmsg << "*** Save beam distribution at turn #" << turnnumber << " ***" << endl;
                        *gmsg << "*** After one revolution, Multi-Bunch Mode will be invorked ***" << endl;

                    }

                    stepsNextCheck += stepsPerTurn;

                    *gmsg << "RLastTurn = " << RLastTurn_m << " [mm]" << endl;
                    *gmsg << "RThisTurn = " << RThisTurn_m << " [mm]" << endl;
                    *gmsg << "    XYrms = " << XYrms    << " [mm]" << endl;

                    RLastTurn_m = RThisTurn_m;
                }
            } else if(SteptoLastInj == stepsPerTurn - 1) {
                if(BunchCount_m < numBunch_m) {

                    // under 4 conditions, following code will be execute
                    // to read new bunch from hdf5 format file for FORCE or AUTO mode
                    // 1.multi-bunch mode
                    // 2.after each revolution
                    // 3.existing bunches is less than the specified bunches
                    // 4.FORCE mode, or AUTO mode with flagTransition = true
                    // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
                    *gmsg << "step " << step_m << ", inject a new bunch... ... ..." << endl;
                    BunchCount_m++;

                    // read initial distribution from h5 file
                    if(multiBunchMode_m == 1) {
                        readOneBunch(BunchCount_m - 1);
                        itsBunch->resetPartBinID2(eta);
                    } else if(multiBunchMode_m == 2) {

                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        itsBunch->resetPartBinID2(eta);
                    }

                    SteptoLastInj = 0;

                    itsBunch->setNumBunch(BunchCount_m);

                    stepsNextCheck += stepsPerTurn;

                    // update  after injection
                    itsBunch->boundp();

                    Ippl::Comm->barrier();
                    *gmsg << BunchCount_m << "'th bunch injected, total particle number = " << itsBunch->getTotalNum() << endl;
                }
            } else if(BunchCount_m == numBunch_m) {
                // After this, numBunch_m is wrong but useless
                numBunch_m--;

            } else {
                SteptoLastInj++;
            }
        }

        // calculate self fields Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            if(step_m % scSolveFreq == 0) {
                //    *gmsg << "Calculate space charge at step " << step_m<<endl;
                // Firstly reset E and B to zero before fill new space charge field data for each track step
                itsBunch->Bf = Vector_t(0.0);
                itsBunch->Ef = Vector_t(0.0);

                Vector_t const meanR = calcMeanR();
                if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
                    IpplTimings::startTimer(TransformTimer_m);
                    double const binsPhi = itsBunch->calcMeanPhi() - 0.5 * pi;
                    angleSpaceChargeSolve = binsPhi;
                    globalToLocal(itsBunch->R, binsPhi, meanR);

                    //scale coordinates
                    itsBunch->R /= Vector_t(1000.0); // mm --> m

                    if((step_m + 1) % boundpDestroyFreq == 0)
                        itsBunch->boundp_destroy();
                    else
                        itsBunch->boundp();

                    IpplTimings::stopTimer(TransformTimer_m);

                    // calcualte gamma for each energy bin
                    itsBunch->calcGammas_cycl();

                    repartition();
                    
                    // calculate space charge field for each energy bin
                    for(int b = 0; b < itsBunch->getLastemittedBin() ; b++) {

                        if(itsBunch->pbin_m->getTotalNumPerBin(b) >= 1000) {
                            //if(itsBunch->getNumPartInBin(b) >= 1000) {
                            itsBunch->setBinCharge(b, itsBunch->getChargePerParticle());
                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                            itsBunch->computeSelfFields_cycl(b);
                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                            INFOMSG("Bin:" << b << ", charge per particle " <<  itsBunch->getChargePerParticle() << endl);
                        } else {
                            INFOMSG("Note: Bin " << b << ": less than 1000 particles, omit space charge fields" << endl);
                        }
                    }

                    itsBunch->Q = itsBunch->getChargePerParticle();

                    IpplTimings::startTimer(TransformTimer_m);

                    //scale coordinates back
                    itsBunch->R *= Vector_t(1000.0); // m --> mm
                    
                    localToGlobal(itsBunch->R, binsPhi, meanR);
                    localToGlobal(itsBunch->Ef, binsPhi);
                    localToGlobal(itsBunch->Bf, binsPhi);
                } else {
                    Vector_t const meanP = calcMeanP();
                    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                    angleSpaceChargeSolve = phi;
                    globalToLocal(itsBunch->R, phi, meanR);

                    //scale coordinates
                    itsBunch->R /= Vector_t(1000.0); // mm --> m

                    if((step_m + 1) % boundpDestroyFreq == 0)
                        itsBunch->boundp_destroy();
                    else
                        itsBunch->boundp();

                    IpplTimings::stopTimer(TransformTimer_m);
                    repartition();
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                    double const meanGamma = sqrt(1.0 + pow(meanP(0), 2.0) + pow(meanP(1), 2.0));
                    itsBunch->computeSelfFields_cycl(meanGamma);
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

                    IpplTimings::startTimer(TransformTimer_m);

                    //scale coordinates back
                    itsBunch->R *= Vector_t(1000.0); // m --> mm
                    
                    localToGlobal(itsBunch->R, phi, meanR);
                    localToGlobal(itsBunch->Ef, phi);
                    localToGlobal(itsBunch->Bf, phi);
                }

                IpplTimings::stopTimer(TransformTimer_m);
            } else {
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                double const deltaPhi = phi - angleSpaceChargeSolve;
                localToGlobal(itsBunch->Ef, deltaPhi);
                localToGlobal(itsBunch->Bf, deltaPhi);
            }
        } else {
            // if field solver is not available , only update bunch, to transfer particles between nodes if needed,
            // reset parameters such as LocalNum, initialTotalNum_m.
            // INFOMSG("No space charge Effects are included!"<<endl;);
            if((step_m % Options::repartFreq * 100) == 0 && initialTotalNum_m >= 1000) {
                Vector_t const meanR = calcMeanR();
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                angleSpaceChargeSolve = phi; // we do not solve anything why set this?
                globalToLocal(itsBunch->R, phi, meanR);

                //scale coordinates
                itsBunch->R /= Vector_t(1000.0); // mm --> m

                if((step_m + 1) % boundpDestroyFreq == 0)
                    itsBunch->boundp_destroy();
                else
                    itsBunch->boundp();
                repartition();

                //scale coordinates back
                itsBunch->R *= Vector_t(1000.0); // m --> mm

                localToGlobal(itsBunch->R, phi, meanR);
            }
        }

        //  kick particles for one step
        IpplTimings::startTimer(IntegrationTimer_m);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            Vector_t externalE, externalB;
            double partR;

            externalB = Vector_t(0.0, 0.0, 0.0);
            externalE = Vector_t(0.0, 0.0, 0.0);

            partR = sqrt(dot(itsBunch->R[i], itsBunch->R[i]));

            beamline_list::iterator sindex = FieldDimensions.begin();

            if(((((*sindex)->second).first)[0] <= partR) && ((((*sindex)->second).first)[1] >= partR))
                (((*sindex)->second).second)->apply(itsBunch->R[i], Vector_t(0.0), itsBunch->getT() * 1e9, externalE, externalB);

	    externalB = externalB / 10.0; // kgauss -> T

            if(itsBunch->hasFieldSolver()) {
                externalE += itsBunch->Ef[i];
                externalB += itsBunch->Bf[i];
            }
            pusher.kick(itsBunch->R[i], itsBunch->P[i], externalE , externalB, dt * 1.0e-9, itsBunch->M[i] * 1.0e9, itsBunch->Q[i] / q_e);
        }
        IpplTimings::stopTimer(IntegrationTimer_m);

        // Push for second half step
        itsBunch->R *= Vector_t(0.001);
        push(0.5 * dt * 1e-9);
        itsBunch->R *= Vector_t(1000.0);
	
	// apply the plugin elements: probe, collimator, stripper, septum 
	bool flagNeedUpdate=false;
	flagNeedUpdate = applyPluginElements(turnnumber, dt);
	if(itsBunch->weHaveBins() && flagNeedUpdate)
	  itsBunch->resetPartBinID2(eta);
	
        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % resetBinFreq == 0)
            itsBunch->resetPartBinID2(eta);

        // dump  data after one push in single particle tracking
        if(initialTotalNum_m == 1) {
            int i = 0;

            // change phase space parameters from local reference frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
            for(int j = 0; j < 3; j++) {
                variable_m[j]   = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
            }

            double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]);//[ -pi ~ pi ]
            if((oldReferenceTheta < initialReferenceTheta - deltaTheta) &&
               (temp_meanTheta >= initialReferenceTheta - deltaTheta)) {
                turnnumber++;
                *gmsg << "Turn " << turnnumber << endl;
                dumpEachTurn = true;
                outfThetaEachTurn_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfThetaEachTurn_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                                    << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                                    << " " << temp_meanTheta / pi * 180
                                    << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                                    << " " << variable_m[2]
                                    << " " << variable_m[5] << endl;
            }
            // FixMe: should be defined elesewhere !
            // define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
            const double azimuth_angle0 = 0.0;
            const double azimuth_angle1 = 22.5 / 180.0 * pi;
            const double azimuth_angle2 = 45.0 / 180.0 * pi;
            if((oldReferenceTheta < azimuth_angle0 - deltaTheta) && (temp_meanTheta >= azimuth_angle0 - deltaTheta)) {
                outfTheta0_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                outfTheta1_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                outfTheta2_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }
            oldReferenceTheta = temp_meanTheta;
        }


        // check whether one turn over for multi-bunch tracking.
        if(doDumpAfterEachTurn && initialTotalNum_m > 2) {
            Vector_t const meanR = calcMeanR();

            // in global Cartesian frame, calculate the location in global frame of bunch
            oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

            // both for single bunch and multi-bunch
            // avoid dump at the first step
            // dumpEachTurn has not been changed in first push
            if((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {
                ++turnnumber;
                dumpEachTurn = true;
                *gmsg << "Turn " << turnnumber << " total particles " << itsBunch->getTotalNum() << endl;
            }
        }

        // dump phase space distribution of bunch
        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2) ||
           (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {

            IpplTimings::startTimer(DumpTimer_m);

            itsBunch->setSteptoLastInj(SteptoLastInj);

            itsBunch->setTrackStep((step_m + 1));

            extE_m = Vector_t(0.0, 0.0, 0.0);
            extB_m = Vector_t(0.0, 0.0, 0.0);

            //--------------------- calculate mean coordinates  of bunch -------------------------------//
            //------------  and calculate the external field at the mass of bunch-----------------------//

            Vector_t const meanR = calcMeanR();            
            *gmsg << "meanR=( " << meanR(0) << " " << meanR(1) << " " << meanR(2) << " ) [mm] " << endl;

            double meanRadius = sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1));

            beamline_list::iterator DumpSindex = FieldDimensions.begin();

            if(((((*DumpSindex)->second).first)[0] <= meanRadius) &&
               ((((*DumpSindex)->second).first)[1] >= meanRadius)) {
                (((*DumpSindex)->second).second)->apply(meanR, Vector_t(0.0), itsBunch->getT() * 1e9, extE_m, extB_m);
            }

            FDext_m[0] = extB_m / 10.0; // kgauss -> T
            FDext_m[1] = extE_m;

            //----------------------------dump in global frame-------------------------------------//
            // Note: Don't dump when
            // 1. after one turn
            // in order to sychronize the dump step for multi-bunch and single bunch for compare
            // with each other during post-process phase.
            if(!(Options::psDumpLocalFrame)) {
                itsBunch->R /= Vector_t(1000.0); // mm --> m

                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);

                //  itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);
                itsBunch->R *= Vector_t(1000.0); // m --> mm
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (global frame) at integration step "
                      << step_m + 1 << " T= " << itsBunch->getT() * 1e9 << " [ns]" << endl;

                //----------------------------dump in local frame-------------------------------------//
            } else {
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                globalToLocal(itsBunch->R, phi, meanR);
                globalToLocal(itsBunch->P, phi, meanP);
                itsBunch->R /= Vector_t(1000.0); // mm --> m
                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
                itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0);
                itsBunch->R *= Vector_t(1000.0); // m --> mm
                localToGlobal(itsBunch->R, phi, meanR);
                localToGlobal(itsBunch->P, phi, meanP);
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (local frame) at integration step "
                      << step_m + 1 << " T= " << itsBunch->getT() * 1e9 << " [ns]" << endl;

            }
            IpplTimings::stopTimer(DumpTimer_m);
        }
    }

    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 0) {
            // FixMe: FinalMomentum2  = dot(itsBunch->P[ii],itsBunch->P[ii]);
            double FinalMomentum2  = pow(itsBunch->P[ii](0), 2.0) +
                                     pow(itsBunch->P[ii](1), 2.0) +
                                     pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsgAll << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsgAll << "* Total phase space dump number including the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsgAll << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
        }
    }

    Ippl::Comm->barrier();

    if(myNode_m == 0) outfTrackOrbit_m.close();

    if(initialTotalNum_m == 1)
        closeFiles();

    *gmsg << *itsBunch << endl;

    // free memory
    if(gmsgAll)
        free(gmsgAll);

}

void ParallelCyclotronTracker::Tracker_RK4() {

    Inform *gmsgAll;
    gmsgAll = new  Inform("CycTracker RK4", INFORM_ALL_NODES);
    // time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();
    // record how many bunches has already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();
    // decide how many energy bins. ONLY FOR MPM
    BinCount_m = BunchCount_m;
 
    beamline_list::iterator sindex = FieldDimensions.begin();
    Cyclotron *elptr = dynamic_cast<Cyclotron *>(((*sindex)->second).second);
    if (elptr == NULL)         
        throw OpalException("ParallelCyclotronTracker::Tracker_LF()",
                            "The first item in the FieldDimensions list does not seem to be a cyclotron element");
    const double harm = elptr-> getCyclHarm();

    // load time
    double t  = itsBunch->getT() * 1.0e9;
    const double dt = itsBunch->getdT() * 1.0e9 * harm; //[s]-->[ns]

    // find the injection time interval
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << stepsPerTurn *dt << "[ns]" << endl;
    }

    initTrackOrbitFile();
    // get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun()) {
        restartStep0_m = itsBunch->getTrackStep();
        step_m = restartStep0_m;
        *gmsg << "Restart at integration step " << restartStep0_m << endl;
    }

    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        t = 0.0;
    }

    *gmsg << "* Beginning of this run is at t= " << t << " [ns]" << endl;
    *gmsg << "* The time step is set to dt= " << dt << " [ns]" << endl;

    // for single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1)
        openFiles(inputFileNameWithoutExtension());

    if(!OpalData::getInstance()->inRestartRun()) {
        PathLength_m = 0.0;
        // Force the initial phase space values of the particle with ID=0 to zero, to set it as a reference particle.
        if(initialTotalNum_m > 2) {  // only for mulit-particle mode
            for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
                if(itsBunch->ID[ii] == 0) {
                    itsBunch->R[ii] = Vector_t(0.0);
                    itsBunch->P[ii] = Vector_t(0.0);
                }
            }
        }

        if(Options::psDumpLocalFrame) {
            // dump the initial distribution
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
            itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0);
            *gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the local frame " << endl;
        }

        itsBunch->R *= Vector_t(1000.0); // m --> mm
        for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
            // rotate the local beam frame for angle of minus theta
            double tempDx = itsBunch->R[ii](0) * cosRefTheta_m - itsBunch->R[ii](1) * sinRefTheta_m;
            double tempDy = itsBunch->R[ii](0) * sinRefTheta_m + itsBunch->R[ii](1) * cosRefTheta_m;

            // move the origin point of global frame
            itsBunch->R[ii](0) = tempDx + referenceR * cosRefTheta_m;
            itsBunch->R[ii](1) = tempDy + referenceR * sinRefTheta_m;

            double tempPx = (itsBunch->P[ii](0) + referencePr) * cosRefTheta_m
                            - (itsBunch->P[ii](1) + referencePt) * sinRefTheta_m;
            double tempPy = (itsBunch->P[ii](0) + referencePr) * sinRefTheta_m
                            + (itsBunch->P[ii](1) + referencePt) * cosRefTheta_m;

            itsBunch->P[ii](0) = tempPx;
            itsBunch->P[ii](1) = tempPy;
            // initialize the bin number of the first bunch to 0
            itsBunch->Bin[ii] = 0;
        }
        if(!(Options::psDumpLocalFrame)) {
            // dump the initial distribution
            itsBunch->R /= Vector_t(1000.0); // mm --> m
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
            itsBunch->R *= Vector_t(1000.0); // m --> mm
            *gmsg << "* meanR=( " << referenceR *cosRefTheta_m << " " << referenceR *sinRefTheta_m << " " << 0.0 << " ) [mm] " << endl;
            *gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the global frame " << endl;
        }

        // AUTO mode
        if(multiBunchMode_m == 2) {
            Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m->mm
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;
            *gmsg << "* Initial Radial position = " << RThisTurn_m << " [mm]" << endl;
        }

        if((initialTotalNum_m > 2) && (numBunch_m > 1) && (multiBunchMode_m == 1)) {
            // backup initial distribution
            initialR_m = new Vector_t[initialLocalNum_m];
            initialP_m = new Vector_t[initialLocalNum_m];

            for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
                initialR_m[ii] = itsBunch->R[ii];
                initialP_m[ii] = itsBunch->P[ii];
            }
        }
    } else {
        PathLength_m = itsBunch->getLPath();
        // AUTO mode
        if(multiBunchMode_m == 2) {
            itsBunch->calcBeamParameters_cycl();
            Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m->mm
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;
            *gmsg << "* Radial position at restart position = " << RThisTurn_m << " [mm]" << endl;
        }
        itsBunch->R *= Vector_t(1000.0); // m --> mm
    }

    //  read in some control parameters
    const int SinglePartDumpFreq = Options::sptDumpFreq;
    const int resetBinFreq = Options::rebinFreq;
    const int scSolveFreq = Options::scSolveFreq;
    const bool doDumpAfterEachTurn = Options::psDumpEachTurn;

    // parameter for reset bin in multi-bunch run, todo: readin from inputfile
    const  double eta = 0.01;

    int boundpDestroyFreq = 10; // todo: is it better treat as a control parameter

    // prepare for dump after each turn
    const double initialReferenceTheta = referenceTheta / 180.0 * pi;
    double oldReferenceTheta = initialReferenceTheta;

    *gmsg << "* Single particle trajectory dump frequency is set to " << SinglePartDumpFreq << endl;
    *gmsg << "* The frequency to solve space charge fields is set to " << scSolveFreq << endl;
    *gmsg << "* The repartition frequency is set to " << Options::repartFreq << endl;

    if(numBunch_m > 1)
        *gmsg << "* The particles energy bin reset frequency is set to " << resetBinFreq << endl;

    // if initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
    vector<double> Ttime, Tdeltr, Tdeltz;
    vector<int> TturnNumber;
    int turnnumber = 1;
    int lastTurn = 1;

    bool flagNoDeletion = false;

    // flag to determine when to transit from single-bunch to multi-bunches mode
    bool flagTransition = false;
    // step point determining the next time point of check for transition
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();

    const double deltaTheta = pi / (stepsPerTurn); // half of the average angle per step
    // record at which angle the space charge are solved
    double angleSpaceChargeSolve = 0.0;
    // local particle number
    int nlp;

    if(initialTotalNum_m == 1) {
        *gmsg << "* ---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," << endl
              << "* The initial distribution file must be specified which should contain only one line for the single particle " << endl
              << "* ------------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
        *gmsg << "* ------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." << endl
              << "* This mode does NOT include any RF cavities. The initial distribution file must be specified" << endl
              << "* In the file the first line is for reference particle and the second line is for offcenter particle." << endl
              << "* The tune is calculated by FFT routines based on these two particles. " << endl
              << "* ------------NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    } else {
        nlp = itsBunch->getLocalNum();
        checkNumPart(string("Before boundp "), nlp);

        itsBunch->boundp();

        nlp = itsBunch->getLocalNum();
        checkNumPart(string("After boundp"), nlp);

    }

    // main integration loop
    *gmsg << "---------------------------- Start tracking ----------------------------" << endl;
    for(; step_m < maxSteps_m; step_m++) {
        bool dumpEachTurn = false;
        if(initialTotalNum_m > 2) {

            // single particle dumping
            if(step_m % SinglePartDumpFreq == 0) { // dump
                IpplTimings::startTimer(DumpTimer_m);

                double x;
                int  id;
                vector<double> tmpr;
                vector<int> tmpi;

                int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);

                    // for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
                    int found[2] = { -1, -1};
                    int counter = 0;

                    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        if(itsBunch->ID[ii] == 0) {
                            found[counter] = ii;
                            counter++;
                        }
                        if(itsBunch->ID[ii] == 1) {
                            found[counter] = ii;
                            counter++;
                        }
                    }
                    // for the regular modes only the space data of particles with ID = 0 and 1 need be transfored
                    if(myNode_m == 0) {
                        // for root node
                        int notReceived =  Ippl::getNodes() - 1;
                        int numberOfPart = 0;

                        while(notReceived > 0) {
                            int node = COMM_ANY_NODE;
                            Message *rmsg =  Ippl::Comm->receive_block(node, tag);
                            if(rmsg == 0)
                                ERRORMSG("Could not receive from client nodes in main." << endl);
                            notReceived--;
                            rmsg->get(&numberOfPart);
                            for(int ii = 0; ii < numberOfPart; ii++) {
                                rmsg->get(&id);
                                tmpi.push_back(id);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                                rmsg->get(&x);
                                tmpr.push_back(x);
                            }
                            delete rmsg;

                        }
                        for(int ii = 0; ii < counter; ii++) {
                            tmpi.push_back(itsBunch->ID[found[ii]]);
                            for(int jj = 0; jj < 3; jj++) {
                                tmpr.push_back(itsBunch->R[found[ii]](jj));
                                tmpr.push_back(itsBunch->P[found[ii]](jj));
                            }
                        }
                        vector<double>::iterator itParameter = tmpr.begin();
                        vector<int>::iterator  itId = tmpi.begin();

                        for(itId = tmpi.begin(); itId != tmpi.end(); itId++) {
                            outfTrackOrbit_m << "ID" << *itId;
                            for(int ii = 0; ii < 6; ii++) {
                                outfTrackOrbit_m << " " << *itParameter;
                                itParameter++;
                            }
                            outfTrackOrbit_m << endl;
                        }
                        // sample frequency = SinglePartDumpFreq
                    } else {
                        // for other nodes
                        Message *smsg = new Message();
                        smsg->put(counter);
                        for(int ii = 0; ii < counter; ii++) {
                            smsg->put(itsBunch->ID[found[ii]]);
                            for(int jj = 0; jj < 3; jj++) {
                                smsg->put(itsBunch->R[found[ii]](jj));
                                smsg->put(itsBunch->P[found[ii]](jj));
                            }
                        }
                        bool res = Ippl::Comm->send(smsg, 0, tag);
                        if(!res)
                            ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
                        // FixMe: why it block at here if I delete smsg? delete smsg;

                    }
                 
                IpplTimings::stopTimer(DumpTimer_m);
            }
            //end dump
            Ippl::Comm->barrier();

            // bunch injection
            if(numBunch_m > 1) {
                if((BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)) {
                    if(step_m == stepsNextCheck) {
                        // under 3 conditions, following code will be execute
                        // to check the distance between two neighborring bunches
                        // 1.multi-bunch mode, AUTO sub-mode
                        // 2.After each revolution
                        // 3.only one bunch exists

                        *gmsg << "checking for automatically injecting new bunch ..." << endl;

                        itsBunch->R /= Vector_t(1000.0); // mm --> m
                        itsBunch->calcBeamParameters_cycl();
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm

                        RThisTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));

                        Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

                        double XYrms =  sqrt(pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0));

                        // if the distance between two nieghbour bunch is less than 5 times of its 2D rms size
                        // start multi-bunch simulation, fill current phase space to initialR and initialP arrays

                        if((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                            // since next turn, start multi-bunches
                            saveOneBunch();
                            flagTransition = true;
                            *gmsg << "*** Save beam distribution at turn #" << turnnumber << " ***" << endl;
                            *gmsg << "*** After one revolution, Multi-Bunch Mode will be invorked ***" << endl;

                        }
                        stepsNextCheck += stepsPerTurn;

                        *gmsg << "RLastTurn = " << RLastTurn_m << " [mm]" << endl;
                        *gmsg << "RThisTurn = " << RThisTurn_m << " [mm]" << endl;
                        *gmsg << "    XYrms = " << XYrms    << " [mm]" << endl;

                        RLastTurn_m = RThisTurn_m;
                    }

                } else if((BunchCount_m < numBunch_m) && (step_m == stepsNextCheck)) {

                    // under 4 conditions, following code will be execute
                    // to read new bunch from hdf5 format file for FORCE or AUTO mode
                    // 1.multi-bunch mode
                    // 2.after each revolution
                    // 3.existing bunches is less than the specified bunches
                    // 4.FORCE mode, or AUTO mode with flagTransition = true
                    // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
                    *gmsg << "step " << step_m << ", inject a new bunch... ... ..." << endl;

                    BunchCount_m++;

                    // read initial distribution from h5 file
                    if(multiBunchMode_m == 1)
                        readOneBunch(BunchCount_m - 1);
                    else if(multiBunchMode_m == 2) {

                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        //itsBunch->resetPartBinID2(eta);
                    }
                    itsBunch->setNumBunch(BunchCount_m);

                    stepsNextCheck += stepsPerTurn;

                    // update  after injection
                    itsBunch->boundp();

                    Ippl::Comm->barrier();
                    *gmsg << BunchCount_m << "th bunch injected, total particle number = " << itsBunch->getTotalNum() << endl;

                } else if(BunchCount_m == numBunch_m) {
                    // After this, numBunch_m is wrong but useless
                    numBunch_m--;
                }
            }

            // Calculate SC field before each time step and keep constant during integration.
            // Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
            if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
                if(step_m % scSolveFreq == 0) {
                    // Firstly reset E and B to zero before fill new space charge field data for each track step
                    itsBunch->Bf = Vector_t(0.0);
                    itsBunch->Ef = Vector_t(0.0);

                    IpplTimings::startTimer(TransformTimer_m);

                    //HERE transform particles coordinates to local frame (rotate and shift)
                    Vector_t const meanR = calcMeanR();

                    // in global Cartesian frame, calculate the location in global frame of bunch
                    oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

                    if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
                        double binsMeanPhi = itsBunch->calcMeanPhi() - 0.5 * pi;
                        angleSpaceChargeSolve = binsMeanPhi;

                        double cosTemp_binsMeanPhi = cos(binsMeanPhi);
                        double sinTemp_binsMeanPhi = sin(binsMeanPhi);


                        // remove mean coordinates
                        itsBunch->R -= meanR;

                        //scale coordinates
                        itsBunch->R /= Vector_t(1000.0); // mm --> m

                        // rotate from global frame to local frame(transverse horizontal,longitudinal,transverse vertical)
                        // For multi-bin, rotate the frame for binsMeanPhi degree
                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_binsMeanPhi
                                                         + itsBunch->R[ii](1) * sinTemp_binsMeanPhi;
                            double temp_RLongitudinal = -itsBunch->R[ii](0) * sinTemp_binsMeanPhi
                                                        + itsBunch->R[ii](1) * cosTemp_binsMeanPhi;

                            itsBunch->R[ii](0) = temp_RHorizontal;
                            itsBunch->R[ii](1) = temp_RLongitudinal;
                        }

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        // FixMe: check this
                        //   IpplTimings::startTimer(BinRepartTimer_m);

                        // todo: overuse binaryRepart() can cause buffer overflow during large scale job.
                        // as a temporal solution, reduce its using frequency can solve such kind of problem.
                        //   if( (step_m+1)%Options::repartFreq == 0 )
                        //    {
                        //        NDIndex<3> ldom = itsBunch->getFieldLayout().getLocalNDIndex();
                        //        double hx, hy, hz;
                        //        int IndexMaxX,IndexMaxY,IndexMaxZ;
                        //        int IndexMinX,IndexMinY,IndexMinZ;
                        //        hx = itsBunch->get_hr()(0);
                        //        hy = itsBunch->get_hr()(1);
                        //        hz = itsBunch->get_hr()(2);
                        //        IndexMaxX =  ldom[0].max();
                        //        IndexMinX =  ldom[0].min();
                        //        IndexMaxY =  ldom[1].max();
                        //        IndexMinY =  ldom[1].min();
                        //        IndexMaxZ =  ldom[2].max();
                        //        IndexMinZ =  ldom[2].min();
                        //        *gmsg<<"hx = "<<hx<<", hy = "<<hy<<", hz = "<<hz<<endl;
                        //        Ippl::Comm->barrier();
                        //        *gmsgAll <<"min/max X:"<< IndexMinX<<"--"<<IndexMaxX
                        //                 << ",Y: "<< IndexMinY<<"--"<<IndexMaxY
                        //                 << ", Z: "<< IndexMinZ<<"--"<<IndexMaxZ<<endl ;
                        //        Ippl::Comm->barrier();
                        //        *gmsgAll <<"After boundp, particle  "<<itsBunch->getLocalNum()<<endl ;
                        //        Ippl::Comm->barrier();
                        //    }
                        //           IpplTimings::stopTimer(BinRepartTimer_m);
                        //           Ippl::Comm->barrier();


                        // calcualte gamma for each energy bin
                        itsBunch->calcGammas_cycl();

                        repartition();
                        
                        // calculate space charge field for each energy bin
                        for(int b = 0; b < itsBunch->getLastemittedBin(); b++) {

                            itsBunch->setBinCharge(b, itsBunch->getChargePerParticle());
                            itsBunch->computeSelfFields_cycl(b);
                        }

                        itsBunch->Q = itsBunch->getChargePerParticle();

                        IpplTimings::startTimer(TransformTimer_m);

                        // HERE transform particles coordinates back to global frame (rotate and shift)
                        // rotate back from local reference frame to global frame
                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_binsMeanPhi
                                                         - itsBunch->R[ii](1) * sinTemp_binsMeanPhi;
                            double temp_RLongitudinal =  itsBunch->R[ii](0) * sinTemp_binsMeanPhi
                                                         + itsBunch->R[ii](1) * cosTemp_binsMeanPhi;

                            itsBunch->R[ii](0) = temp_RHorizontal;
                            itsBunch->R[ii](1) = temp_RLongitudinal;
                        }

                        //scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // retrieve mean coordinates
                        itsBunch->R += meanR;

                        // HERE transform self field back to global frame (rotate)

                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_binsMeanPhi
                                              - itsBunch->Ef[ii](1) * sinTemp_binsMeanPhi;
                            double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_binsMeanPhi
                                              + itsBunch->Ef[ii](1) * cosTemp_binsMeanPhi;

                            itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                            itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                        }

                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_binsMeanPhi
                                              - itsBunch->Bf[ii](1) * sinTemp_binsMeanPhi;
                            double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_binsMeanPhi
                                              + itsBunch->Bf[ii](1) * cosTemp_binsMeanPhi;

                            itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                            itsBunch->Bf[ii](1) = temp_E2;  // By,T
                        }
                    } else {
                        Vector_t const meanP = calcMeanP();

                        // in global Cartesian frame, calculate the direction of longitudinal angle of bunch
                        double meanPhi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

                        angleSpaceChargeSolve = meanPhi;

                        double cosTemp_meanPhi = cos(meanPhi);
                        double sinTemp_meanPhi = sin(meanPhi);

                        double meanPLongitudinal2 = pow(meanP(0), 2.0) + pow(meanP(1), 2.0);
                        double temp_meangamma = sqrt(1.0 + meanPLongitudinal2);

                        // remove mean coordinates
                        itsBunch->R -= meanR;

                        //scale coordinates
                        itsBunch->R /= Vector_t(1000.0); // mm --> m

                        // rotate from global frame to local frame(transverse horizontal,longitudinal, transverse vertical)
                        // For single bin, rotate the frame for meanPhi degree
                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_meanPhi
                                                         + itsBunch->R[ii](1) * sinTemp_meanPhi;
                            double temp_RLongitudinal = -itsBunch->R[ii](0) * sinTemp_meanPhi
                                                        + itsBunch->R[ii](1) * cosTemp_meanPhi;

                            itsBunch->R[ii](0) = temp_RHorizontal;
                            itsBunch->R[ii](1) = temp_RLongitudinal;
                        }

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        // FixMe:
                        //    IpplTimings::startTimer(BinRepartTimer_m);

                        // todo: overuse binaryRepart() can cause baffer overflow during large scale job.
                        // as a temporal solution, reduce its using frequency can solve such kind of problem.

                        //  if( (step_m+1)%Options::repartFreq == 0 )
                        //    {
                        //        NDIndex<3> ldom = itsBunch->getFieldLayout().getLocalNDIndex();
                        //        double hx, hy, hz;
                        //        int IndexMaxX,IndexMaxY,IndexMaxZ;
                        //        int IndexMinX,IndexMinY,IndexMinZ;
                        //        hy = itsBunch->get_hr()(1);
                        //        hz = itsBunch->get_hr()(2);
                        //        IndexMaxX =  ldom[0].max();
                        //        IndexMinX =  ldom[0].min();
                        //        IndexMaxY =  ldom[1].max();
                        //        IndexMinY =  ldom[1].min();
                        //        IndexMaxZ =  ldom[2].max();
                        //        IndexMinZ =  ldom[2].min();
                        //        *gmsg<<"hx = "<<hx<<", hy = "<<hy<<", hz = "<<hz<<endl;
                        //        Ippl::Comm->barrier();
                        //        *gmsgAll<<"min/max X:"<< IndexMinX<<"--"<<IndexMaxX<<", Y: "
                        //                << IndexMinY<<"--"<<IndexMaxY<<", Z: "<< IndexMinZ<<"--"<<IndexMaxZ<<endl ;
                        //        Ippl::Comm->barrier();
                        //        *gmsgAll <<"After boundp, particle  "<<itsBunch->getLocalNum()<<endl ;
                        //        Ippl::Comm->barrier();
                        //        *gmsg << "***---------------------------- DO REPART   ----------------------*** "<<endl;
                        //        itsBunch->do_binaryRepart();
                        //        *gmsg << "***---------------------------- REPART DONE ----------------------*** "<<endl;
                        //    }
                        //        IpplTimings::stopTimer(BinRepartTimer_m);
                        //           Ippl::Comm->barrier();

                        repartition();
                        itsBunch->computeSelfFields_cycl(temp_meangamma);

                        IpplTimings::startTimer(TransformTimer_m);
                        // HERE transform particles coordinates back to global frame (rotate and shift)
                        // rotate back from local reference frame to global frame
                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_meanPhi
                                                         - itsBunch->R[ii](1) * sinTemp_meanPhi;
                            double temp_RLongitudinal =  itsBunch->R[ii](0) * sinTemp_meanPhi
                                                         + itsBunch->R[ii](1) * cosTemp_meanPhi;

                            itsBunch->R[ii](0) = temp_RHorizontal;
                            itsBunch->R[ii](1) = temp_RLongitudinal;
                        }

                        //scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // retrieve mean coordinates
                        itsBunch->R += meanR;

                        // HERE transform self field back to global frame (rotate)
                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_meanPhi  - itsBunch->Ef[ii](1) * sinTemp_meanPhi;
                            double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_meanPhi  + itsBunch->Ef[ii](1) * cosTemp_meanPhi;
                            itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                            itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                        }

                        for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                            double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_meanPhi  - itsBunch->Bf[ii](1) * sinTemp_meanPhi;
                            double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_meanPhi  + itsBunch->Bf[ii](1) * cosTemp_meanPhi;
                            itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                            itsBunch->Bf[ii](1) = temp_E2;  // By,T
                        }
                    }

                    IpplTimings::stopTimer(TransformTimer_m);

                } else {
                    //HERE transform particles coordinates to local frame (rotate and shift)
                    Vector_t const meanR = calcMeanR();
                    Vector_t const meanP = calcMeanP();

                    // in global Cartesian frame, calculate the location in global frame of bunch
                    oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

                    // in global Cartesian frame, calculate the direction of longitudinal angle of bunch
                    double meanPhi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                    
                    double deltaPhi = meanPhi - angleSpaceChargeSolve;

                    double cosTemp_deltaPhi = cos(deltaPhi);
                    double sinTemp_deltaPhi = sin(deltaPhi);


                    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_deltaPhi  - itsBunch->Ef[ii](1) * sinTemp_deltaPhi;
                        double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_deltaPhi  + itsBunch->Ef[ii](1) * cosTemp_deltaPhi;
                        itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                        itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                    }

                    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_deltaPhi  - itsBunch->Bf[ii](1) * sinTemp_deltaPhi;
                        double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_deltaPhi  + itsBunch->Bf[ii](1) * cosTemp_deltaPhi;
                        itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                        itsBunch->Bf[ii](1) = temp_E2;  // By,T
                    }
                }

            } else {
                // if field solver is not available , only update bunch, to transfer particles between nodes if needed,
                // reset parameters such as LocalNum, initialTotalNum_m.

                // in global Cartesian frame, calculate the location in global frame of bunch
                Vector_t const meanR = calcMeanR();
                oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

                // if no space charge effects are included, don't need to call update() in local frame
            }
            // track all particles for one step
            IpplTimings::startTimer(IntegrationTimer_m);


            for(size_t i = 0; i < (itsBunch->getLocalNum()); i++) {
                flagNoDeletion = true;
                // change phase space parameters from localframe of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
                for(int j = 0; j < 3; j++) {
                    variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                    variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
                    rold_m[j] = variable_m[j]; // used for gap cross checking
                    pold_m[j] = variable_m[j+3]; // used for gap cross
                }

                // integrate for one step in the lab Cartesian frame (absulate value ).
                // IpplTimings::startTimer(IntegrationTimer_m);
                flagNoDeletion = rk4(variable_m, t, dt, i);

                // IpplTimings::stopTimer(IntegrationTimer_m);
                for(int j = 0; j < 3; j++) {
                    itsBunch->R[i](j) = variable_m[j];    //[x,y,z]  units: [mm]
                    itsBunch->P[i](j) = variable_m[j+3];  //[px,py,pz]  units: dimensionless, beta*gama
                }

                //If gap crossing happens, do momenta kicking
                for(beamline_list::iterator sindex = ++(FieldDimensions.begin());
                    sindex != FieldDimensions.end(); sindex++) {

                    bool tag_crossing = false;
                    double DistOld = 0.0; //mm
                    RFCavity * rfcav;
                    if(((*sindex)->first) == "CAVITY") {
                        // here check gap cross in the list, if do , set tag_crossing to TRUE
                        for(int j = 0; j < 3; j++)
                            rnew_m[j] = variable_m[j];
                        rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
                        tag_crossing = checkGapCross(rold_m, rnew_m, rfcav, DistOld);
                    }
                    if(tag_crossing) {
                        double oldMomentum2  = dot(pold_m, pold_m);
                        double oldBetgam = sqrt(oldMomentum2);
                        double oldGamma = sqrt(1.0 + oldMomentum2);
                        double oldBeta = oldBetgam / oldGamma;
                        double dt1 = DistOld / (c * oldBeta * 1.0e-6); // ns
                        double dt2 = dt - dt1;
			
                        // retrack particle from the old postion to cavity gap point
                        // restore the old coordinates and momenta
                        for(int j = 0; j < 3; j++) {
                            variable_m[j] = rold_m[j];
                            variable_m[j+3] = pold_m[j];
                        }

                        if(dt / dt1 < 1.0e9) rk4(variable_m, t, dt1, i);

                        for(int j = 0; j < 3; j++) {
                            itsBunch->R[i](j) = variable_m[j] ;    //[x,y,z]  units: [mm]
                            itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units:[] beta*gama
                        }

                        //momentum kick
                        RFkick(rfcav, t, dt1, i);

                        // retrack particle  from cavity gap point for the left time to finish the entire timestep
                        for(int j = 0; j < 3; j++) {
                            variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                            variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
                        }

                        if(dt / dt2 < 1.0e9) rk4(variable_m, t, dt2, i);

                        for(int j = 0; j < 3; j++) {
                            itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
                            itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [] beta*gama
                        }
                    } // end if: gap-crossing monentum kicking at certain cavity
                } //end for: finish checking for all cavities
            } //end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode
	    
	    // apply the plugin elements: probe, collimator, stripper, septum 
	    bool flagNeedUpdate=false;
	    flagNeedUpdate = applyPluginElements(turnnumber, dt);
	    if(itsBunch->weHaveBins() && flagNeedUpdate)
	      itsBunch->resetPartBinID2(eta);
	    
	    // recalculate bingamma and reset the BinID for each particles according to its current gamma
	    if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % resetBinFreq == 0)
	      itsBunch->resetPartBinID2(eta);

            if((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {
                turnnumber++;
                dumpEachTurn = true;
                *gmsg << "Turn " << turnnumber << " total particles " << itsBunch->getTotalNum() << endl;
            }

            IpplTimings::stopTimer(IntegrationTimer_m);
            Ippl::Comm->barrier();

        } else if(initialTotalNum_m == 2) {
            // initialTotalNum_m == 2
            // trigger SEO mode (swith off cavity) and calculate betatron osciliation tuning.
            double r_tuning[2], z_tuning[2] ;

            for(size_t i = 0; i < (itsBunch->getLocalNum()); i++) {

                // change phase space parameters from local frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
                for(int j = 0; j < 3; j++)variable_m[j] = itsBunch->R[i](j); //[x,y,z]  units: [mm]
                for(int j = 0; j < 3; j++)variable_m[j+3] = itsBunch->P[i](j); //[px,py,pz]  units: []

                if((step_m % SinglePartDumpFreq == 0)) {
                    outfTrackOrbit_m << "ID" << (itsBunch->ID[i]);
                    outfTrackOrbit_m << " " << variable_m[0] << " " << variable_m[3] << " " << variable_m[1] << " "
                         << variable_m[4] << " " << variable_m[2] << " " << variable_m[5] << endl;
                }

                double OldTheta = 0.0;

                OldTheta = calculateAngle(variable_m[0], variable_m[1]);
                r_tuning[i] = variable_m[0] * cos(OldTheta) + variable_m[1] * sin(OldTheta);
                z_tuning[i] = variable_m[2];
                turnnumber = lastTurn;

                // integrate for one step in the lab Cartesian frame (absulate value ).

                rk4(variable_m, t, dt, i);

                double NewTheta = 0.0;

                NewTheta = calculateAngle(variable_m[0], variable_m[1]);
                if((i == 0) && (NewTheta < OldTheta)) lastTurn++;

                for(int j = 0; j < 3; j++) itsBunch->R[i](j) = variable_m[j] ; //[x,y,z]  units: [mm]
                for(int j = 0; j < 3; j++) itsBunch->P[i](j) = variable_m[j+3] ; //[px,py,pz]  units: dimensionless, beta*gama

            }//end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode

            // store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
            if(step_m % SinglePartDumpFreq == 0) {
                Ttime.push_back(t * 1.0e-9);
                Tdeltz.push_back(z_tuning [1]);
                Tdeltr.push_back(r_tuning[1]  - r_tuning[0]);
                TturnNumber.push_back(turnnumber);
            }
        } else if(initialTotalNum_m == 1) {
            // initialTotalNum_m == 1 trigger single particle mode
	  IpplTimings::startTimer(IntegrationTimer_m);
	    int i = 0;
            flagNoDeletion = true;

            if((step_m % SinglePartDumpFreq == 0)) {
	      for(unsigned int i = 0; i < itsBunch->getLocalNum(); i++) {
		  outfTrackOrbit_m << "ID" <<itsBunch->ID[i];
		  outfTrackOrbit_m << " " << itsBunch->R[i](0) << " " <<itsBunch->P[i](0) << " " <<itsBunch->R[i](1) 
		       << " " << itsBunch->P[i](1) << " " <<itsBunch->R[i](2) << " " <<itsBunch->P[i](2)<< endl; 
	      }
	    }
       
            // change phase space parameters from local frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
            for(int j = 0; j < 3; j++) {
                variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
                rold_m[j] = variable_m[j]; // used for gap cross checking
                pold_m[j] = variable_m[j+3]; // used for gap cross
            }

            double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]);//[ -pi ~ pi ]

            if((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {
                turnnumber++;
                dumpEachTurn = true;
                *gmsg << "Turn " << turnnumber << endl;

                outfThetaEachTurn_m << "#Turn number = " << turnnumber << ", Time = " << t << " [ns]" << endl;
                outfThetaEachTurn_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                                    << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                                    << " " << temp_meanTheta / pi * 180
                                    << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                                    << " " << variable_m[2]
                                    << " " << variable_m[5] << endl;
            }

            //define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
            const double azimuth_angle0 = 0.0;
            const double azimuth_angle1 = 22.5 / 180.0 * pi;
            const double azimuth_angle2 = 45.0 / 180.0 * pi;

            if((oldReferenceTheta < azimuth_angle0 - deltaTheta) && (temp_meanTheta >= azimuth_angle0 - deltaTheta)) {
                outfTheta0_m << "#Turn number = " << turnnumber << ", Time = " << t << " [ns]" << endl;
                outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                outfTheta1_m << "#Turn number = " << turnnumber << ", Time = " << t << " [ns]" << endl;
                outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                outfTheta2_m << "#Turn number = " << turnnumber << ", Time = " << t << " [ns]" << endl;
                outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            oldReferenceTheta = temp_meanTheta;

            // integrate for one step in the lab Cartesian frame (absulate value ).

            flagNoDeletion = rk4(variable_m, t, dt, i);

            if(!flagNoDeletion) {
                *gmsg << "particle" << "is lost at " << step_m << "st step!" << endl;
                throw OpalException("ParallelCyclotronTracker", "the particle is out of the region of interest.");
            }
            
            for(int j = 0; j < 3; j++) itsBunch->R[i](j) = variable_m[j] ; //[x,y,z]  units: [mm]
            for(int j = 0; j < 3; j++) itsBunch->P[i](j) = variable_m[j+3] ; //[px,py,pz]  units: [] beta*gama

            //If gap crossing happens, do momenta kicking

            for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++) {
                bool tag_crossing = false;
                double DistOld = 0.0; //mm
                RFCavity * rfcav;
                if(((*sindex)->first) == "CAVITY") {
                    // here check gap cross in the list, if do , set tag_crossing to TRUE
                    for(int j = 0; j < 3; j++)
                        rnew_m[j] = variable_m[j];
                    rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
                    tag_crossing = checkGapCross(rold_m, rnew_m, rfcav, DistOld);
                }
                if(tag_crossing) {
                    double oldMomentum2  = dot(pold_m, pold_m);
                    double oldBetgam = sqrt(oldMomentum2);
                    double oldGamma = sqrt(1.0 + oldMomentum2);
                    double oldBeta = oldBetgam / oldGamma;
                    double dt1 = DistOld / (c * oldBeta * 1.0e-6); // ns
                    double dt2 = dt - dt1;

                    // retrack particle from the old postion to cavity gap point
                    // restore the old coordinates and momenta
                    for(int j = 0; j < 3; j++) {
                        variable_m[j] = rold_m[j];
                        variable_m[j+3] = pold_m[j];
                    }

                    if(dt / dt1 < 1.0e9) rk4(variable_m, t, dt1, i);

                    for(int j = 0; j < 3; j++) {
                        itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
                        itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [] beta*gama
                    }

                    //momentum kick
                    RFkick(rfcav, t, dt1, i);
		    
                    // retrack particle  from cavity gap point for the left time to finish the entire timestep
                    for(int j = 0; j < 3; j++) {
                        variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                        variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: []
                    }

                    if(dt / dt2 < 1.0e9) rk4(variable_m, t, dt2, i);

                    for(int j = 0; j < 3; j++) {
                        itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
                        itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [], beta*gama
                    }
                }// end if: gap-crossing monentum kicking at certain cavity
            }//end for: finish checking for all cavities

	    // apply the plugin elements: probe, collimator, stripper, septum 
	    applyPluginElements(turnnumber, dt);
	    
	    IpplTimings::stopTimer(IntegrationTimer_m);
        }//end if: finish one step tracking either for initialTotalNum_m==2 || initialTotalNum_m==2 || initialTotalNum_m==1 mode

        // update bunch and some parameters and output some info. after one time step.
        // reference parameters
        // pow(itsBunch->P[0](0),2.0) + pow(itsBunch->P[0](1),2.0) + pow(itsBunch->P[0](2),2.0);
        double tempP2 = dot(itsBunch->P[0], itsBunch->P[0]);
        double tempGamma = sqrt(1.0 + tempP2);
        double tempBeta = sqrt(tempP2) / tempGamma;

        PathLength_m += c_mmtns * dt / 1000.0 * tempBeta; // unit: m

        t += dt;
        itsBunch->setT((t) * 1.0e-9);
        itsBunch->setLPath(PathLength_m);
        // Here is global frame, don't do itsBunch->boundp();

        // FixMe:
        // dump phase space distribution of bunch
        // if ((((step_m+1) % Options::psDumpFreq == 0 ) && initialTotalNum_m > 2) || (dumpEachTurn && initialTotalNum_m > 2) )
        // if(step_m>(1185*(Options::rDump))) {
        //  itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);
        //}
        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2)
           || (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {
            IpplTimings::startTimer(DumpTimer_m);

            itsBunch->setTrackStep((step_m + 1));

            extE_m = Vector_t(0.0, 0.0, 0.0);
            extB_m = Vector_t(0.0, 0.0, 0.0);

            //--------------------- calculate mean coordinates  of bunch -------------------------------//
            //------------  and calculate the external field at the mass of bunch-----------------------//

            Vector_t meanR = calcMeanR();
            Vector_t meanP = calcMeanP();

            // define longitudinal direction of the bunch
            // x:[0] transverse horizontal, y:[1] longitudinal, z:[2] transverse vertical
            double phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

            double temp_cosPhi = cos(phi);
            double temp_sinPhi = sin(phi);

            double meanRadius = sqrt(meanR(0) * meanR(0) +  meanR(1) * meanR(1));

            beamline_list::iterator DumpSindex = FieldDimensions.begin();

            if(((((*DumpSindex)->second).first)[0] <= meanRadius)
               && ((((*DumpSindex)->second).first)[1] >= meanRadius)) {
                (((*DumpSindex)->second).second)->apply(meanR, Vector_t(0.0), t, extE_m, extB_m);
            }

            FDext_m[0] = extB_m / 10.0; // kgauss -> T
            FDext_m[1] = extE_m;

            //----------------------------dump in global frame-------------------------------------//
            // Note: Don't dump when
            // 1. after one turn
            // in order to sychronize the dump step for multi-bunch and single
            // bunch for compare with each other during post-process phase.
            if(!(Options::psDumpLocalFrame)) {

                itsBunch->R /= Vector_t(1000.0); // mm --> m
                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
                itsBunch->R *= Vector_t(1000.0); // m --> mm
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (global frame) at integration step "
                      << step_m + 1 << " T= " << t << " [nS]" << endl;
                //----------------------------dump in local frame-------------------------------------//
            } else {
                itsBunch->R -= meanR;
                itsBunch->P -= meanP;

                for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {

                    double temp_RHorizontal     =   itsBunch->R[ii](0) * temp_cosPhi  + itsBunch->R[ii](1) * temp_sinPhi;
                    double temp_RLongitudinal   =  -itsBunch->R[ii](0) * temp_sinPhi  + itsBunch->R[ii](1) * temp_cosPhi;

                    itsBunch->R[ii](0) = temp_RHorizontal;
                    itsBunch->R[ii](1) = temp_RLongitudinal;

                    double temp_PHorizontal     =   itsBunch->P[ii](0) * temp_cosPhi  + itsBunch->P[ii](1) * temp_sinPhi;
                    double temp_PLongitudinal   =  -itsBunch->P[ii](0) * temp_sinPhi  + itsBunch->P[ii](1) * temp_cosPhi;

                    itsBunch->P[ii](0) = temp_PHorizontal;
                    itsBunch->P[ii](1) = temp_PLongitudinal;
                }

                itsBunch->R /= Vector_t(1000.0); // mm --> m
                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
                itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0);
                itsBunch->R *= Vector_t(1000.0); // m --> mm


                for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {

                    double temp_Rx   =  -itsBunch->R[ii](1) * temp_sinPhi  + itsBunch->R[ii](0) * temp_cosPhi;
                    double temp_Ry   =   itsBunch->R[ii](1) * temp_cosPhi  + itsBunch->R[ii](0) * temp_sinPhi;

                    itsBunch->R[ii](0) = temp_Rx;
                    itsBunch->R[ii](1) = temp_Ry;

                    double temp_Px   =  -itsBunch->P[ii](1) * temp_sinPhi  + itsBunch->P[ii](0) * temp_cosPhi;
                    double temp_Py   =   itsBunch->P[ii](1) * temp_cosPhi  + itsBunch->P[ii](0) * temp_sinPhi;

                    itsBunch->P[ii](0) = temp_Px;
                    itsBunch->P[ii](1) = temp_Py;

                }

                itsBunch->R += meanR;
                itsBunch->P += meanP;
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (local frame) at integration step "
                      << step_m + 1 << " T= " << t << " [nS]" << endl;
            }
            IpplTimings::stopTimer(DumpTimer_m);
        }
        if(!(step_m + 1 % 1000))
            *gmsg << "Step " << step_m + 1 << endl;

    }// end for: the integration is DONE after maxSteps_m steps!

    // some post-integration works
    *gmsg << "* *---------------------------- PARTICLES TRACK DONE------ ----------------------------*** " << endl;

    // calculate tunes after tracking.

    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 0) {
            double FinalMomentum2 = pow(itsBunch->P[ii](0), 2.0) + pow(itsBunch->P[ii](1), 2.0) + pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsgAll << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsgAll << "* Total phase space dump number(includes the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsgAll << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
        }
    }

    Ippl::Comm->barrier();

    if(initialTotalNum_m == 2) {
        *gmsg << "* ******** The result for tune calulation (NO space charge) ********** " << endl
              << "* Number of tracked turns: " << TturnNumber.back() << endl;
        double nur, nuz;
        getTunes(Ttime, Tdeltr, Tdeltz, TturnNumber.back(), nur, nuz);
	
    } else {
        // not for multibunch
        if(!(itsBunch->weHaveBins()))
            *gmsg << "* Total finished turn number (not correct for restart mode) = " << turnnumber << endl;
    }

    Ippl::Comm->barrier();

    if(myNode_m == 0) outfTrackOrbit_m.close();

    if(initialTotalNum_m == 1)
        closeFiles();

    *gmsg << *itsBunch << endl;

    // free memory
    if(gmsgAll)
        free(gmsgAll);
}


/**
 *
 *
 * @param y
 * @param t
 * @param yp
 * @param Pindex
 *
 * @return
 */
bool ParallelCyclotronTracker::derivate(double *y, double t, double *yp, int Pindex) {
    Vector_t externalE, externalB, tempR;

    double partR;

    // a flag of visit flag, if cyclotron is not visited, set false.
    // bool visitflag = false;


    externalB = Vector_t(0.0, 0.0, 0.0);
    externalE = Vector_t(0.0, 0.0, 0.0);

    for(int i = 0; i < 3; i++)
        tempR(i) = y[i];

    partR = sqrt(tempR(0) * tempR(0) + tempR(1) * tempR(1));

    beamline_list::iterator sindex = FieldDimensions.begin();

    //debug
    if(((*sindex)->first) != "CYCLOTRON") {
        ERRORMSG("Error in ParallelCyclotronTracker, the CYCLOTRON object is not the first element in the beamline_list! " << endl);
        exit(1);
    }

    if(((((*sindex)->second).first)[0] <= partR) && ((((*sindex)->second).first)[1] >= partR))
        (((*sindex)->second).second)->apply(tempR, Vector_t(0.0), t, externalE, externalB);
    else
        return false;

    for(int i = 0; i < 3; i++) externalB(i) = externalB(i) * 0.10; //[kGauss] -> [T]
    for(int i = 0; i < 3; i++) externalE(i) = externalE(i) * 1.0e6; //[kV/mm ] -> [V/m]

    // for working modes without space charge effects, override this step to save time
    if(itsBunch->hasFieldSolver()) {
        if(itsBunch->ID[Pindex] != 0) {
            // add external Field and self space charge field
            externalE += itsBunch->Ef[Pindex];
            externalB += itsBunch->Bf[Pindex];
        }

    }

    double qtom = itsBunch->Q[Pindex] / (itsBunch->M[Pindex] * mass_coeff);   // m^2/s^2/GV

    double tempgamma = sqrt(1 + (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]));

    /*
       d2x/dt2 = q/m * ( E + v x B )
       dx/dt =vx
       units: mm, and []
    */

    yp[0] = c_mmtns / tempgamma * y[3]; // [mn/ns]
    yp[1] = c_mmtns / tempgamma * y[4]; // [mm/ns]
    yp[2] = c_mmtns / tempgamma * y[5]; // [mm/ns]

    yp[3] = (externalE(0) / c  + (externalB(2) * y[4] - externalB(1) * y[5]) / tempgamma) * qtom; // [1/ns]
    yp[4] = (externalE(1) / c  - (externalB(2) * y[3] - externalB(0) * y[5]) / tempgamma) * qtom; // [1/ns];
    yp[5] = (externalE(2) / c  + (externalB(1) * y[3] - externalB(0) * y[4]) / tempgamma) * qtom; // [1/ns];

    return true;
}


/**
 *
 *
 * @param x
 * @param t
 * @param tau
 * @param Pindex
 *
 * @return
 */
bool ParallelCyclotronTracker::rk4(double x[], double t, double tau, int Pindex) {
    // Forth order Runge-Kutta integrator
    // arguments:
    //   x          Current value of dependent variable
    //   t          Independent variable (usually time)
    //   tau        Step size (usually time step)
    //   Pindex     index of particel, not used yet

    bool visitflag ;
    double  deriv1[PSdim];
    double  deriv2[PSdim];
    double  deriv3[PSdim];
    double  deriv4[PSdim];
    double  xtemp[PSdim];

    visitflag = true;

    // Evaluate f1 = f(x,t).

    visitflag = derivate(x, t, deriv1 , Pindex);

    if(!visitflag) return false;

    // Evaluate f2 = f( x+tau*f1/2, t+tau/2 ).
    double half_tau = 0.5 * tau;
    double t_half = t + half_tau;

    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + half_tau * deriv1[i];

    visitflag = derivate(xtemp, t_half, deriv2 , Pindex);

    if(!visitflag) return false;

    // Evaluate f3 = f( x+tau*f2/2, t+tau/2 ).
    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + half_tau * deriv2[i];

    visitflag = derivate(xtemp, t_half, deriv3 , Pindex);

    if(!visitflag) return false;

    // Evaluate f4 = f( x+tau*f3, t+tau ).
    double t_full = t + tau;
    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + tau * deriv3[i];

    visitflag = derivate(xtemp, t_full, deriv4 , Pindex);

    if(!visitflag) return false;

    // Return x(t+tau) computed from fourth-order R-K.
    for(int i = 0; i < PSdim; i++)
        x[i] += tau / 6.*(deriv1[i] + deriv4[i] + 2.*(deriv2[i] + deriv3[i]));

    return true;

}

/**
 *
 *
 * @param Rold
 * @param Rnew
 * @param elptr
 * @param Dold
 *
 * @return
 */

bool ParallelCyclotronTracker::checkGapCross(Vector_t Rold, Vector_t Rnew, RFCavity * rfcavity, double &Dold) {
    bool flag = false;
    double sinx = rfcavity->getSinAzimuth();
    double cosx = rfcavity->getCosAzimuth();
    double PerpenDistance = rfcavity->getPerpenDistance();
    double distNew = (Rnew[0] * sinx - Rnew[1] * cosx) - PerpenDistance;
    double distOld = (Rold[0] * sinx - Rold[1] * cosx) - PerpenDistance;
    if(distOld > 0.0 && distNew <= 0.0) flag = true;
    // This parameter is used correct cavity phase
    Dold = distOld;
    return flag;
}

bool ParallelCyclotronTracker::RFkick(RFCavity * rfcavity, const double t, const double dt, const int Pindex){
  double radius = sqrt(pow(itsBunch->R[Pindex](0), 2.0) + pow(itsBunch->R[Pindex](1), 2.0)
                  - pow(rfcavity->getPerpenDistance() , 2.0));
  double rmin = rfcavity->getRmin();
  double rmax = rfcavity->getRmax();
  double nomalRadius = (radius - rmin) / (rmax - rmin);
  double tempP[3];
  if(nomalRadius <= 1.0 && nomalRadius >= 0.0) {

    for(int j = 0; j < 3; j++)
      tempP[j] = itsBunch->P[Pindex](j);  //[px,py,pz]  units: dimensionless
    
    // here evaluate voltage and conduct momenta kicking;
    rfcavity -> getMomentaKick2(nomalRadius, tempP, t, dt, itsBunch->ID[Pindex], itsBunch->getM(), itsBunch->getQ()); // t : ns
    
    for(int j = 0; j < 3; j++)
      itsBunch->P[Pindex](j) = tempP[j];
    return true;
  }
  return false;
}


struct adder : public unary_function<double, void> {
    adder() : sum(0) {}
    double sum;
    void operator()(double x) { sum += x; }
};

/**
 *
 *
 * @param t
 * @param r
 * @param z
 * @param lastTurn
 * @param nur
 * @param nuz
 *
 * @return
 */
bool ParallelCyclotronTracker::getTunes(vector<double> &t, vector<double> &r, vector<double> &z,
                                        int lastTurn, double &nur, double &nuz) {
    TUNE_class *tune;

    int Ndat = t.size();

    /*
      remove mean
    */
    double rsum =  for_each(r.begin(), r.end(), adder()).sum;

    for(int i = 0; i < Ndat; i++)
        r[i] -= rsum;

    double zsum =  for_each(z.begin(), z.end(), adder()).sum;

    for(int i = 0; i < Ndat; i++)
        z[i] -= zsum;
    double ti = *(t.begin());
    double tf = t[t.size()-1];
    double T = (tf - ti);

    t.clear();
    double dt = T / Ndat;
    double time = 0.0;
    for(int i = 0; i < Ndat; i++) {
        t.push_back(i);
        time += dt;
    }

    T = t[Ndat-1];

    *gmsg << "* *************** nuR ***************" << endl;
    *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

    int nhis_lomb = 10;
    int  stat = 0;
    // book tune class
    tune = new TUNE_class();
    stat = tune->LombAnalysis(t, r, nhis_lomb, T / lastTurn);
    if(stat != 0)
        *gmsg << "* TUNE: Lomb analysis failed" << endl;
    *gmsg << "* ***********************************" << endl << endl;

    delete tune;
    tune = NULL;
    // FixMe: need to come from the inputfile
    nhis_lomb = 100;

    if(zsum != 0.0) {
        *gmsg << "* *************** nuZ ***************" << endl;
        *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

        // book tune class
        tune = new TUNE_class();
        stat = tune->LombAnalysis(t, z, nhis_lomb, T / lastTurn);
        if(stat != 0)
            *gmsg << "* TUNE: Lomb analysis failed" << endl;
        *gmsg << " ***********************************" << endl << endl;

        delete tune;
        tune = NULL;
    }
    return true;
}

void ParallelCyclotronTracker::saveOneBunch() {
    *gmsg << "---------------- Clone Beam----------------" << endl;

    npart_mb = itsBunch->getLocalNum();

    r_mb.create(npart_mb);
    r_mb = itsBunch->R;

    p_mb.create(npart_mb);
    p_mb = itsBunch->P;

    m_mb.create(npart_mb);
    m_mb = itsBunch->M;

    q_mb.create(npart_mb);
    q_mb = itsBunch->Q;

    ptype_mb.create(npart_mb);
    ptype_mb = itsBunch->PType;

    string fn_appendix = "-onebunch";
    itsDataSink->storeOneBunch(*itsBunch, fn_appendix);

}

/**
 *
 * @param BinID
 * @param step
 */
bool ParallelCyclotronTracker::readOneBunch(const size_t BinID) {
    *gmsg << "---------------- Copy Beam----------------" << endl;

    const size_t LocalNum = itsBunch->getLocalNum();
    const size_t NewLocalNum = LocalNum + npart_mb;

    itsBunch->create(npart_mb);

    for(size_t ii = LocalNum; ii < NewLocalNum; ii++) {
        size_t i = ii - LocalNum;
        itsBunch->R[ii] = r_mb[i];
        itsBunch->P[ii] = p_mb[i];
        itsBunch->M[ii] = m_mb[i];
        itsBunch->Q[ii] = q_mb[i];
        itsBunch->PType[ii] = ptype_mb[i];
        itsBunch->Bin[ii] = BinID;

    }

    // update statistics parameters of PartBunch
    // allocate ID for new particles
    itsBunch->boundp();

    return true;

}

bool ParallelCyclotronTracker::readOneBunchFromFile(const size_t BinID) {
    //todo
    static bool restartflag = true;

    if(restartflag) {
        *gmsg << "----------------Read beam from hdf5 file----------------" << endl;
        const size_t LocalNum = itsBunch->getLocalNum();

        string fn_appendix = "-onebunch";
        itsDataSink->readOneBunch(*itsBunch, fn_appendix, BinID);
        restartflag = false;

        const size_t NewLocalNum = itsBunch->getLocalNum();

        npart_mb = NewLocalNum - LocalNum;

        r_mb.create(npart_mb);
        p_mb.create(npart_mb);
        m_mb.create(npart_mb);
        q_mb.create(npart_mb);
        ptype_mb.create(npart_mb);

        for(size_t ii = 0; ii < npart_mb; ii++) {
            r_mb[ii] = itsBunch->R[ii+LocalNum];
            p_mb[ii] = itsBunch->P[ii+LocalNum];
            m_mb[ii] = itsBunch->M[ii+LocalNum];
            q_mb[ii] = itsBunch->Q[ii+LocalNum];
            ptype_mb[ii] = itsBunch->PType[ii+LocalNum];
        }

    } else
        readOneBunch(BinID);

    itsBunch->boundp();

    return true;
}

void ParallelCyclotronTracker::Tracker_MTS() {
    Inform *gmsgAll;
    gmsgAll = new Inform("CycTracker MTS", INFORM_ALL_NODES);
    beamline_list::iterator sindex = FieldDimensions.begin();
    Cyclotron *elptr = dynamic_cast<Cyclotron *>(((*sindex)->second).second);
    if (elptr == NULL)         
        throw OpalException("ParallelCyclotronTracker::Tracker_LF()",
                            "The first item in the FieldDimensions list does not seem to be a cyclotron element");

    const double harm = elptr->getCyclHarm();
    const double dt = itsBunch->getdT() * harm;
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << itsBunch->getStepsPerTurn() * dt * 1.0e9 << "[ns]" << endl;
    }

    initTrackOrbitFile();

    // parameter for reset bin in multi-bunch run, todo: readin from inputfile
    const double eta = 0.01;

    // get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun() && numBunch_m > 1) {
        restartStep0_m = itsBunch->getTrackStep();
        step_m = restartStep0_m;
        itsBunch->resetPartBinID2(eta);
        *gmsg << "Restart at integration step " << restartStep0_m << endl;
    }
    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        itsBunch->setT(0.0);
    }
    *gmsg << "* Beginning of this run is at t= " << itsBunch->getT() * 1e9 << " [ns]" << endl;
    *gmsg << "* The time step is set to dt= " << dt * 1e9 << " [ns]" << endl;

    // for single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1) {
        openFiles(inputFileNameWithoutExtension());
    }

	// In this method, we use [R] = m
	referenceR *= 0.001;

    double const initialReferenceTheta = referenceTheta / 180.0 * pi;
    if(!OpalData::getInstance()->inRestartRun()) {
        PathLength_m = 0.0;
        // Force the initial phase space values of the particle with ID=0 to zero, to set it as a reference particle.
        if(initialTotalNum_m > 2) {  // only for mulit-particle mode
            for(unsigned int ii = 0; ii < initialLocalNum_m; ii++) {
                if(itsBunch->ID[ii] == 0) {
                    itsBunch->R[ii] = Vector_t(0.0);
                    itsBunch->P[ii] = Vector_t(0.0);
                }
            }
        }

        if(Options::psDumpLocalFrame) {
            // dump the initial distribution
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
            itsDataSink->writeStatData(*itsBunch, FDext_m, 0, 0, 0);
            *gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the local frame " << endl;
        }

        // Initialize global R
        Vector_t const initMeanR = Vector_t(referenceR * cosRefTheta_m, referenceR * sinRefTheta_m, 0.0);
        localToGlobal(itsBunch->R, initialReferenceTheta, initMeanR);

        // Initialize global P
        for(unsigned int ii = 0; ii < initialLocalNum_m; ++ii) {
            itsBunch->P[ii](0) += referencePr;
            itsBunch->P[ii](1) += referencePt;
        }
        localToGlobal(itsBunch->P, initialReferenceTheta);

        // Initialize the bin number of the first bunch to 0
        for(unsigned int ii = 0; ii < initialLocalNum_m; ++ii) {
            itsBunch->Bin[ii] = 0;
        }

        if(!(Options::psDumpLocalFrame)) {

            // dump the initial distribution
            lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);

            *gmsg << "meanR=( " << initMeanR(0) * 1000.0 << " " << initMeanR(1) * 1000.0 << " " << initMeanR(2) * 1000.0 << " ) [mm] " << endl;
            *gmsg << "phase space dumping in glabol frame at initial position " << endl;
            *gmsg << "Dump step = " << lastDumpedStep_m << endl;
        }

        // AUTO mode
        if(multiBunchMode_m == 2) {
            Vector_t Rmean = itsBunch->get_centroid();
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;
            *gmsg << "initial Radial position = " << RThisTurn_m * 1000.0 << " [mm]" << endl;
        }

        if((initialTotalNum_m > 2) && (numBunch_m > 1) && (multiBunchMode_m == 1)) {

            // backup initial distribution
            initialR_m = new Vector_t[initialLocalNum_m];
            initialP_m = new Vector_t[initialLocalNum_m];
            for(unsigned int ii = 0; ii < initialLocalNum_m; ++ii) {
                initialR_m[ii] = itsBunch->R[ii];
                initialP_m[ii] = itsBunch->P[ii];
            }
        }

    } else {
        PathLength_m = itsBunch->getLPath();

        // AUTO mode
        if(multiBunchMode_m == 2) {
            itsBunch->calcBeamParameters_cycl();
			Vector_t Rmean = itsBunch->get_centroid();
            RLastTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
            RThisTurn_m = RLastTurn_m;
            *gmsg << "Radial position at restart position = " << RThisTurn_m * 1000.0 << " [mm]" << endl;
        }
    }
    *gmsg << "single particle trajectory dump frequency is set to " << Options::sptDumpFreq << endl;
    *gmsg << "particles repartition frequency is set to " << Options::repartFreq << endl;
    if(numBunch_m > 1) {
        *gmsg << "particles energy bin ID reset frequency is set to " << Options::rebinFreq << endl;
    }
    if(initialTotalNum_m == 1) {
        *gmsg << "* *---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," << endl
              << "* The initial distribution file must be specified which should contain only one line for the single particle " << endl
              << "* *------------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
        *gmsg << "* *------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." << endl
              << "* This mode does NOT include any RF cavities. The initial distribution file must be specified" << endl
              << "* In the file the first line is for reference particle and the second line is for offcenter particle." << endl
              << "* The tunes are calculated by FFT routines based on these two particles. " << endl
              << "* *------------NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    } else {
        Vector_t const meanP = calcMeanP();
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
        Vector_t const meanR = calcMeanR();
        globalToLocal(itsBunch->R, phi, meanR);
        itsBunch->boundp();
        double const meanGamma = sqrt(1.0 + pow(meanP(0), 2.0) + pow(meanP(1), 2.0));
        itsBunch->Bf = Vector_t(0.0);
        itsBunch->Ef = Vector_t(0.0);
        itsBunch->computeSelfFields_cycl(meanGamma);
        localToGlobal(itsBunch->R, phi, meanR);
    }

    // do repartation for restart run
    if(OpalData::getInstance()->inRestartRun() && numBunch_m > 1) {
        Vector_t const meanP = calcMeanP();
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
        Vector_t const meanR = calcMeanR();
        globalToLocal(itsBunch->R, phi, meanR);
        itsBunch->boundp();
        *gmsg << "Do binaryRepartation for restart run" << endl;
        checkNumPart(string("Before repartation: "), itsBunch->getLocalNum());
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        IpplTimings::stopTimer(BinRepartTimer_m);
        Ippl::Comm->barrier();
        checkNumPart(string("After repartation: "), itsBunch->getLocalNum());
        *gmsg << "Done binaryRepartation for restart run" << endl;
        localToGlobal(itsBunch->R, phi, meanR);
    }

    int const numSubsteps = std::max(Options::mtsSubsteps, 1);
    *gmsg << "MTS: Number of substeps per step is " << numSubsteps << endl;
    double const dt_inner = dt / double(numSubsteps);
    *gmsg << "MTS: The inner time step is therefore " << dt_inner << endl;
    int SteptoLastInj = itsBunch->getSteptoLastInj();
    double oldReferenceTheta = initialReferenceTheta;
    int turnnumber = 1;
    bool flagTransition = false; // flag to determine when to transit from single-bunch to multi-bunches mode
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn(); // step point determining the next time point of check for transition
    const double deltaTheta = pi / itsBunch->getStepsPerTurn();
    *gmsg << "---------------------------- Start tracking ----------------------------" << endl;
    if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
        evaluateSpaceChargeField();
    }
    for(; step_m < maxSteps_m; step_m++) {
        bool dumpEachTurn = false;
        if(step_m % Options::sptDumpFreq == 0) {
            itsBunch->R *= Vector_t(1000.0);
            singleParticleDump();
            itsBunch->R *= Vector_t(0.001);
        }
        Ippl::Comm->barrier();

        // First half kick from space charge force
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            kick(0.5 * dt);
        }

        // Substeps for external field integration
        for(int n = 0; n < numSubsteps; ++n) {
            borisExternalFields(dt_inner);
        }

        // bunch injection
        // TODO: Where is correct location for this piece of code? Beginning/end of step? Before field solve?
        if(numBunch_m > 1) {
            if((BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)) {
                if(step_m == stepsNextCheck) {
                    // under 3 conditions, following code will be execute
                    // to check the distance between two neighborring bunches
                    // 1.multi-bunch mode, AUTO sub-mode
                    // 2.After each revolution
                    // 3.only one bunch exists
                    *gmsg << "checking for automatically injecting new bunch ..." << endl;
                    itsBunch->calcBeamParameters_cycl();
                    Vector_t Rmean = itsBunch->get_centroid();
                    RThisTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));
                    Vector_t Rrms = itsBunch->get_rrms();
                    double XYrms = sqrt(pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0));

                    // if the distance between two neighbour bunch is less than CoeffDBunches_m times of its 2D rms size
                    // start multi-bunch simulation, fill current phase space to initialR and initialP arrays
                    if((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                        // since next turn, start multi-bunches
                        saveOneBunch();
                        flagTransition = true;
                        *gmsg << "*** Save beam distribution at turn #" << turnnumber << " ***" << endl;
                        *gmsg << "*** After one revolution, Multi-Bunch Mode will be invorked ***" << endl;
                    }

                    stepsNextCheck += itsBunch->getStepsPerTurn();
                    *gmsg << "RLastTurn = " << RLastTurn_m * 1000.0 << " [mm]" << endl;
                    *gmsg << "RThisTurn = " << RThisTurn_m * 1000.0 << " [mm]" << endl;
                    *gmsg << "    XYrms = " << XYrms * 1000.0 << " [mm]" << endl;
                    RLastTurn_m = RThisTurn_m;
                }
            } else if(SteptoLastInj == itsBunch->getStepsPerTurn() - 1) {
                if(BunchCount_m < numBunch_m) {
                    // under 4 conditions, following code will be execute
                    // to read new bunch from hdf5 format file for FORCE or AUTO mode
                    // 1.multi-bunch mode
                    // 2.after each revolution
                    // 3.existing bunches is less than the specified bunches
                    // 4.FORCE mode, or AUTO mode with flagTransition = true
                    // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
                    *gmsg << "step " << step_m << ", inject a new bunch... ... ..." << endl;
                    BunchCount_m++;

                    // read initial distribution from h5 file
                    if(multiBunchMode_m == 1) {
                        readOneBunch(BunchCount_m - 1);
                        itsBunch->resetPartBinID2(eta);
                    } else if(multiBunchMode_m == 2) {
                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        itsBunch->resetPartBinID2(eta);
                    }
                    SteptoLastInj = 0;
                    itsBunch->setNumBunch(BunchCount_m);
                    stepsNextCheck += itsBunch->getStepsPerTurn();

                    // update  after injection
                    itsBunch->boundp();

                    Ippl::Comm->barrier();
                    *gmsg << BunchCount_m << "'th bunch injected, total particle number = " << itsBunch->getTotalNum() << endl;
                }
            } else if(BunchCount_m == numBunch_m) {
                // After this, numBunch_m is wrong but useless
                numBunch_m--;
            } else {
                SteptoLastInj++;
            }
        }

        // calculate self fields Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            evaluateSpaceChargeField();
        } else {
            // if field solver is not available , only update bunch, to transfer particles between nodes if needed,
            // reset parameters such as LocalNum, initialTotalNum_m.
            // INFOMSG("No space charge Effects are included!"<<endl;);
            if((step_m % Options::repartFreq * 100) == 0 && initialTotalNum_m >= 1000) { //TODO: why * 100?
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                Vector_t const meanR = calcMeanR();
                globalToLocal(itsBunch->R, phi, meanR);
                itsBunch->boundp();
                repartition();
                localToGlobal(itsBunch->R, phi, meanR);
            }
        }

        // Second half kick from space charge force
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            kick(0.5 * dt);
        }

        // apply the plugin elements: probe, collimator, stripper, septum
        itsBunch->R *= Vector_t(1000.0); // applyPluginElements expects [R] = mm
        bool const flagNeedUpdate = applyPluginElements(turnnumber, dt * 1e9);
        itsBunch->R *= Vector_t(0.001);
        if(itsBunch->weHaveBins() && flagNeedUpdate)
            itsBunch->resetPartBinID2(eta);

        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
            if(step_m % Options::rebinFreq == 0) {
                itsBunch->resetPartBinID2(eta);
            }
        }

        // dump some data after one push in single particle tracking
        if(initialTotalNum_m == 1) {
            int i = 0;

            // change phase space parameters from local reference frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
            for(int j = 0; j < 3; j++) {
                variable_m[j]   = itsBunch->R[i](j) * 1000.0;  //[x,y,z]  units: [mm]
                variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
            }

            double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]);//[ -pi ~ pi ]
            if((oldReferenceTheta < initialReferenceTheta - deltaTheta) &&
               (temp_meanTheta >= initialReferenceTheta - deltaTheta)) {
                turnnumber++;
                *gmsg << "Turn " << turnnumber << endl;
                dumpEachTurn = true;
                outfThetaEachTurn_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfThetaEachTurn_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                                    << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                                    << " " << temp_meanTheta / pi * 180.0
                                    << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                                    << " " << variable_m[2]
                                    << " " << variable_m[5] << endl;
            }
            // FixMe: should be defined elesewhere !
            // define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
            const double azimuth_angle0 = 0.0;
            const double azimuth_angle1 = 22.5 / 180.0 * pi;
            const double azimuth_angle2 = 45.0 / 180.0 * pi;
            if((oldReferenceTheta < azimuth_angle0 - deltaTheta) && (temp_meanTheta >= azimuth_angle0 - deltaTheta)) {
                outfTheta0_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                outfTheta1_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                outfTheta2_m << "#Turn number = " << turnnumber << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }
            oldReferenceTheta = temp_meanTheta;
        }

        // check whether one turn over for multi-bunch tracking.
        if(Options::psDumpEachTurn && initialTotalNum_m > 2) {
            Vector_t const meanR = calcMeanR();

            // in global Cartesian frame, calculate the location in global frame of bunch
            oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

            // both for single bunch and multi-bunch
            // avoid dump at the first step
            // dumpEachTurn has not been changed in first push
            if((step_m > 10) && ((step_m + 1) % itsBunch->getStepsPerTurn()) == 0) {
                ++turnnumber;
                dumpEachTurn = true;
                *gmsg << "Turn " << turnnumber << " total particles " << itsBunch->getTotalNum() << endl;
            }
        }
        // reset Bin ID for each particle
        if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % Options::rebinFreq == 0)
            itsBunch->resetPartBinID2(eta);

        // dump phase space distribution of bunch
        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2) ||
           (Options::psDumpEachTurn && dumpEachTurn && initialTotalNum_m != 2))
        {
            IpplTimings::startTimer(DumpTimer_m);
            itsBunch->setSteptoLastInj(SteptoLastInj);
            itsBunch->setTrackStep((step_m + 1));

            //--------------------- calculate mean coordinates of bunch -------------------------------//
            //------------ and calculate the external field at the mass of bunch -----------------------//

            Vector_t const meanR = calcMeanR();
            Vector_t const meanP = calcMeanP();

            // define longitudinal direction of the bunch
            double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

            *gmsg << "meanR=( " << meanR(0) * 1000.0 << " " << meanR(1) * 1000.0 << " " << meanR(2) * 1000.0 << " ) [mm] " << endl;
            double meanRadius = sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1)) * 1000.0;
			extE_m = Vector_t(0.0, 0.0, 0.0);
            extB_m = Vector_t(0.0, 0.0, 0.0);
            beamline_list::iterator DumpSindex = FieldDimensions.begin();
            if(((((*DumpSindex)->second).first)[0] <= meanRadius) &&
               ((((*DumpSindex)->second).first)[1] >= meanRadius)) {
                (((*DumpSindex)->second).second)->apply(meanR * 1000.0, Vector_t(0.0), itsBunch->getT() * 1e9, extE_m, extB_m);
            }
            FDext_m[0] = extB_m * 0.1; // kgauss -> T
            FDext_m[1] = extE_m;

            // Note: Don't dump when
            // 1. after one turn
            // in order to sychronize the dump step for multi-bunch and single bunch for compare
            // with each other during post-process phase.
            if(!(Options::psDumpLocalFrame)) {

                // dump in global frame
                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
                //  itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);
                // TODO: why no stat data in global frame?
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (global frame) at integration step "
                      << step_m + 1 << " T= " << itsBunch->getT() * 1e9 << " [ns]" << endl;

            } else {
                globalToLocal(itsBunch->R, phi, meanR);
                globalToLocal(itsBunch->P, phi, meanP);

                // dump in local frame
                lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
                itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0);
                *gmsg << "* Phase space dump " << lastDumpedStep_m << " (local frame) at integration step "
                      << step_m + 1 << " T= " << itsBunch->getT() * 1e9 << " [ns]" << endl;

                localToGlobal(itsBunch->R, phi, meanR);
                localToGlobal(itsBunch->P, phi, meanP);
            }
            IpplTimings::stopTimer(DumpTimer_m);
        }
    }

    for(size_t ii = 0; ii < itsBunch->getLocalNum(); ++ii) {
        if(itsBunch->ID[ii] == 0) {
            // FixMe: FinalMomentum2  = dot(itsBunch->P[ii],itsBunch->P[ii]);
            double FinalMomentum2  = pow(itsBunch->P[ii](0), 2.0) +
                                     pow(itsBunch->P[ii](1), 2.0) +
                                     pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsgAll << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsgAll << "* Total phase space dump number including the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsgAll << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
        }
    }
    Ippl::Comm->barrier();
    if(myNode_m == 0) outfTrackOrbit_m.close();
    if(initialTotalNum_m == 1) closeFiles();
    *gmsg << *itsBunch << endl;

    // free memory
    if(gmsgAll) free(gmsgAll);
}

Vector_t ParallelCyclotronTracker::calcMeanR() const {
    Vector_t meanR(0.0, 0.0, 0.0);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanR(d) += itsBunch->R[i](d);
        }
    }
    reduce(meanR, meanR, OpAddAssign());
    return meanR / Vector_t(itsBunch->getTotalNum());
}

Vector_t ParallelCyclotronTracker::calcMeanP() const {
    Vector_t meanP(0.0, 0.0, 0.0);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanP(d) += itsBunch->P[i](d);
        }
    }
    reduce(meanP, meanP, OpAddAssign());
    return meanP / Vector_t(itsBunch->getTotalNum());
}

void ParallelCyclotronTracker::repartition() {
    if((step_m % Options::repartFreq) == 0) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
    }
}

void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors, double phi, Vector_t const translationToGlobal) {
    if(translationToGlobal != 0) {
        particleVectors -= translationToGlobal;
    }
    Tenzor<double, 3> const rotation( cos(phi), sin(phi), 0, 
                                     -sin(phi), cos(phi), 0, 
                                             0,        0, 1); // clockwise rotation
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors, double phi, Vector_t const translationToGlobal) {
    Tenzor<double, 3> const rotation(cos(phi), -sin(phi), 0, 
                                     sin(phi),  cos(phi), 0, 
                                            0,         0, 1); // counter-clockwise rotation
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
    if(translationToGlobal != 0) {
        particleVectors += translationToGlobal;
    }
}

void ParallelCyclotronTracker::push(double h) {
    IpplTimings::startTimer(IntegrationTimer_m);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        Vector_t const oldR = itsBunch->R[i];
        double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
        double const c_gamma = c / gamma;
        Vector_t const v = itsBunch->P[i] * c_gamma;
        itsBunch->R[i] += h * v;
        for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++) {
            if(((*sindex)->first) == "CAVITY") {
                double distOld;
                RFCavity * rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
                bool const tagCrossing = checkGapCross(oldR * 1000.0, itsBunch->R[i] * 1000.0, rfcav, distOld);
                distOld *= 0.001; // checkGapCross gives [distOld] == mm
                if(tagCrossing) {
                    double const dt1 = distOld / sqrt(dot(v, v));
                    double const dt2 = h - dt1;
                    
                    // Retrack particle from the old postion to cavity gap point
                    itsBunch->R[i] = oldR + dt1 * v;
                    
                    // Momentum kick
                    itsBunch->R[i] *= 1000.0; // RFkick uses [itsBunch->R] == mm
                    RFkick(rfcav, itsBunch->getT() * 1.0e9, dt1 * 1.0e9, i);
                    itsBunch->R[i] *= 0.001;
                    
                    itsBunch->R[i] += dt2 * itsBunch->P[i] * c_gamma;
                }
            }
        }
    }
    itsBunch->setT(itsBunch->getT() + h);
    
    // Path length update
    double const dotP = dot(itsBunch->P[0], itsBunch->P[0]);
    double const gamma = sqrt(1.0 + dotP);
    PathLength_m += h * sqrt(dotP) * c / gamma;
    itsBunch->setLPath(PathLength_m);
    
    IpplTimings::stopTimer(IntegrationTimer_m);
}

void ParallelCyclotronTracker::kick(double h) {
    IpplTimings::startTimer(IntegrationTimer_m);
    double const q = itsBunch->Q[0] / q_e; // For now all particles have the same charge
    double const M = itsBunch->M[0] * 1.0e9; // For now all particles have the same rest energy
    double const h12Halfqc_M = 0.5 * h * q * c / M;
    double const h12Halfqcc_M = h12Halfqc_M * c;
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        // Half step E
        itsBunch->P[i] += h12Halfqc_M * itsBunch->Ef[i];

        // Full step B
        double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
        Vector_t const r = h12Halfqcc_M * itsBunch->Bf[i] / gamma;
        Vector_t const w = itsBunch->P[i] + cross(itsBunch->P[i], r);
        Vector_t const s = 2.0 / (1.0 + dot(r, r)) * r;
        itsBunch->P[i] += cross(w, s);

        // Half step E
        itsBunch->P[i] += h12Halfqc_M * itsBunch->Ef[i];

    }
    IpplTimings::stopTimer(IntegrationTimer_m);
}

void ParallelCyclotronTracker::borisExternalFields(double h) {

    // push particles for first half step
    push(0.5 * h);

    // Evaluate external fields
    IpplTimings::startTimer(IntegrationTimer_m);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        double const partR = sqrt(dot(itsBunch->R[i], itsBunch->R[i])) * 1000.0;
        itsBunch->Ef[i] = Vector_t(0.0, 0.0, 0.0);
        itsBunch->Bf[i] = Vector_t(0.0, 0.0, 0.0);
        beamline_list::iterator sindex = FieldDimensions.begin();
        if(((((*sindex)->second).first)[0] <= partR) && ((((*sindex)->second).first)[1] >= partR)) {
            (((*sindex)->second).second)->apply(itsBunch->R[i] * 1000.0, Vector_t(0.0), itsBunch->getT() * 1e9, itsBunch->Ef[i], itsBunch->Bf[i]);
        }
        itsBunch->Bf[i] *= 0.1; // kgauss -> T
    }
    IpplTimings::stopTimer(IntegrationTimer_m);

    // Kick particles for full step
    kick(h);

    // push particles for second half step
    push(0.5 * h);
}

bool ParallelCyclotronTracker::applyPluginElements(const int turnnumber, const double dt) {

  bool flagNeedUpdate = false;
  for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++) {
    if(((*sindex)->first) == "SEPTUM")    {
      (static_cast<Septum *>(((*sindex)->second).second))->checkSeptum(*itsBunch);
    }
    
    if(((*sindex)->first) == "PROBE")    {
      (static_cast<Probe *>(((*sindex)->second).second))->checkProbe(*itsBunch, turnnumber, itsBunch->getT() * 1e9, dt);
    }
    
    if(((*sindex)->first) == "STRIPPER")    {
      bool flag_stripper = (static_cast<Stripper *>(((*sindex)->second).second))
	-> checkStripper(*itsBunch, turnnumber, itsBunch->getT() * 1e9, dt);
      if(flag_stripper) {
	itsBunch->boundp();
	*gmsg << "total particle after stripping =" << itsBunch->getTotalNum() << endl;
      }
    }
    
    if(((*sindex)->first) == "CCOLLIMATOR") {
        Collimator * collim;
	collim = static_cast<Collimator *>(((*sindex)->second).second);
	if(collim->hasSurfacePhysics()) {
	  sphys = collim->getSurfacePhysics();
	  sphys->apply(*itsBunch);
	} else {
	  collim->checkCollimator(*itsBunch, turnnumber, itsBunch->getT() * 1e9, dt);
	}
    }
    
    bool globPartOutOfBounds = (min(itsBunch->Bin) < 0);
    
    size_t lostParticleNum = 0;
    if(globPartOutOfBounds) {
      lostParticleNum=itsBunch->boundp_destroyT();
      *gmsg << "Step " << step_m << ", " << lostParticleNum << " particles lost on stripper, collimator or septum" << endl;
      flagNeedUpdate = true;
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    
    // update immediately if some particle are lost during this step
    if(flagNeedUpdate) {
      Vector_t const meanR = calcMeanR();
      Vector_t const meanP = calcMeanP();
      double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
      globalToLocal(itsBunch->R, phi, meanR);
      
      //scale coordinates
      itsBunch->R /= Vector_t(1000.0); // mm --> m
      
      // now destroy particles and update pertinent parameters in local frame
      itsBunch->boundp();
      repartition();
      
      //scale coordinates back
      itsBunch->R *= Vector_t(1000.0); // m --> mm
      
      localToGlobal(itsBunch->R, phi, meanR);
    }
  }
  return flagNeedUpdate;
}
  
void ParallelCyclotronTracker::initTrackOrbitFile() {
    std::string f = inputFileNameWithoutExtension() + string("-trackOrbit.dat");
    outfTrackOrbit_m.setf(ios::scientific, ios::floatfield);
    outfTrackOrbit_m.precision(8);
    if(myNode_m == 0) {
        if(OpalData::getInstance()->inRestartRun()) {
            outfTrackOrbit_m.open(f.c_str(), ios::app);
            outfTrackOrbit_m << "# Restart at " << step_m << " dumping step" << endl;
        } else {
            outfTrackOrbit_m.open(f.c_str());
            outfTrackOrbit_m << "# ID   x [mm]          px [rad]       y [mm]          py [rad]        z [mm]          pz [rad]" << endl;
        }
    }
}

std::string ParallelCyclotronTracker::inputFileNameWithoutExtension() {
    std::string f = OpalData::getInstance()->getInputFn();
    int const pdot = f.find(string("."), 0);
    return f.erase(pdot, f.size() - pdot);
}

void ParallelCyclotronTracker::singleParticleDump() {
    IpplTimings::startTimer(DumpTimer_m);
    if(Ippl::getNodes() > 1 ) {
        double x;
        int id;
        vector<double> tmpr;
        vector<int> tmpi;
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);

        // for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
        int found[2] = {-1, -1};
        int counter = 0;
        for(size_t i = 0; i < itsBunch->getLocalNum(); ++i) {
            if(itsBunch->ID[i] == 0) {
                found[counter] = i;
                counter++;
            }
            if(itsBunch->ID[i] == 1) {
                found[counter] = i;
                counter++;
            }
        }

        if(myNode_m == 0) {
            int notReceived = Ippl::getNodes() - 1;
            int numberOfPart = 0;
            while(notReceived > 0) {
                int node = COMM_ANY_NODE;
                Message *rmsg =  Ippl::Comm->receive_block(node, tag);
                if(rmsg == 0) ERRORMSG("Could not receive from client nodes in main." << endl);
                --notReceived;
                rmsg->get(&numberOfPart);
                for(int i = 0; i < numberOfPart; ++i) {
                    rmsg->get(&id);
                    tmpi.push_back(id);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                    rmsg->get(&x);
                    tmpr.push_back(x);
                }
                delete rmsg;
            }
            for(int i = 0; i < counter; ++i) {
                tmpi.push_back(itsBunch->ID[found[i]]);
                for(int j = 0; j < 3; ++j) {
                    tmpr.push_back(itsBunch->R[found[i]](j));
                    tmpr.push_back(itsBunch->P[found[i]](j));
                }
            }
            vector<double>::iterator itParameter = tmpr.begin();
            vector<int>::iterator itId = tmpi.begin();
            for(itId = tmpi.begin(); itId != tmpi.end(); itId++) {
                outfTrackOrbit_m << "ID" << *itId;
                for(int ii = 0; ii < 6; ii++) {
                    outfTrackOrbit_m << " " << *itParameter;
                    itParameter++;
                }
                outfTrackOrbit_m << endl;
            }
        } else {
            // for other nodes
            Message *smsg = new Message();
            smsg->put(counter);
            for(int i = 0; i < counter; ++i) {
                smsg->put(itsBunch->ID[found[i]]);
                for(int j = 0; j < 3; ++j) {
                    smsg->put(itsBunch->R[found[i]](j));
                    smsg->put(itsBunch->P[found[i]](j));
                }
            }
            if(!Ippl::Comm->send(smsg, 0, tag)) {
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
            }
        }
    } else {
        for(size_t i = 0; i < itsBunch->getLocalNum(); ++i) {
            if(itsBunch->ID[i] == 0 || itsBunch->ID[i] == 1) {
                outfTrackOrbit_m << "ID" << itsBunch->ID[i] << " ";
                outfTrackOrbit_m << itsBunch->R[i](0) << " " <<itsBunch->P[i](0) << " ";
                outfTrackOrbit_m << itsBunch->R[i](1) << " " << itsBunch->P[i](1) << " ";
                outfTrackOrbit_m << itsBunch->R[i](2) << " " << itsBunch->P[i](2) << endl; 
            }
        }
    }
    IpplTimings::stopTimer(DumpTimer_m);
}

void ParallelCyclotronTracker::evaluateSpaceChargeField() {
    Vector_t const meanR = calcMeanR();
    itsBunch->Bf = Vector_t(0.0);
    itsBunch->Ef = Vector_t(0.0);
    if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
        double const binsPhi = itsBunch->calcMeanPhi() - 0.5 * pi;
        globalToLocal(itsBunch->R, binsPhi, meanR);
        itsBunch->boundp();
        itsBunch->calcGammas_cycl();
        repartition();
        for(int b = 0; b < itsBunch->getLastemittedBin(); ++b) {
            if(itsBunch->pbin_m->getTotalNumPerBin(b) >= 1000) {
                itsBunch->setBinCharge(b, itsBunch->getChargePerParticle());
                itsBunch->computeSelfFields_cycl(b);
                INFOMSG("Bin:" << b << ", charge per particle " <<  itsBunch->getChargePerParticle() << endl);
            } else {
                INFOMSG("Note: Bin " << b << ": less than 1000 particles, omit space charge fields" << endl);
            }
        }
        itsBunch->Q = itsBunch->getChargePerParticle();
        localToGlobal(itsBunch->R, binsPhi, meanR);
        localToGlobal(itsBunch->Ef, binsPhi);
        localToGlobal(itsBunch->Bf, binsPhi);
    } else {
        Vector_t const meanP = calcMeanP();
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
        globalToLocal(itsBunch->R, phi, meanR);
        itsBunch->boundp();
        repartition();
        double const meanGamma = sqrt(1.0 + pow(meanP(0), 2.0) + pow(meanP(1), 2.0));
        itsBunch->computeSelfFields_cycl(meanGamma);
        localToGlobal(itsBunch->Ef, phi);
        localToGlobal(itsBunch->Bf, phi);
        localToGlobal(itsBunch->R, phi, meanR);
    }
}
