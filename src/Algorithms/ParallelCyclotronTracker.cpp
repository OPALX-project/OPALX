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
// $Author: adelmann, yang, winklehner $
//
// ------------------------------------------------------------------------
#include <cfloat>
#include <iostream>
#include <fstream>
#include <vector>
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ParallelCyclotronTracker.h"

#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/VariableRFCavity.h"

#include "Elements/OpalBeamline.h"
#include "Elements/OpalRing.h"

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

#include "Structure/BoundaryGeometry.h"

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
const double PIOVER180 = pi / 180.0;

Vector_t const ParallelCyclotronTracker::xaxis = Vector_t(1.0, 0.0, 0.0);
Vector_t const ParallelCyclotronTracker::yaxis = Vector_t(0.0, 1.0, 0.0);
Vector_t const ParallelCyclotronTracker::zaxis = Vector_t(0.0, 0.0, 1.0);

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
    eta_m(0.01),
    myNode_m(Ippl::myNode()),
    initialLocalNum_m(0),
    initialTotalNum_m(0),
    opalRing_m(NULL) {
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
    maxSteps_m(maxSTEPS),
    timeIntegrator_m(timeIntegrator),
    eta_m(0.01),
    myNode_m(Ippl::myNode()),
    initialLocalNum_m(bunch.getLocalNum()),
    initialTotalNum_m(bunch.getTotalNum()),
    opalRing_m(NULL) {
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
    if (opalRing_m != NULL) {
        // delete opalRing_m;
    }
}

/**
 * AAA
 *
 * @param none
 */
void ParallelCyclotronTracker::initializeBoundaryGeometry() {
  for(list<Component *>::iterator compindex = myElements.begin(); compindex != myElements.end(); compindex++) {
    bgf_m = dynamic_cast<ElementBase *>(*compindex)->getBoundaryGeometry();
    if(!bgf_m)
      continue;
    else
      break;
  }
  if (bgf_m) {
    itsDataSink->writeGeomToVtk(*bgf_m, std::string("data/testGeometry-00000.vtk"));
    OpalData::getInstance()->setGlobalGeometry(bgf_m);
    *gmsg << "* Boundary geometry initialized " << endl;
  }
}

/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::bgf_main_collision_test() {
  if(!bgf_m) return;

  Inform msg("bgf_main_collision_test ");

  /**
   *Here we check if a particles is outside the domain, flag it for deletion
   */

  Vector_t intecoords = 0.0;

  // This has to match the dT in the rk4 pusher! -DW
  //double dtime = 0.5 * itsBunch->getdT();  // Old
  double dtime = itsBunch->getdT() * getHarmonicNumber();  // New

  int triId = 0;
  size_t Nimpact = 0;
  for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
    int res = bgf_m->PartInside(itsBunch->R[i]*1.0e-3, itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId);
    if(res >= 0) {
      itsBunch->Bin[i] = -1;
      Nimpact++;
    }
  }
}


/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::openFiles(std::string SfileName) {

    std::string  SfileName2 = SfileName + std::string("-Angle0.dat");

    outfTheta0_m.precision(8);
    outfTheta0_m.setf(ios::scientific, ios::floatfield);
    outfTheta0_m.open(SfileName2.c_str());
    outfTheta0_m << "#  r [mm]      beta_r*gamma       theta [mm]      beta_theta*gamma        z [mm]          beta_z*gamma" << endl;

    SfileName2 = SfileName + std::string("-Angle1.dat");
    outfTheta1_m.precision(8);
    outfTheta1_m.setf(ios::scientific, ios::floatfield);
    outfTheta1_m.open(SfileName2.c_str());
    outfTheta1_m << "#  r [mm]      beta_r*gamma       theta [mm]      beta_theta*gamma        z [mm]          beta_z*gamma"  << endl;

    SfileName2 = SfileName + std::string("-Angle2.dat");
    outfTheta2_m.precision(8);
    outfTheta2_m.setf(ios::scientific, ios::floatfield);
    outfTheta2_m.open(SfileName2.c_str());
    outfTheta2_m << "#  r [mm]      beta_r*gamma       theta [mm]      beta_theta*gamma        z [mm]          beta_z*gamma"  << endl;

    // for single Particle Mode, output after each turn, to define matched initial phase ellipse.

    SfileName2 = SfileName + std::string("-afterEachTurn.dat");

    outfThetaEachTurn_m.precision(8);
    outfThetaEachTurn_m.setf(ios::scientific, ios::floatfield);

    outfThetaEachTurn_m.open(SfileName2.c_str());
    outfTheta2_m << "#  r [mm]      beta_r*gamma       theta [mm]      beta_theta*gamma        z [mm]          beta_z*gamma"  << endl;
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
 * @param ring
 */
void ParallelCyclotronTracker::visitOpalRing(const OpalRing &ring) {

    *gmsg << "Adding OpalRing" << endl;

    if (opalRing_m != NULL)
        delete opalRing_m;

    opalRing_m = dynamic_cast<OpalRing*>(ring.clone());

    myElements.push_back(opalRing_m);

    opalRing_m->initialise(itsBunch);

    referenceR = opalRing_m->getBeamRInit();
    referencePr = opalRing_m->getBeamPRInit();
    referenceTheta = opalRing_m->getBeamPhiInit();

    if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException("Error in ParallelCyclotronTracker::visitOpalRing",
                            "PHIINIT is out of [-180, 180)!");
    }

    referenceZ = 0.0;
    referencePz = 0.0;

    referencePtot = itsReference.getGamma() * itsReference.getBeta();
    referencePt = sqrt(referencePtot * referencePtot - referencePr * referencePr);

    if(referencePtot < 0.0)
        referencePt *= -1.0;

    sinRefTheta_m = sin(referenceTheta * PIOVER180);
    cosRefTheta_m = cos(referenceTheta * PIOVER180);

    double BcParameter[8];

    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;

    buildupFieldList(BcParameter, "OPALRING", opalRing_m);

    // Finally print some diagnostic
    *gmsg << "* Initial beam radius = " << referenceR << " [mm] " << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta = " << itsReference.getBeta() << endl;
    *gmsg << "* Total reference momentum   = " << referencePtot * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* Reference azimuthal momentum  = " << referencePt * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* Reference radial momentum     = " << referencePr * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* " << opalRing_m->getSymmetry() << " fold field symmetry "
          << endl;
    *gmsg << "* Harmonic number h= " << opalRing_m->getHarmonicNumber() << " "
          << endl;
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

    // Fresh run (no restart):
    if(!OpalData::getInstance()->inRestartRun()) {

        // Get reference values from cyclotron element
        referenceR     = elptr->getRinit();
        referenceTheta = elptr->getPHIinit();
        referenceZ     = elptr->getZinit();
        referencePr    = elptr->getPRinit();
        referencePz    = elptr->getPZinit();

        if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron", "PHIINIT is out of [-180, 180)!");
        }

        referencePtot =  itsReference.getGamma() * itsReference.getBeta();

        // Calculate reference azimuthal (tangential) momentum from total-, z- and radial momentum:
        float insqrt = referencePtot * referencePtot - referencePr * referencePr - referencePz * referencePz;

        if(insqrt < 0) {

            if(insqrt > -1.0e-10) {

	        referencePt = 0.0;

            } else {

	        throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron", "Pt imaginary!");
            }

        } else {

            referencePt = sqrt(insqrt);
        }

        if(referencePtot < 0.0)
            referencePt *= -1.0;
        // End calculate referencePt

    // Restart a run:
    } else {

        if(referenceTheta <= -180.0 || referenceTheta > 180.0) {

            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron", "PHIINIT is out of [-180, 180)!");

        }

        referencePtot =  bega;
        // Note: Nothing else has to be set, b/c everything comes from the h5 file -DW
    }

    sinRefTheta_m = sin(referenceTheta * PIOVER180);
    cosRefTheta_m = cos(referenceTheta * PIOVER180);

    *gmsg << endl;
    *gmsg << "* Bunch global starting position:" << endl;
    *gmsg << "* RINIT = " << referenceR  << " [mm]" << endl;
    *gmsg << "* PHIINIT = " << referenceTheta << " [deg]" << endl;
    *gmsg << "* ZINIT = " << referenceZ << " [mm]" << endl;
    *gmsg << endl;
    *gmsg << "* Bunch global starting momenta:" << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta = " << itsReference.getBeta() << endl;
    *gmsg << "* Reference total momentum (beta * gamma) = " << referencePtot * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference azimuthal momentum (Pt) = " << referencePt * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference radial momentum (Pr) = " << referencePr * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference axial momentum (Pz) = " << referencePz * 1000.0 << " [MCU]" << endl;
    *gmsg << endl;

    double sym = elptr->getSymmetry();
    *gmsg << "* " << sym << "-fold field symmerty " << endl;

    // ckr: this just returned the default value as defined in Component.h
    // double rff = elptr->getRfFrequ();
    // *gmsg << "* Rf frequency= " << rff << " [MHz]" << endl;

    std::string fmfn = elptr->getFieldMapFN();
    *gmsg << "* Field map file name = " << fmfn << " " << endl;

    std::string type = elptr->getType();
    *gmsg << "* Type of cyclotron = " << type << " " << endl;

    double rmin = elptr->getMinR();
    double rmax = elptr->getMaxR();
    *gmsg << "* Radial aperture = " << rmin << " ... " << rmax<<" [mm] "<< endl;

    double zmin = elptr->getMinZ();
    double zmax = elptr->getMaxZ();
    *gmsg << "* Vertical aperture = " << zmin << " ... " << zmax<<" [mm]"<< endl;

    /**
    bool Sflag = elptr->getSuperpose();
    std::string flagsuperposed;
    if (Sflag)
      flagsuperposed="yes";
    else
      flagsuperposed="no";
    *gmsg << "* Electric field maps are superposed? " << flagsuperposed << " " << endl;
    */

    double h = elptr->getCyclHarm();
    *gmsg << "* Harmonic number h = " << h << " " << endl;

    /**
    if (elptr->getSuperpose())
        *gmsg << "* Fields are superposed " << endl;
    */

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
    if(type == std::string("CARBONCYCL")) {
        fieldflag = 2;
    } else if(type == std::string("CYCIAE")) {
        fieldflag = 3;
    } else if(type == std::string("AVFEQ")) {
        fieldflag = 4;
    } else if(type == std::string("FFAG")) {
        fieldflag = 5;
    } else if(type == std::string("BANDRF")) {
        fieldflag = 6;
    } else //(type == "RING")
        fieldflag = 1;

    // read field map on the  middle plane of cyclotron.
    // currently scalefactor is set to 1.0
    // TEMP changed 1.0 to getBScale() to test if we can scale the midplane field -DW
    elptr->initialise(itsBunch, fieldflag, elptr->getBScale());

    double BcParameter[8];

    for(int i = 0; i < 8; i++)
      BcParameter[i] = 0.0;

    BcParameter[0] = elptr->getRmin();
    BcParameter[1] = elptr->getRmax();

    // Store inner radius and outer radius of cyclotron field map in the list
    buildupFieldList(BcParameter, "CYCLOTRON", elptr);
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

    *gmsg << "* --------- Collimator -----------------------------" << endl;

    Collimator* elptr = dynamic_cast<Collimator *>(coll.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXStart();
    *gmsg << "* Xstart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* Xend= " << xend << " [mm]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* Ystart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* Yend= " << yend << " [mm]" << endl;

    double zstart = elptr->getZStart();
    *gmsg << "* Zstart= " << zstart << " [mm]" << endl;

    double zend = elptr->getZEnd();
    *gmsg << "* Zend= " << zend << " [mm]" << endl;

    double width = elptr->getWidth();
    *gmsg << "* Width= " << width << " [mm]" << endl;

    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    std::string ElementType = "CCOLLIMATOR";
    BcParameter[0] = xstart ;
    BcParameter[1] = xend;
    BcParameter[2] = ystart ;
    BcParameter[3] = yend;
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
 * @param degrader
 */
void ParallelCyclotronTracker::visitDegrader(const Degrader &deg) {
    *gmsg << "In Degrader; L= " << deg.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Degrader *>(deg.clone()));

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

void ParallelCyclotronTracker::visitOffset(const Offset & off) {
    if (opalRing_m == NULL)
        throw OpalException(
                        "ParallelCylcotronTracker::visitOffset",
                        "Attempt to place an offset when OpalRing not defined");
    Offset* offNonConst = const_cast<Offset*>(&off);
    offNonConst->updateGeometry(opalRing_m->getNextPosition(),
                       opalRing_m->getNextNormal());
    opalRing_m->appendElement(off);
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
    *gmsg << "* -----------  Probe -------------------------------" << endl;
    Probe *elptr = dynamic_cast<Probe *>(prob.clone());
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


    // initialise, do nothing
    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    std::string ElementType = "PROBE";
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

void ParallelCyclotronTracker::visitSBend3D(const SBend3D &bend) {
    *gmsg << "Adding SBend3D" << endl;
    if (opalRing_m != NULL)
        opalRing_m->appendElement(bend);
    else
        throw OpalException("ParallelCyclotronTracker::visitSBend3D",
                      "Need to define a RINGDEFINITION to use SBend3D element");
}

void ParallelCyclotronTracker::visitVariableRFCavity(const VariableRFCavity &cav) {
    *gmsg << "Adding Variable RF Cavity" << endl;
    if (opalRing_m != NULL)
        opalRing_m->appendElement(cav);
    else
        throw OpalException("ParallelCyclotronTracker::visitVariableRFCavity",
            "Need to define a RINGDEFINITION to use VariableRFCavity element");
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

    std::string fmfn = elptr->getFieldMapFN();
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
    std::string ElementType = "CAVITY";
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
    *gmsg << "In RFQuadrupole; L = " << rfq.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RFQuadrupole *>(rfq.clone()));
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitSBend(const SBend &bend) {
    *gmsg << "In SBend; L = " << bend.getElementLength() << " however the element is missing " << endl;
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

    *gmsg << endl << "* -----------------------------  Septum ------------------------------- *" << endl;

    Septum *elptr = dynamic_cast<Septum *>(sept.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXstart();
    *gmsg << "XStart = " << xstart << " [mm]" << endl;

    double xend = elptr->getXend();
    *gmsg << "XEnd = " << xend << " [mm]" << endl;

    double ystart = elptr->getYstart();
    *gmsg << "YStart = " << ystart << " [mm]" << endl;

    double yend = elptr->getYend();
    *gmsg << "YEnd = " << yend << " [mm]" << endl;

    double width = elptr->getWidth();
    *gmsg << "Width = " << width << " [mm]" << endl;


    // initialise, do nothing
    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    std::string ElementType = "SEPTUM";
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
    *gmsg << "Charge of outcoming particle = +e * " << opcharge << endl;

    double opmass = elptr->getOPMass();
    *gmsg << "* Mass of the outcoming particle = " << opmass << " [GeV/c^2]" << endl;

    elptr->initialise(itsBunch, 1.0);

    double BcParameter[8];
    for(int i = 0; i < 8; i++)
        BcParameter[i] = 0.0;
    std::string ElementType = "STRIPPER";
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
void ParallelCyclotronTracker::buildupFieldList(double BcParameter[], std::string ElementType, Component *elptr) {
    beamline_list::iterator sindex;

    type_pair *localpair = new type_pair();
    localpair->first = ElementType;

    for(int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;

    // always put cyclotron as the first element in the list.
    if(ElementType == "OPALRING" || ElementType == "CYCLOTRON") { // CYCLOTRON is TEMP -DW
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);

}

/**
 *
 *
 * @param bl
 */
void ParallelCyclotronTracker::visitBeamline(const Beamline &bl) {
    itsBeamline->iterate(*dynamic_cast<BeamlineVisitor *>(this), false);
}

void ParallelCyclotronTracker::checkNumPart(std::string s) {
    int nlp = itsBunch->getLocalNum();
    int minnlp = 0;
    int maxnlp = 111111;
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

    // Record how many bunches have already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();

    // For the time being, we set bin number equal to bunch number. FixMe: not used
    BinCount_m = BunchCount_m;

    itsBeamline->accept(*this);
    if (opalRing_m != NULL)
        opalRing_m->lockRing();

    // Display the selected elements
    *gmsg << "* -------------------------------------" << endl;
    *gmsg << "* The selected Beam line elements are :" << endl;

    for(beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++)
        *gmsg << "* -> " <<  ((*sindex)->first) << endl;

    *gmsg << "* -------------------------------------" << endl;

    // Don't initializeBoundaryGeometry()
    // Get BoundaryGeometry that is already initialized
    bgf_m = OpalData::getInstance()->getGlobalGeometry();

    // External field arrays for dumping
    for(int k = 0; k < 2; k++)
        FDext_m[k] = Vector_t(0.0, 0.0, 0.0);

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    if(timeIntegrator_m == 0) {
        *gmsg << "* 4th order Runge-Kutta integrator" << endl;
#ifdef GENERICTRACKER
        Tracker_Generic();
#else
        Tracker_RK4();
#endif
    } else if(timeIntegrator_m == 1) {
        *gmsg << "* 2nd order Leap-Frog integrator" << endl;
#ifdef GENERICTRACKER
        Tracker_Generic();
#else
        Tracker_LF();
#endif
    } else if(timeIntegrator_m == 2) {
        *gmsg << "* Multiple time stepping (MTS) integrator" << endl;
        Tracker_MTS();
    } else {
        *gmsg << "ERROR: Invalid name of TIMEINTEGRATOR in Track command" << endl;
        exit(1);
    }
    *gmsg << "* ----------------------------------------------- *" << endl;
    *gmsg << "* Finalizing i.e. write data and close files :" << endl;
    for(beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++) {
        (((*sindex)->second).second)->finalise();
    }
    *gmsg << "* ----------------------------------------------- *" << endl;
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

    BorisPusher pusher;

    // time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();

    const double harm = getHarmonicNumber();

    // load time
    const double dt = itsBunch->getdT() * 1.0e9 * harm; //[s]-->[ns]

    // find the injection time interval
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << stepsPerTurn *dt << "[ns]" << endl;
    }

    initTrackOrbitFile();

    int SteptoLastInj = itsBunch->getSteptoLastInj();

    // get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun()) {
        restartStep0_m = itsBunch->getLocalTrackStep();
        step_m = restartStep0_m;
        if (numBunch_m > 1) itsBunch->resetPartBinID2(eta_m);
        *gmsg << "* Restart at integration step " << restartStep0_m << endl;
    }

    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        itsBunch->setT(0.0);
    }

    *gmsg << "* Beginning of this run is at t= " << itsBunch->getT() * 1e9 << " [ns]" << endl;
    *gmsg << "* The time step is set to dt= " << dt << " [ns]" << endl;

    // for single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1)
        openFiles(OpalData::getInstance()->getInputBasename());

    double const initialReferenceTheta = referenceTheta * PIOVER180;

    initDistInGlobalFrame();

    //  read in some control parameters
    const int SinglePartDumpFreq = Options::sptDumpFreq;
    const int resetBinFreq = Options::rebinFreq;
    const int scSolveFreq = Options::scSolveFreq;
    const bool doDumpAfterEachTurn = Options::psDumpEachTurn;
    const int boundpDestroyFreq = Options::boundpDestroyFreq;

    // prepare for dump after each turn
    double oldReferenceTheta = initialReferenceTheta;

    *gmsg << "* Single particle trajectory dump frequency is set to " << SinglePartDumpFreq << endl;
    *gmsg << "* Repartition frequency is set to " << Options::repartFreq << endl;
    if(numBunch_m > 1)
        *gmsg << "particles energy bin ID reset frequency is set to " << resetBinFreq << endl;

    // if initialTotalNum_m = 2, trigger SEO mode
    // prepare for transverse tuning calculation
    vector<double> Ttime, Tdeltr, Tdeltz;
    // prepare for transverse tuning calculation
    vector<int> TturnNumber;
    turnnumber_m = 1;


    // flag to determine when to transit from single-bunch to multi-bunches mode
    bool flagTransition = false;
    // step point determining the next time point of check for transition
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();

    const  double deltaTheta = pi / (stepsPerTurn);
    // record at which angle the space charge are solved
    double angleSpaceChargeSolve = 0.0;

    if(initialTotalNum_m == 1) {
        *gmsg << endl;
        *gmsg << "* *---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," << endl
              << "* The initial distribution file must be specified which should contain only one line for the single particle " << endl
              << "* *------------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
	*gmsg << endl;
        *gmsg << "* *------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." << endl
              << "* This mode does NOT include any RF cavities. The initial distribution file must be specified" << endl
              << "* In the file the first line is for reference particle and the second line is for offcenter particle." << endl
              << "* The tunes are calculated by FFT routines based on these two particles. " << endl
              << "* *------------NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    }

    // apply the plugin elements: probe, collimator, stripper, septum
    // make sure that we apply elements even on first step
    applyPluginElements(dt);

    // *****************II***************
    // main integration loop
    // *****************II***************
    *gmsg << "* ---------------------------- Start tracking ----------------------------" << endl;
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

                    itsBunch->R *= Vector_t(0.001); // mm --> m
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

                        *gmsg << "*** Save beam distribution at turn #" << turnnumber_m << " ***" << endl;
                        *gmsg << "*** After one revolution, Multi-Bunch Mode will be invoked ***" << endl;

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
                        itsBunch->resetPartBinID2(eta_m);
                    } else if(multiBunchMode_m == 2) {

                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        itsBunch->resetPartBinID2(eta_m);
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
                    itsBunch->R *= Vector_t(0.001); // mm --> m

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
                    itsBunch->R *= Vector_t(0.001); // mm --> m

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
                itsBunch->R *= Vector_t(0.001); // mm --> m

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

            externalB = Vector_t(0.0, 0.0, 0.0);
            externalE = Vector_t(0.0, 0.0, 0.0);

            beamline_list::iterator sindex = FieldDimensions.begin();
            (((*sindex)->second).second)->apply(i, itsBunch->getT() * 1e9, externalE, externalB);
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

        // apply the plugin elements: probe, colilmator, stripper, septum
        applyPluginElements(dt);
        // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global aperture
        bool flagNeedUpdate = deleteParticle();
        if(itsBunch->weHaveBins() && flagNeedUpdate)
          itsBunch->resetPartBinID2(eta_m);

        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % resetBinFreq == 0)
            itsBunch->resetPartBinID2(eta_m);

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
                ++turnnumber_m;

                *gmsg << endl;
                *gmsg << "*** Turn " << turnnumber_m << endl;

                dumpEachTurn = true;
                outfThetaEachTurn_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfThetaEachTurn_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                                    << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                                    << " " << temp_meanTheta / PIOVER180
                                    << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                                    << " " << variable_m[2]
                                    << " " << variable_m[5] << endl;
            }
            // FixMe: should be defined elesewhere !
            // define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
            const double azimuth_angle0 = 0.0;
            const double azimuth_angle1 = 22.5 * PIOVER180;
            const double azimuth_angle2 = 45.0 * PIOVER180;
            if((oldReferenceTheta < azimuth_angle0 - deltaTheta) && (temp_meanTheta >= azimuth_angle0 - deltaTheta)) {
                outfTheta0_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / PIOVER180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                outfTheta1_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / PIOVER180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                outfTheta2_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / PIOVER180
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
                ++turnnumber_m;
                dumpEachTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch->getTotalNum() << endl;
            }
        }

        // dump phase space distribution of bunch
        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2) ||
           (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {

            IpplTimings::startTimer(DumpTimer_m);

            itsBunch->setSteptoLastInj(SteptoLastInj);
            itsBunch->setLocalTrackStep((step_m + 1));

            // Write Phase Space and Statistics Data to h5 and dat files
            bunchDumpPhaseSpaceStatData();

            // extE_m = Vector_t(0.0, 0.0, 0.0);
            // extB_m = Vector_t(0.0, 0.0, 0.0);

            // //--------------------- calculate mean coordinates  of bunch -------------------------------//
            // //------------  and calculate the external field at the mass of bunch-----------------------//

            // Vector_t const meanR = calcMeanR();
	    // // *gmsg << "* meanR=( " << meanR(0) << " " << meanR(1) << " " << meanR(2) << " ) [mm] " << endl;

            // beamline_list::iterator DumpSindex = FieldDimensions.begin();
            // (((*DumpSindex)->second).second)->apply(meanR, Vector_t(0.0), itsBunch->getT() * 1e9, extE_m, extB_m);
            // FDext_m[0] = extB_m / 10.0; // kgauss -> T
            // FDext_m[1] = extE_m;

            // // --------------------------- Particle dumping ------------------------------------------ //
            // // Note: Don't dump when
            // // 1. after one turn
            // // in order to sychronize the dump step for multi-bunch and single bunch for compare
            // // with each other during post-process phase.
            // double const E = itsBunch->get_meanEnergy();

            // if(!(Options::psDumpLocalFrame)) {
            //     // --------------------------- Dump in global frame ---------------------------------- //

            //     itsBunch->R /= Vector_t(1000.0); // mm --> m

            //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m, E, referencePr, referenceR, referenceTheta);
            //     itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);

            //     itsBunch->R *= Vector_t(1000.0); // m --> mm

            //     *gmsg << "* Phase space dump " << lastDumpedStep_m << " (global frame) at integration step "
            //           << step_m + 1 << " T = " << itsBunch->getT() * 1e9 << " [ns]" 	    << " E = " << itsBunch->get_meanEnergy()  << endl;

            // } else {
            //     // --------------------------- Dump in local frame ---------------------------------- //

	    //     Vector_t const meanP = calcMeanP();
	    //     double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
	    //     Vektor<double, 4> quaternionToXAxis;
            //     Vector_t const xaxis = Vector_t(1.0, 0.0, 0.0);

            //     getQuaternionTwoVectors(meanP, xaxis, quaternionToXAxis);

            //     globalToLocal(itsBunch->R, quaternionToXAxis, meanR);
            //     globalToLocal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //globalToLocal(itsBunch->R, phi, meanR);
            //     //globalToLocal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

	    //     itsBunch->R /= Vector_t(1000.0); // mm --> m

	    //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m, E, referencePr, referenceR, referenceTheta);
	    //     itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0, E);

	    //     itsBunch->R *= Vector_t(1000.0); // m --> mm

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
            //     localToGlobal(itsBunch->R, quaternionToXAxis, meanR);
            //     localToGlobal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //localToGlobal(itsBunch->R, phi, meanR);
            //     //localToGlobal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

	    //     *gmsg << "* Phase space dump " << lastDumpedStep_m << " (local frame) at integration step "
	    // 	      << step_m + 1 << " T = " << itsBunch->getT() * 1e9 << " [ns]" << " E = " << E  << " phi= " << phi/pi*180.0 << endl;
            // }
            // // -------------------- Done with particle dumping ------------------------------------ //
            IpplTimings::stopTimer(DumpTimer_m);
        }
    }

    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 0) {

            double FinalMomentum2  = pow(itsBunch->P[ii](0), 2.0) +
                                     pow(itsBunch->P[ii](1), 2.0) +
                                     pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " MeV" << endl;
            *gmsg << "* Total phase space dump number including the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
        }
    }

    Ippl::Comm->barrier();

    if(myNode_m == 0) outfTrackOrbit_m.close();

    if(initialTotalNum_m == 1)
        closeFiles();

    *gmsg << *itsBunch << endl;
}

void ParallelCyclotronTracker::Tracker_RK4() {

    // time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();

    // record how many bunches have already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();

    // decide how many energy bins. ONLY FOR MPM
    BinCount_m = BunchCount_m;

    const double harm = getHarmonicNumber();

    // load time
    double t  = itsBunch->getT() * 1.0e9;
    const double dt = itsBunch->getdT() * 1.0e9 * harm; //[s]-->[ns]

    // Find the injection time interval
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << stepsPerTurn * dt << "[ns]" << endl;
    }

    initTrackOrbitFile();

    // Get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun()) {
        restartStep0_m = itsBunch->getLocalTrackStep();
        step_m = restartStep0_m;

        if (numBunch_m > 1) itsBunch->resetPartBinID2(eta_m);
        *gmsg << "* Restart at integration step " << restartStep0_m << endl;
    }

    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        t = 0.0;
    }

    *gmsg << "* Beginning of this run is at t = " << t << " [ns]" << endl;
    *gmsg << "* The time step is set to dt = " << dt << " [ns]" << endl;

    // For single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1)
        openFiles(OpalData::getInstance()->getInputBasename());

    initDistInGlobalFrame(); // AAA

    //  Read in some control parameters
    const int SinglePartDumpFreq = Options::sptDumpFreq;
    const int resetBinFreq = Options::rebinFreq;
    const int scSolveFreq = Options::scSolveFreq;
    const bool doDumpAfterEachTurn = Options::psDumpEachTurn;
    const int boundpDestroyFreq = Options::boundpDestroyFreq;

    // Set up some angles in the midplane
    const double initialReferenceTheta = referenceTheta / 180.0 * pi;
    double oldReferenceTheta = initialReferenceTheta; // init here, reset each timestep
    const double deltaTheta = pi / stepsPerTurn;      // half of the average angle per step

    *gmsg << "* Single particle trajectory dump frequency is set to " << SinglePartDumpFreq << endl;
    *gmsg << "* The frequency to solve space charge fields is set to " << scSolveFreq << endl;
    *gmsg << "* The repartition frequency is set to " << Options::repartFreq << endl;

    if(numBunch_m > 1)
        *gmsg << "* The particles energy bin reset frequency is set to " << resetBinFreq << endl;

    // If initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
    vector<double> Ttime, Tdeltr, Tdeltz;
    vector<int> TturnNumber;
    turnnumber_m = 1;

    bool flagNoDeletion = false;

    // Flag to determine when to transit from single-bunch to multi-bunches mode
    bool flagTransition = false;

    // Step point determining the next time point of check for transition
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();

    if(initialTotalNum_m == 1) {
        *gmsg << endl;
        *gmsg << "* ---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," << endl
              << "* The initial distribution file must be specified which should contain only one line for the single particle " << endl
              << "* ------------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
        *gmsg << endl;
        *gmsg << "* ------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** " << endl;
        *gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." << endl
              << "* This mode does NOT include any RF cavities. The initial distribution file must be specified" << endl
              << "* In the file the first line is for reference particle and the second line is for offcenter particle." << endl
              << "* The tune is calculated by FFT routines based on these two particles. " << endl
              << "* ------------NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** " << endl;
        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    }


    // apply the plugin elements: probe, collimator, stripper, septum
    // make sure that we apply elements even on first step
    applyPluginElements(dt);

    // main integration loop
    *gmsg << "* ---------------------------- Start tracking ----------------------------" << endl;
    for(; step_m < maxSteps_m; step_m++) {
	//*gmsg << "step_m= " << step_m << endl;
        /*
        Vector_t B(0., 0., 0.), E(0., 0., 0.);
        Vector_t Pt = itsBunch->P[0];
        (*FieldDimensions.begin())->second.second->apply(0, itsBunch->getT()*1e9, E, B);
        double PTot = sqrt(Pt(0)*Pt(0)+Pt(1)*Pt(1)+Pt(2)*Pt(2));
        std::cerr << step_m << " P: " << PTot << " ";
        std::cerr << "R: ";
        for (int jj = 0; jj < 3; jj++)
            std::cerr << itsBunch->R[0](jj) << " ";
        std::cerr << "B: ";
        for (int jj = 0; jj < 3; jj++)
            std::cerr << B(jj) << " ";
        std::cerr << "E: ";
        for (int jj = 0; jj < 3; jj++)
            std::cerr << E(jj) << " ";
        std::cerr << "P: ";
        for (int jj = 0; jj < 3; jj++)
            std::cerr << itsBunch->P[0](jj) << " ";
        std::cerr << std::endl;
        */
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
                        for(int ii = 0; ii < 1; ii++) {
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
                            for(int ii = 0; ii < 12; ii++) {
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

                        itsBunch->R *= Vector_t(0.001); // mm --> m
                        itsBunch->calcBeamParameters_cycl();
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm

                        RThisTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));

                        Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

                        double XYrms =  sqrt(pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0));

                        // if the distance between two neighbour bunch is less than 5 times of its 2D rms size
                        // start multi-bunch simulation, fill current phase space to initialR and initialP arrays

                        if((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                            // since next turn, start multi-bunches
                            saveOneBunch();
                            flagTransition = true;
                            *gmsg << "*** Save beam distribution at turn #" << turnnumber_m << " ***" << endl;
                            *gmsg << "*** After one revolution, Multi-Bunch Mode will be invoked ***" << endl;

                        }
                        stepsNextCheck += stepsPerTurn;

                        *gmsg << "* RLastTurn = " << RLastTurn_m << " [mm]" << endl;
                        *gmsg << "* RThisTurn = " << RThisTurn_m << " [mm]" << endl;
                        *gmsg << "* XYrms = " << XYrms    << " [mm]" << endl;

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

                        //itsBunch->resetPartBinID2(eta_m);
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

            Vector_t const meanR = calcMeanR();
            oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

            // Calculate SC field before each time step and keep constant during integration.
            // Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
            if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {

                if(step_m % scSolveFreq == 0) {

                    // Firstly reset E and B to zero before fill new space charge field data for each track step
                    itsBunch->Bf = Vector_t(0.0);
                    itsBunch->Ef = Vector_t(0.0);

                    IpplTimings::startTimer(TransformTimer_m);

                    PreviousMeanP = calcMeanP();

                    // // in global Cartesian frame, calculate the location in global frame of bunch
                    // oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

                    if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
                        // --- Multibunche mode --- //

                        // Since calcMeanP takes into account all particles of all bins (TODO: Check this! -DW)
                        // Using the quaternion method with PreviousMeanP and yaxis should give the correct result

			Quaternion_t quaternionToYAxis;

                        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

                        globalToLocal(itsBunch->R, quaternionToYAxis, meanR);

                        itsBunch->R *= Vector_t(0.001); // mm --> m

                        // double binsMeanPhi = itsBunch->calcMeanPhi() - 0.5 * pi;
                        // angleSpaceChargeSolve = binsMeanPhi;

                        // double cosTemp_binsMeanPhi = cos(binsMeanPhi);
                        // double sinTemp_binsMeanPhi = sin(binsMeanPhi);

                        // // remove mean coordinates
                        // itsBunch->R -= meanR;

                        // //scale coordinates
                        // itsBunch->R /= Vector_t(1000.0); // mm --> m

                        // // rotate from global frame to local frame(transverse horizontal,longitudinal,transverse vertical)
                        // // For multi-bin, rotate the frame for binsMeanPhi degree
                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_binsMeanPhi
                        //                                  + itsBunch->R[ii](1) * sinTemp_binsMeanPhi;
                        //     double temp_RLongitudinal = -itsBunch->R[ii](0) * sinTemp_binsMeanPhi
                        //                                 + itsBunch->R[ii](1) * cosTemp_binsMeanPhi;

                        //     itsBunch->R[ii](0) = temp_RHorizontal;
                        //     itsBunch->R[ii](1) = temp_RLongitudinal;
                        // }

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        // Calcualte gamma for each energy bin
                        itsBunch->calcGammas_cycl();

                        repartition();

                        // calculate space charge field for each energy bin
                        // TODO update this to include rotation of the BoundaryGeometry -Dw
                        for(int b = 0; b < itsBunch->getLastemittedBin(); b++) {

                            itsBunch->setBinCharge(b, itsBunch->getChargePerParticle());
                            itsBunch->computeSelfFields_cycl(b);
                        }

                        itsBunch->Q = itsBunch->getChargePerParticle();

                        IpplTimings::startTimer(TransformTimer_m);

                        //scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // Transform coordinates back to global
                        localToGlobal(itsBunch->R, quaternionToYAxis, meanR);

                        // Transform self field back to global frame (rotate only)
                        localToGlobal(itsBunch->Ef, quaternionToYAxis);
                        localToGlobal(itsBunch->Bf, quaternionToYAxis);

                        // // HERE transform particles coordinates back to global frame (rotate and shift)
                        // // rotate back from local reference frame to global frame
                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_binsMeanPhi
                        //                                  - itsBunch->R[ii](1) * sinTemp_binsMeanPhi;
                        //     double temp_RLongitudinal =  itsBunch->R[ii](0) * sinTemp_binsMeanPhi
                        //                                  + itsBunch->R[ii](1) * cosTemp_binsMeanPhi;

                        //     itsBunch->R[ii](0) = temp_RHorizontal;
                        //     itsBunch->R[ii](1) = temp_RLongitudinal;
                        // }

                        // //scale coordinates back
                        // itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // // retrieve mean coordinates
                        // itsBunch->R += meanR;

                        // HERE transform self field back to global frame (rotate)

                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_binsMeanPhi
                        //                       - itsBunch->Ef[ii](1) * sinTemp_binsMeanPhi;
                        //     double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_binsMeanPhi
                        //                       + itsBunch->Ef[ii](1) * cosTemp_binsMeanPhi;

                        //     itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                        //     itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                        // }

                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_binsMeanPhi
                        //                       - itsBunch->Bf[ii](1) * sinTemp_binsMeanPhi;
                        //     double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_binsMeanPhi
                        //                       + itsBunch->Bf[ii](1) * cosTemp_binsMeanPhi;

                        //     itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                        //     itsBunch->Bf[ii](1) = temp_E2;  // By,T
                        // }
                    } else {
                        // --- Single bunch mode --- //
                        double temp_meangamma = sqrt(1.0 + dot(PreviousMeanP, PreviousMeanP));

	                Quaternion_t quaternionToYAxis;

                        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

                        globalToLocal(itsBunch->R, quaternionToYAxis, meanR);

                        itsBunch->R *= Vector_t(0.001); // mm --> m

                        // // in global Cartesian frame, calculate the direction of longitudinal angle of bunch
                        // double meanPhi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

                        // angleSpaceChargeSolve = meanPhi;

                        // double cosTemp_meanPhi = cos(meanPhi);
                        // double sinTemp_meanPhi = sin(meanPhi);

                        // double meanPLongitudinal2 = pow(meanP(0), 2.0) + pow(meanP(1), 2.0);
                        // double temp_meangamma = sqrt(1.0 + meanPLongitudinal2);

                        // // remove mean coordinates
                        // itsBunch->R -= meanR;

                        //scale coordinates
                        //itsBunch->R /= Vector_t(1000.0); // mm --> m

                        // // rotate from global frame to local frame(transverse horizontal,longitudinal, transverse vertical)
                        // // For single bin, rotate the frame for meanPhi degree
                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_meanPhi
                        //                                  + itsBunch->R[ii](1) * sinTemp_meanPhi;
                        //     double temp_RLongitudinal = -itsBunch->R[ii](0) * sinTemp_meanPhi
                        //                                 + itsBunch->R[ii](1) * cosTemp_meanPhi;

                        //     itsBunch->R[ii](0) = temp_RHorizontal;
                        //     itsBunch->R[ii](1) = temp_RLongitudinal;
                        // }

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        repartition();

                        itsBunch->computeSelfFields_cycl(temp_meangamma, 0.001 * meanR, quaternionToYAxis);

                        IpplTimings::startTimer(TransformTimer_m);

                        // HERE transform particles coordinates back to global frame (rotate and shift)
                        // rotate back from local reference frame to global frame
                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_RHorizontal   =  itsBunch->R[ii](0) * cosTemp_meanPhi
                        //                                  - itsBunch->R[ii](1) * sinTemp_meanPhi;
                        //     double temp_RLongitudinal =  itsBunch->R[ii](0) * sinTemp_meanPhi
                        //                                  + itsBunch->R[ii](1) * cosTemp_meanPhi;

                        //     itsBunch->R[ii](0) = temp_RHorizontal;
                        //     itsBunch->R[ii](1) = temp_RLongitudinal;
                        // }

                        //scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // Transform coordinates back to global
                        localToGlobal(itsBunch->R, quaternionToYAxis, meanR);

                        // Transform self field back to global frame (rotate only)
                        localToGlobal(itsBunch->Ef, quaternionToYAxis);
                        localToGlobal(itsBunch->Bf, quaternionToYAxis);

                        // retrieve mean coordinates
                        //itsBunch->R += meanR;

                        // Transform self field back to global frame (rotate only)
                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_meanPhi  - itsBunch->Ef[ii](1) * sinTemp_meanPhi;
                        //     double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_meanPhi  + itsBunch->Ef[ii](1) * cosTemp_meanPhi;
                        //     itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                        //     itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                        // }

                        // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                        //     double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_meanPhi  - itsBunch->Bf[ii](1) * sinTemp_meanPhi;
                        //     double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_meanPhi  + itsBunch->Bf[ii](1) * cosTemp_meanPhi;
                        //     itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                        //     itsBunch->Bf[ii](1) = temp_E2;  // By,T
                        // }
                    }

                    IpplTimings::stopTimer(TransformTimer_m);

                } else {
		    // If we are not solving for the space charge fields at this time step
                    // we will apply the fields from the previous step and have to rotate them
                    // accordingly. For this we find the quaternion between the previous mean momentum (PreviousMeanP)
                    // and the current mean momentum (meanP) and rotate the fields with this quaternion.

                    // Transform particles coordinates to local frame (rotate and shift)
                    //Vector_t const meanR = calcMeanR();
                    Vector_t const meanP = calcMeanP();

                    Quaternion_t quaternionToNewMeanP;

                    getQuaternionTwoVectors(PreviousMeanP, meanP, quaternionToNewMeanP);

                    // Reset PreviousMeanP. Cave: This HAS to be after the quaternion is calculated!
                    PreviousMeanP = calcMeanP();

                    // Rotate the fields
                    globalToLocal(itsBunch->Ef, quaternionToNewMeanP);
                    globalToLocal(itsBunch->Bf, quaternionToNewMeanP);

                    // // in global Cartesian frame, calculate the location in global frame of bunch
                    //oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

                    // // in global Cartesian frame, calculate the direction of longitudinal angle of bunch
                    // double meanPhi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

                    // double deltaPhi = meanPhi - angleSpaceChargeSolve;

                    // double cosTemp_deltaPhi = cos(deltaPhi);
                    // double sinTemp_deltaPhi = sin(deltaPhi);


                    // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                    //     double temp_E1 =  itsBunch->Ef[ii](0) * cosTemp_deltaPhi  - itsBunch->Ef[ii](1) * sinTemp_deltaPhi;
                    //     double temp_E2 =  itsBunch->Ef[ii](0) * sinTemp_deltaPhi  + itsBunch->Ef[ii](1) * cosTemp_deltaPhi;
                    //     itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
                    //     itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
                    // }

                    // for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
                    //     double temp_E1 =  itsBunch->Bf[ii](0) * cosTemp_deltaPhi  - itsBunch->Bf[ii](1) * sinTemp_deltaPhi;
                    //     double temp_E2 =  itsBunch->Bf[ii](0) * sinTemp_deltaPhi  + itsBunch->Bf[ii](1) * cosTemp_deltaPhi;
                    //     itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
                    //     itsBunch->Bf[ii](1) = temp_E2;  // By,T
                }
            }
            // If field solver is not available or number of particles is less than 1000
            // we don't need to do anything space charge related.

            // ------------------------ Track all particles for one step --------------------------- //
            IpplTimings::startTimer(IntegrationTimer_m);

            // *** This was moved here b/c collision should be tested before the **********************
            // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
            // apply the plugin elements: probe, collimator, stripper, septum
            applyPluginElements(dt);  // New

	    // check if we loose particles at the boundary
	    bgf_main_collision_test();  // New
            // ****************************************************************************************

            for(size_t i = 0; i < (itsBunch->getLocalNum()); i++) {
                flagNoDeletion = true;
                // change phase space parameters from localframe of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
                for(int j = 0; j < 3; j++) {
                    variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
                    variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
                    rold_m[j] = variable_m[j]; // used for gap cross checking
                    pold_m[j] = variable_m[j+3]; // used for gap cross
                }
                // integrate for one step in the lab Cartesian frame (absolute value).
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

            // *** This has to go before the actual timestep -DW ******************************
            // apply the plugin elements: probe, collimator, stripper, septum
            //applyPluginElements(dt);  // Old

	    // check if we loose particles at the boundary
	    //bgf_main_collision_test();  // Old
            // ********************************************************************************

            // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global aperture
            bool flagNeedUpdate = deleteParticle();

            if(itsBunch->weHaveBins() && flagNeedUpdate)
              itsBunch->resetPartBinID2(eta_m);

	    // recalculate bingamma and reset the BinID for each particles according to its current gamma
	    if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % resetBinFreq == 0)
	      itsBunch->resetPartBinID2(eta_m);

            if((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {
                ++turnnumber_m;
                dumpEachTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch->getTotalNum() << endl;
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

                // integrate for one step in the lab Cartesian frame (absulate value ).
                rk4(variable_m, t, dt, i);

		if( (i == 0) && (step_m > 10) && ((step_m%stepsPerTurn) == 0))   ++turnnumber_m;

                for(int j = 0; j < 3; j++) itsBunch->R[i](j) = variable_m[j] ; //[x,y,z]  units: [mm]
                for(int j = 0; j < 3; j++) itsBunch->P[i](j) = variable_m[j+3] ; //[px,py,pz]  units: dimensionless, beta*gama

            }//end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode

            // store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
            if(step_m % SinglePartDumpFreq == 0) {
                Ttime.push_back(t * 1.0e-9);
                Tdeltz.push_back(z_tuning[1]);
                Tdeltr.push_back(r_tuning[1]  - r_tuning[0]);
                TturnNumber.push_back(turnnumber_m);
            }
        }
        else if(initialTotalNum_m == 1) {
          // initialTotalNum_m == 1 trigger single particle mode

          IpplTimings::startTimer(IntegrationTimer_m);
          flagNoDeletion = true;

          // *** This was moved here b/c collision should be tested before the **********************
          // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
          // apply the plugin elements: probe, collimator, stripper, septum
          applyPluginElements(dt);  // New

	  // check if we loose particles at the boundary
	  bgf_main_collision_test();  // New
          // ****************************************************************************************

          // track for one step
          for( unsigned int i = 0; i < itsBunch->getLocalNum(); i++) {

            if((step_m % SinglePartDumpFreq == 0)) {
              outfTrackOrbit_m << "ID" <<itsBunch->ID[i];
              outfTrackOrbit_m << " " << itsBunch->R[i](0) << " " <<itsBunch->P[i](0) << " " <<itsBunch->R[i](1)
                               << " " << itsBunch->P[i](1) << " " <<itsBunch->R[i](2) << " " <<itsBunch->P[i](2)<< endl;
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
              ++turnnumber_m;
              dumpEachTurn = true;
              *gmsg << "Turn " << turnnumber_m << endl;

              outfThetaEachTurn_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
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
              outfTheta0_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
              outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                           << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                           << " " << temp_meanTheta / pi * 180
                           << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                           << " " << variable_m[2]
                           << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
              outfTheta1_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
              outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                           << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                           << " " << temp_meanTheta / pi * 180
                           << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                           << " " << variable_m[2]
                           << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
              outfTheta2_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
              outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                           << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                           << " " << temp_meanTheta / pi * 180
                           << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                           << " " << variable_m[2]
                           << " " << variable_m[5] << endl;
            }

            oldReferenceTheta = temp_meanTheta;

            // integrate for one step in the lab Cartesian frame (absolute value).
            flagNoDeletion = rk4(variable_m, t, dt, i);

            if(!flagNoDeletion) {
              *gmsg << "particle" << "is lost at " << step_m << "th step!" << endl;
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
          }

          // *** This has to go before the actual timestep -DW ******************************
          // apply the plugin elements: probe, collimator, stripper, septum
          //applyPluginElements(dt);  // Old

          // check if we loose particles at the boundary
          //bgf_main_collision_test();  // Old
          // ********************************************************************************

          // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global apeture
          deleteParticle();

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

        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2)
           || (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {

            IpplTimings::startTimer(DumpTimer_m);

            itsBunch->setLocalTrackStep((step_m + 1));

            // Write Phase Space and Statistics Data to h5 and dat files
            bunchDumpPhaseSpaceStatData();

            // extE_m = Vector_t(0.0, 0.0, 0.0);
            // extB_m = Vector_t(0.0, 0.0, 0.0);

            // //--------------------- calculate mean coordinates  of bunch -------------------------------//
            // //------------  and calculate the external field at the mass of bunch-----------------------//

            // Vector_t meanR = calcMeanR();

            // // define longitudinal direction of the bunch
            // // x:[0] transverse horizontal, y:[1] longitudinal, z:[2] transverse vertical

            // beamline_list::iterator DumpSindex = FieldDimensions.begin();
            // (((*DumpSindex)->second).second)->apply(meanR, Vector_t(0.0), t, extE_m, extB_m);
            // FDext_m[0] = extB_m / 10.0; // kgauss -> T
            // FDext_m[1] = extE_m;

            // // ---------------------------- Particle dumping ----------------------------------------- //
            // // Note: Don't dump when
            // // 1. after one turn
            // // in order to sychronize the dump step for multi-bunch and single
            // // bunch for compare with each other during post-process phase.
            // double const E = itsBunch->get_meanEnergy();

            // if(!(Options::psDumpLocalFrame)) {
            //     // --------------------------- Dump in global frame ----------------------------------- //
	    //     // Fixme: ROGERS: BUG - THIS IO ROUTINE CAUSES FACTOR 100 SLOW DOWN IN PROCESSING TIME!!!

            //     itsBunch->R /= Vector_t(1000.0); // mm --> m

            //     itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0, E);
            //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m, E, referencePr, referenceR, referenceTheta);

            //     itsBunch->R *= Vector_t(1000.0); // m --> mm

            //     *gmsg << "* Phase space dump " << lastDumpedStep_m
            //           << " (global frame) at integration step " << step_m + 1
            //           << " T = " << t << " [ns]" << endl;

            // } else {
            //     // --------------------------- Dump in local frame ----------------------------------- //

            //     double phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
	    //     Quaternion_t quaternionToXAxis;
            //     Vector_t const xaxis = Vector_t(1.0, 0.0, 0.0);

            //     getQuaternionTwoVectors(meanP, xaxis, quaternionToXAxis);

            //     globalToLocal(itsBunch->R, quaternionToXAxis, meanR);
            //     globalToLocal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //globalToLocal(itsBunch->R, phi, meanR);
            //     //globalToLocal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

            //     // dump in local frame
            //     itsBunch->R /= Vector_t(1000.0); // mm --> m

            //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m, E, referencePr, referenceR, referenceTheta);
            //     itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0, E);

            //     itsBunch->R *= Vector_t(1000.0); // m --> mm

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
            //     localToGlobal(itsBunch->R, quaternionToXAxis, meanR);
            //     localToGlobal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //localToGlobal(itsBunch->R, phi, meanR);
            //     //localToGlobal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

            //     *gmsg << "* Phase space dump " << lastDumpedStep_m << " (local frame) at integration step "
            //           << step_m + 1 << " T = " << itsBunch->getT() * 1e9 << " [ns], phi = " << phi/pi*180.0 <<" [deg]" <<endl;
            // }
            // // ------------- Done with dumping ------------------------------------------------------- //
            IpplTimings::stopTimer(DumpTimer_m);
        }
        if(!(step_m + 1 % 1000))
            *gmsg << "Step " << step_m + 1 << endl;
    } // end for: the integration is DONE after maxSteps_m steps!

    // some post-integration works
    *gmsg << endl;
    *gmsg << "* *---------------------------- DONE TRACKING PARTICLES --------------------------------* " << endl;

    // calculate tunes after tracking.

    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 0) {
            double FinalMomentum2 = pow(itsBunch->P[ii](0), 2.0) + pow(itsBunch->P[ii](1), 2.0) + pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsg << "* Total phase space dump number(includes the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
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
          *gmsg << "*" << endl;
	  *gmsg << "* Finished during turn " << turnnumber_m << " (" << turnnumber_m - 1 << " turns completed)" << endl;
          *gmsg << "* Cave: Turn number is not correct for restart mode"<< endl;
    }

    Ippl::Comm->barrier();

    if(myNode_m == 0) outfTrackOrbit_m.close();

    if(initialTotalNum_m == 1)
        closeFiles();

    *gmsg << *itsBunch << endl;

}

#ifdef GENERICTRACKER
/**
 *
 *
 */
void ParallelCyclotronTracker::Tracker_Generic() {
    // Generic Tracker that has three modes defined by timeIntegrator_m:
    // 0 --- RK-4 (default)
    // 1 --- LF-2
    // (2 --- MTS ... not yet implemented in generic tracker)
    // numBunch_m determines the number of bunches in multibunch mode (MBM, 1 for OFF)
    // Total number of particles determines single particle mode (SPM, 1 particle) or
    // Static Equilibrium Orbit Mode (SEO, 2 particles)

    //if (timeIntegrator_m == 1) BorisPusher pusher; // Create a BorisPusher only for LF-2 method
    BorisPusher pusher; // TEMP: OPAL will not compile without a pusher. -DW

    // Read in some control parameters
    const int SinglePartDumpFreq = Options::sptDumpFreq;
    const int resetBinFreq = Options::rebinFreq;
    const int scSolveFreq = Options::scSolveFreq;
    const int boundpDestroyFreq = Options::boundpDestroyFreq;
    const bool doDumpAfterEachTurn = Options::psDumpEachTurn;

    // Read in some more parameters
    const int stepsPerTurn = itsBunch->getStepsPerTurn();
    const double harm = getHarmonicNumber();
    const double dt = itsBunch->getdT() * 1.0e9 * harm; // time step size (s --> ns)
    double t  = itsBunch->getT() * 1.0e9;               // current time   (s --> ns)

    // Define 3 special azimuthal angles where dump particle's six parameters
    // at each turn into 3 ASCII files. only important in single particle tracking
    const double azimuth_angle0 = 0.0;
    const double azimuth_angle1 = 22.5 * PIOVER180;
    const double azimuth_angle2 = 45.0 * PIOVER180;

    const double initialReferenceTheta = referenceTheta * PIOVER180;
    double oldReferenceTheta = initialReferenceTheta; // init here, reset each step
    const double deltaTheta = pi / (stepsPerTurn);    // half of the average angle per step

    // Set some flags
    bool flagNoDeletion = false;
    bool flagTransition = false; // Flag for transition from single-bunch to multi-bunches mode
    bool outOfBound = false;     // Flag whether a particle is out of field

    int stepsNextCheck = step_m + stepsPerTurn; // Steps to next check for transition
    //int stepToLastInj = itsBunch->getSteptoLastInj(); // TODO: Do we need this? -DW

    // Record how many bunches have already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();
    BinCount_m = BunchCount_m;

    initTrackOrbitFile();

    // Get data from h5 file for restart run and reset current step
    // to last step of previous simulation
    if(OpalData::getInstance()->inRestartRun()) {

        restartStep0_m = itsBunch->getLocalTrackStep();
        step_m = restartStep0_m;

        if (numBunch_m > 1)
            itsBunch->resetPartBinID2(eta_m);

        *gmsg << "* Restart at integration step " << restartStep0_m << endl;
    }

    // TODO: Comment about what scan option does -DW
    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {

        lastDumpedStep_m = 0;
        itsBunch->setT(0.0);
        t = 0.0;
    }

    // For single particle mode open output files
    if(initialTotalNum_m == 1)
        openFiles(OpalData::getInstance()->getInputBasename());

    initDistInGlobalFrame();

    // If initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
    // Where is the IF here? -DW
    vector<double> Ttime, Tdeltr, Tdeltz;
    vector<int> TturnNumber;
    turnnumber_m = 1;

    // --- Output to user --- //
    *gmsg << "* Beginning of this run is at t = " << t << " [ns]" << endl;
    *gmsg << "* The time step is set to dt = " << dt << " [ns]" << endl;

    if(numBunch_m > 1) {

        *gmsg << "* MBM: Time interval between neighbour bunches is set to "
              << stepsPerTurn * dt << "[ns]" << endl;
        *gmsg << "* MBM: The particles energy bin reset frequency is set to "
              << resetBinFreq << endl;
    }

    *gmsg << "* Single particle trajectory dump frequency is set to " << SinglePartDumpFreq << endl;
    *gmsg << "* The frequency to solve space charge fields is set to " << scSolveFreq << endl;
    *gmsg << "* The repartition frequency is set to " << Options::repartFreq << endl;

    if(initialTotalNum_m == 1) {
        *gmsg << endl;
        *gmsg << "* ------------------------------ SINGLE PARTICLE MODE---------------------------------- *" << endl
              << "* Instruction: When the total particle number is equal to 1, single particle mode is    *" << endl
              << "* triggered automatically. The initial distribution file must be specified which should *" << endl
              << "* contain only one line for the single particle                                         *" << endl
              << "* ---------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE -------------- *" << endl;

        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

    } else if(initialTotalNum_m == 2) {
        *gmsg << endl;
        *gmsg << "* ------------------------- STATIC EQUILIBRIUM ORBIT MODE ----------------------------- *" << endl
              << "* Instruction: When the total particle number is equal to 2, SEO mode is triggered      *" << endl
              << "* automatically. This mode does NOT include any RF cavities. The initial distribution   *" << endl
              << "* file must be specified. In the file the first line is for reference particle and the  *" << endl
              << "* second line is for off-center particle. The tune is calculated by FFT routines based  *" << endl
              << "* on these two particles.                                                               *" << endl
              << "* ---------------- NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------ *" << endl;

        if(Ippl::getNodes() != 1)
            throw OpalException("Error in ParallelCyclotronTracker::execute", "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
    }

    // Apply the plugin elements: probe, collimator, stripper, septum once before first step
    applyPluginElements(dt);

    /********************************
     *     Main integration loop    *
     ********************************/
    *gmsg << endl;
    *gmsg << "* --------------------------------- Start tracking ------------------------------------ *" << endl;

    for(; step_m < maxSteps_m; step_m++) {

        bool dumpEachTurn = false;

        // Start Tracking for number of particles > 2 (i.e. not single and not SEO mode)
        if(initialTotalNum_m > 2) {

            // single particle dumping
            if(step_m % SinglePartDumpFreq == 0)
                singleParticleDump();

            //Ippl::Comm->barrier(); TEMP moved into singleParticleDump() -DW

            // If we are using the 2nd order Leap Frog method, push half a step
            /* TEMP move this right before the second half step for testing -DW
            if (timeIntegrator_m == 1) {

                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }
            */

            // Find out if we need to do bunch injection
            if (numBunch_m > 1) {

                if ((BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)) {

                    if (step_m == stepsNextCheck) {
                        // If all of the following conditions are met, this code will be executed
                        // to check the distance between two neighboring bunches:
                        // 1. Only one bunch exists (BunchCount_m == 1)
                        // 2. We are in multi-bunch mode, AUTO sub-mode (multiBunchMode_m == 2)
                        // 3. It has been a full revolution since the last check (stepsNextCheck)

                        *gmsg << "* MBM: Checking for automatically injecting new bunch ..." << endl;

                        itsBunch->R *= Vector_t(0.001); // mm --> m
                        itsBunch->calcBeamParameters_cycl();
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm

                        RThisTurn_m = sqrt(pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0));

                        Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

                        double XYrms =  sqrt(pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0));

                        // If the distance between two neighboring bunches is less than 5 times of its 2D rms size
                        // start multi-bunch simulation, fill current phase space to initialR and initialP arrays

                        if ((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                            // since next turn, start multi-bunches
                            saveOneBunch();
                            flagTransition = true;
                            *gmsg << "* MBM: Saving beam distribution at turn " << turnnumber_m << endl;
                            *gmsg << "* MBM: After one revolution, Multi-Bunch Mode will be invoked" << endl;

                        }

                        stepsNextCheck += stepsPerTurn;

                        *gmsg << "* MBM: RLastTurn = " << RLastTurn_m << " [mm]" << endl;
                        *gmsg << "* MBM: RThisTurn = " << RThisTurn_m << " [mm]" << endl;
                        *gmsg << "* MBM: XYrms = " << XYrms    << " [mm]" << endl;

                        RLastTurn_m = RThisTurn_m;
                    }

                }
                else if ((BunchCount_m < numBunch_m) && (step_m == stepsNextCheck)) {

                    // If all of the following conditions are met, this code will be executed
                    // to read new bunch from hdf5 format file:
                    // 1. We are in multi-bunch mode (numBunch_m > 1)
                    // 2. It has been a full revolution since the last check
                    // 3. Number of existing bunches is less than the desired number of bunches
                    // 4. FORCE mode, or AUTO mode with flagTransition = true
                    // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
                    *gmsg << "* MBM: Step " << step_m << ", injecting a new bunch... ... ..." << endl;

                    BunchCount_m++;

                    // read initial distribution from h5 file
                    if (multiBunchMode_m == 1) {

                        readOneBunch(BunchCount_m - 1);

                        if (timeIntegrator_m == 0) itsBunch->resetPartBinID2(eta_m);

                    } else if (multiBunchMode_m == 2) {

                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        if (timeIntegrator_m == 0) itsBunch->resetPartBinID2(eta_m);
                    }

                    itsBunch->setNumBunch(BunchCount_m);

                    stepsNextCheck += stepsPerTurn;

                    // update  after injection
                    itsBunch->boundp();

                    Ippl::Comm->barrier();

                    *gmsg << "* MBM: Bunch " << BunchCount_m
                          << " injected, total particle number = "
                          << itsBunch->getTotalNum() << endl;

                } else if (BunchCount_m == numBunch_m) {
                    // After this, numBunch_m is wrong but not needed anymore...
                    numBunch_m--;
                }
            }

            Vector_t const meanR = calcMeanR();
            oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

            // Calculate SC field before each time step and keep constant during integration.
            // Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
            if (itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {

                if (step_m % scSolveFreq == 0) {

                    IpplTimings::startTimer(TransformTimer_m);

                    // Firstly reset E and B to zero before fill new space charge field data for each track step
                    itsBunch->Bf = Vector_t(0.0);
                    itsBunch->Ef = Vector_t(0.0);

                    PreviousMeanP = calcMeanP();

                    // --- Multibunche mode --- //
                    if ((itsBunch->weHaveBins()) && BunchCount_m > 1) {

                        // Since calcMeanP takes into account all particles of all bins (TODO: Check this! -DW)
                        // Using the quaternion method with PreviousMeanP and yaxis should give the correct result

			Quaternion_t quaternionToYAxis;

                        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

                        globalToLocal(itsBunch->R, quaternionToYAxis, meanR);

                        itsBunch->R *= Vector_t(0.001); // mm --> m

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        // Calcualte gamma for each energy bin
                        itsBunch->calcGammas_cycl();

                        repartition();

                        // Calculate space charge field for each energy bin
                        // The relativistic effects will be slightly wrong with this method
                        // Also, I don't understand how this affects any interaction between
                        // the bunches? -DW
                        for(int b = 0; b < itsBunch->getLastemittedBin(); b++) {

                            itsBunch->setBinCharge(b, itsBunch->getChargePerParticle());
                            itsBunch->computeSelfFields_cycl(b, 0.001 * meanR, quaternionToYAxis);
                        }

                        itsBunch->Q = itsBunch->getChargePerParticle();

                        IpplTimings::startTimer(TransformTimer_m);

                        // Scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // Transform coordinates back to global
                        localToGlobal(itsBunch->R, quaternionToYAxis, meanR);

                        // Transform self field back to global frame (rotate only)
                        localToGlobal(itsBunch->Ef, quaternionToYAxis);
                        localToGlobal(itsBunch->Bf, quaternionToYAxis);

                    } else {
                        // --- Single bunch mode --- //
                        double temp_meangamma = sqrt(1.0 + dot(PreviousMeanP, PreviousMeanP));

	                Quaternion_t quaternionToYAxis;

                        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

                        globalToLocal(itsBunch->R, quaternionToYAxis, meanR);

                        itsBunch->R *= Vector_t(0.001); // mm --> m

                        if((step_m + 1) % boundpDestroyFreq == 0)
                            itsBunch->boundp_destroy();
                        else
                            itsBunch->boundp();

                        IpplTimings::stopTimer(TransformTimer_m);

                        repartition();

                        itsBunch->computeSelfFields_cycl(temp_meangamma, 0.001 * meanR, quaternionToYAxis);

                        IpplTimings::startTimer(TransformTimer_m);

                        //scale coordinates back
                        itsBunch->R *= Vector_t(1000.0); // m --> mm

                        // Transform coordinates back to global
                        localToGlobal(itsBunch->R, quaternionToYAxis, meanR);

                        // Transform self field back to global frame (rotate only)
                        localToGlobal(itsBunch->Ef, quaternionToYAxis);
                        localToGlobal(itsBunch->Bf, quaternionToYAxis);
                    }

                    IpplTimings::stopTimer(TransformTimer_m);

                } else {
		    // If we are not solving for the space charge fields at this time step
                    // we will apply the fields from the previous step and have to rotate them
                    // accordingly. For this we find the quaternion between the previous mean momentum (PreviousMeanP)
                    // and the current mean momentum (meanP) and rotate the fields with this quaternion.

                    Vector_t const meanP = calcMeanP();

                    Quaternion_t quaternionToNewMeanP;

                    getQuaternionTwoVectors(PreviousMeanP, meanP, quaternionToNewMeanP);

                    // Reset PreviousMeanP. Cave: This HAS to be after the quaternion is calculated!
                    PreviousMeanP = calcMeanP();

                    // Rotate the fields
                    globalToLocal(itsBunch->Ef, quaternionToNewMeanP);
                    globalToLocal(itsBunch->Bf, quaternionToNewMeanP);
                }
            }
            // If field solver is not available or number of particles is less than 1000
            // we don't need to do anything space charge related.

            // IpplTimings::startTimer(IntegrationTimer_m);

            // *** This was moved here b/c collision should be tested before the **********************
            // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
            // Apply the plugin elements: probe, collimator, stripper, septum
            applyPluginElements(dt);

	    // check if we loose particles at the boundary
	    bgf_main_collision_test();
            // ****************************************************************************************

            // ------------------------ Track all particles for one step --------------------------- //
            // --- Cave: For now, the LF-2 method push es all particles at one in push() so it   --- //
            // --- has to be done before entering the loop over all particles and again after    --- //
            // ------------------------------------------------------------------------------------- //

            if (timeIntegrator_m == 1) { // LF-2 method

                // Push all particles for first LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            for(size_t i = 0; i < (itsBunch->getLocalNum()); i++) {

                flagNoDeletion = true;

                if (timeIntegrator_m == 1) { // LF-2 method

                    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
                    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);

                    outOfBound = getFieldsAtPoint(t, i, externalE, externalB);
                    if (outOfBound) flagNoDeletion = false;

                    // Do LF2 momentum kick
                    pusher.kick(itsBunch->R[i],
                                itsBunch->P[i],
                                externalE,
                                externalB,
                                dt * 1.0e-9,
                                itsBunch->M[i] * 1.0e9,
                                itsBunch->Q[i] / q_e);

                }

                else if (timeIntegrator_m == 0) { // RK-4 method

                    for(int j = 0; j < 3; j++) {

                        variable_m[j] = itsBunch->R[i](j);    // [x,y,z] (mm)
                        variable_m[j+3] = itsBunch->P[i](j);  // [px,py,pz] (beta*gamma)
                        rold_m[j] = variable_m[j];            // used for gap crossing checking
                        pold_m[j] = variable_m[j+3];          // used for gap crossing
                    }

                    // Integrate for one RK4 step in the lab Cartesian frame (absolute value).
                    flagNoDeletion = rk4(variable_m, t, dt, i);

                    for(int j = 0; j < 3; j++) {

                        itsBunch->R[i](j) = variable_m[j];    // [x,y,z] (mm)
                        itsBunch->P[i](j) = variable_m[j+3];  // [px,py,pz] (beta*gamma)
                    }

                    // --- Cave: LF-2 checks for cavities during the "push" command --- //
                    // --- Maybe some future TODO to consolidate this? -DW          --- //

                    // If gap crossing happens, do momenta kicking
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
                                itsBunch->R[i](j) = variable_m[j] ;    // [x,y,z]    (mm)
                                itsBunch->P[i](j) = variable_m[j+3] ;  // [px,py,pz] (beta*gamma)
                            }

                            //momentum kick
                            RFkick(rfcav, t, dt1, i);

                            // retrack particle  from cavity gap point for the left time to finish the entire timestep
                            for(int j = 0; j < 3; j++) {
                                variable_m[j] = itsBunch->R[i](j);    // [x,y,z]    (mm)
                                variable_m[j+3] = itsBunch->P[i](j);  // [px,py,pz] (beta*gamma)
                            }

                            if(dt / dt2 < 1.0e9) rk4(variable_m, t, dt2, i);

                            for(int j = 0; j < 3; j++) {
                                itsBunch->R[i](j) = variable_m[j] ;    // [x,y,z]    (mm)
                                itsBunch->P[i](j) = variable_m[j+3] ;  // [px,py,pz] (beta*gamma)
                            }
                        } // end if: gap-crossing monentum kicking at certain cavity
                    } // end for: finish checking for all cavities
                } // end if: which mode are we in LF2 or RK4
            } // end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode

            if (timeIntegrator_m == 1) { // LF-2 method

                // Push all particles for second LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            // Destroy particles if they are marked as Bin = -1 in the plugin elements
            // or out of global aperture
            bool flagNeedUpdate = deleteParticle();

            // If particles were deleted, recalculate bingamma and reset BinID for remaining particles
            if(itsBunch->weHaveBins() && flagNeedUpdate)
                itsBunch->resetPartBinID2(eta_m);

	    // Recalculate bingamma and reset the BinID for each particles according to its current gamma
	    if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % resetBinFreq == 0)
	        itsBunch->resetPartBinID2(eta_m);

            // Some status output for user after each turn
            if((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {
                turnnumber_m++;
                dumpEachTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch->getTotalNum() << endl;
            }

            // IpplTimings::stopTimer(IntegrationTimer_m);
            Ippl::Comm->barrier();

        } // end if: more than 2 particles -> bunch or multi-bunch tracking

        else if (initialTotalNum_m == 2) {

            // 2 particles: Trigger SEO mode
            // (Switch off cavity and calculate betatron osciliation tuning)
            double r_tuning[2], z_tuning[2] ;

            if (timeIntegrator_m == 1) { // LF-2 method

                // Push for first LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            for(size_t i = 0; i < (itsBunch->getLocalNum()); i++) {

                for(int j = 0; j < 3; j++) variable_m[j]   = itsBunch->R[i](j); //[x,y,z]    (mm)
                for(int j = 0; j < 3; j++) variable_m[j+3] = itsBunch->P[i](j); //[px,py,pz] (beta*gamma)

                if((step_m % SinglePartDumpFreq == 0)) {

                    outfTrackOrbit_m << "ID" << (itsBunch->ID[i]);
                    outfTrackOrbit_m << " " << variable_m[0]
                                     << " " << variable_m[3]
                                     << " " << variable_m[1]
                                     << " " << variable_m[4]
                                     << " " << variable_m[2]
                                     << " " << variable_m[5] << endl;
                }

                double OldTheta = 0.0;

                OldTheta = calculateAngle(variable_m[0], variable_m[1]);
                r_tuning[i] = variable_m[0] * cos(OldTheta) + variable_m[1] * sin(OldTheta);
                z_tuning[i] = variable_m[2];

                //IpplTimings::startTimer(IntegrationTimer_m);

                if (timeIntegrator_m == 1) { // LF-2 method

                    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
                    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);

                    getFieldsAtPoint(t, i, externalE, externalB);

                    IpplTimings::startTimer(IntegrationTimer_m);

                    // Do LF2 momentum kick
                    pusher.kick(itsBunch->R[i],
                                itsBunch->P[i],
                                externalE,
                                externalB,
                                dt * 1.0e-9,
                                itsBunch->M[i] * 1.0e9,
                                itsBunch->Q[i] / q_e);

                    IpplTimings::stopTimer(IntegrationTimer_m);

                }

                else if (timeIntegrator_m == 0) { // RK-4 method

                    // Integrate for one step in the lab Cartesian frame (absulate value).
                    rk4(variable_m, t, dt, i);

                    for(int j = 0; j < 3; j++) itsBunch->R[i](j) = variable_m[j] ;   //[x,y,z]    (mm)
                    for(int j = 0; j < 3; j++) itsBunch->P[i](j) = variable_m[j+3] ; //[px,py,pz] (beta*gamma)
                }

                //IpplTimings::stopTimer(IntegrationTimer_m);

		if( (i == 0) && (step_m > 10) && ((step_m%stepsPerTurn) == 0)) turnnumber_m++;
            } // end for: finish one step tracking for all particles for initialTotalNum_m == 2 mode

            if (timeIntegrator_m == 1) { // LF-2 method

                // Push for first LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            // store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
            if(step_m % SinglePartDumpFreq == 0) {

                Ttime.push_back(t * 1.0e-9);
                Tdeltz.push_back(z_tuning[1]);
                Tdeltr.push_back(r_tuning[1] - r_tuning[0]);
                TturnNumber.push_back(turnnumber_m);
            }
        } // end if: SEO mode for 2 particles

        else if(initialTotalNum_m == 1) {

            // 1 particle: Trigger single particle mode

            //IpplTimings::startTimer(IntegrationTimer_m);

            flagNoDeletion = true;

            // *** This was moved here b/c collision should be tested before the **********************
            // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
            // apply the plugin elements: probe, collimator, stripper, septum
            applyPluginElements(dt);

	    // check if we loose particles at the boundary
	    bgf_main_collision_test();
            // ****************************************************************************************

            // Track for one step
            if (timeIntegrator_m == 1) { // LF-2 method

                // Push for first LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            for( unsigned int i = 0; i < itsBunch->getLocalNum(); i++) {

                if((step_m % SinglePartDumpFreq == 0)) {

                    outfTrackOrbit_m << "ID" <<itsBunch->ID[i];
                    outfTrackOrbit_m << " " << itsBunch->R[i](0) << " " << itsBunch->P[i](0) << " " << itsBunch->R[i](1)
                                     << " " << itsBunch->P[i](1) << " " << itsBunch->R[i](2) << " " << itsBunch->P[i](2) << endl;
                }

                for(int j = 0; j < 3; j++) {

                    variable_m[j] = itsBunch->R[i](j);    // [x,y,z] (mm)
                    variable_m[j+3] = itsBunch->P[i](j);  // [px,py,pz] (beta*gamma)
                    rold_m[j] = variable_m[j];            // used for gap crossing checking
                    pold_m[j] = variable_m[j+3];          // used for gap crossing checking
                }

                double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]); // [-pi, pi]

                if ((step_m > 10) && ((step_m + 1) % stepsPerTurn) == 0) {

                    ++turnnumber_m;

                    dumpEachTurn = true;

                    *gmsg << "* SPT: Finished turn " << turnnumber_m - 1 << endl;

                    outfThetaEachTurn_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
                    outfThetaEachTurn_m << " " << sqrt(variable_m[0] * variable_m[0] + variable_m[1] * variable_m[1])
                                        << " " << variable_m[3] * cos(temp_meanTheta) + variable_m[4] * sin(temp_meanTheta)
                                        << " " << temp_meanTheta / PIOVER180
                                        << " " << -variable_m[3] * sin(temp_meanTheta) + variable_m[4] * cos(temp_meanTheta)
                                        << " " << variable_m[2]
                                        << " " << variable_m[5] << endl;
                }

                if ((oldReferenceTheta < azimuth_angle0 - deltaTheta) && (temp_meanTheta >= azimuth_angle0 - deltaTheta)) {
                    outfTheta0_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
                    outfTheta0_m << " " << sqrt(variable_m[0] * variable_m[0] + variable_m[1] * variable_m[1])
                                 << " " << variable_m[3] * cos(temp_meanTheta) + variable_m[4] * sin(temp_meanTheta)
                                 << " " << temp_meanTheta / PIOVER180
                                 << " " << -variable_m[3] * sin(temp_meanTheta) + variable_m[4] * cos(temp_meanTheta)
                                 << " " << variable_m[2]
                                 << " " << variable_m[5] << endl;
                }

                if ((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                    outfTheta1_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
                    outfTheta1_m << " " << sqrt(variable_m[0] * variable_m[0] + variable_m[1] * variable_m[1])
                                 << " " << variable_m[3] * cos(temp_meanTheta) + variable_m[4] * sin(temp_meanTheta)
                                 << " " << temp_meanTheta / PIOVER180
                                 << " " << -variable_m[3] * sin(temp_meanTheta) + variable_m[4] * cos(temp_meanTheta)
                                 << " " << variable_m[2]
                                 << " " << variable_m[5] << endl;
                }

                if ((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                    outfTheta2_m << "#Turn number = " << turnnumber_m << ", Time = " << t << " [ns]" << endl;
                    outfTheta2_m << " " << sqrt(variable_m[0] * variable_m[0] + variable_m[1] * variable_m[1])
                                 << " " << variable_m[3] * cos(temp_meanTheta) + variable_m[4] * sin(temp_meanTheta)
                                 << " " << temp_meanTheta / PIOVER180
                                 << " " << -variable_m[3] * sin(temp_meanTheta) + variable_m[4] * cos(temp_meanTheta)
                                 << " " << variable_m[2]
                                 << " " << variable_m[5] << endl;
                }

                oldReferenceTheta = temp_meanTheta;

                if (timeIntegrator_m == 1) { // LF-2 method

                    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
                    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);

                    outOfBound = getFieldsAtPoint(t, i, externalE, externalB);
                    if (outOfBound) flagNoDeletion = false;

                    // Do LF2 momentum kick
                    IpplTimings::startTimer(IntegrationTimer_m);

                    pusher.kick(itsBunch->R[i], itsBunch->P[i], externalE, externalB, dt * 1.0e-9, itsBunch->M[i] * 1.0e9, itsBunch->Q[i] / q_e);

                    IpplTimings::stopTimer(IntegrationTimer_m);

                    if(!flagNoDeletion) {

                        *gmsg << "* SPT: The particle was lost at step " << step_m << "!" << endl;
                        throw OpalException("ParallelCyclotronTracker", "The particle is out of the region of interest.");
                    }
                }

                else if (timeIntegrator_m == 0) { // RK-4 method

                    // integrate for one step in the lab Cartesian frame (absolute value).
                    flagNoDeletion = rk4(variable_m, t, dt, i);

                    if(!flagNoDeletion) {

                        *gmsg << "* SPT: The particle was lost at step " << step_m << "!" << endl;
                        throw OpalException("ParallelCyclotronTracker", "The particle is out of the region of interest.");
                    }

                    for(int j = 0; j < 3; j++) itsBunch->R[i](j) = variable_m[j] ; //[x,y,z]  units: [mm]
                    for(int j = 0; j < 3; j++) itsBunch->P[i](j) = variable_m[j+3] ; //[px,py,pz]  units: [] beta*gama

                    // --- Cave: LF-2 checks for cavities during the "push" command --- //
                    // --- Maybe some future TODO to consolidate this? -DW          --- //

                    // If gap crossing happens, do momenta kicking
                    for (beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++) {

                        bool tag_crossing = false;
                        double DistOld = 0.0; //mm
                        RFCavity * rfcav;

                        if (((*sindex)->first) == "CAVITY") {

                            // here check gap cross in the list, if do , set tag_crossing to TRUE
                            for(int j = 0; j < 3; j++)
                                rnew_m[j] = variable_m[j];

                            rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
                            tag_crossing = checkGapCross(rold_m, rnew_m, rfcav, DistOld);
                        }

                        if (tag_crossing) {

                            double oldMomentum2  = dot(pold_m, pold_m);
                            double oldBetgam = sqrt(oldMomentum2);
                            double oldGamma = sqrt(1.0 + oldMomentum2);
                            double oldBeta = oldBetgam / oldGamma;
                            double dt1 = DistOld / (c * oldBeta * 1.0e-6); // ns
                            double dt2 = dt - dt1;

                            // retrack particle from the old postion to cavity gap point
                            // restore the old coordinates and momenta
                            for(int j = 0; j < 3; j++) {

                                variable_m[j]   = rold_m[j]; // [x,y,z]    (mm)
                                variable_m[j+3] = pold_m[j]; // [px,py,pz] (beta*gamma)
                            }

                            if(dt / dt1 < 1.0e9) rk4(variable_m, t, dt1, i);

                            for(int j = 0; j < 3; j++) {

                                itsBunch->R[i](j) = variable_m[j] ;   // [x,y,z]    (mm)
                                itsBunch->P[i](j) = variable_m[j+3] ; // [px,py,pz] (beta*gamma)
                            }

                            // Momentum kick
                            RFkick(rfcav, t, dt1, i);

                            // Retrack particle from cavity gap point for the left time to finish the entire timestep
                            for (int j = 0; j < 3; j++) {

                                variable_m[j]   = itsBunch->R[i](j); // [x,y,z]    (mm)
                                variable_m[j+3] = itsBunch->P[i](j); // [px,py,pz] (beta*gamma)
                            }

                            if (dt / dt2 < 1.0e9) rk4(variable_m, t, dt2, i);

                            for (int j = 0; j < 3; j++) {

                                itsBunch->R[i](j) = variable_m[j] ;   // [x,y,z]    (mm)
                                itsBunch->P[i](j) = variable_m[j+3] ; // [px,py,pz] (beta*gamma)
                            }
                        } // end if: gap-crossing monentum kicking at certain cavity
                    } // end for: finish checking for all cavities
                } // end if: LF2 or RK4 method
            } // end for: track one particle one step

            if (timeIntegrator_m == 1) { // LF-2 method

                // Push for second LF2 half step
                itsBunch->R *= Vector_t(0.001);  // mm --> m
                push(0.5 * dt * 1.0e-9);         // ns --> s
                itsBunch->R *= Vector_t(1000.0); // m  --> mm

            }

            // Destroy particles if they are marked as Bin = -1 in the plugin elements
            // or out of the global aperture
            deleteParticle();

            //IpplTimings::stopTimer(IntegrationTimer_m);

        } // end if: finish one step tracking either for
          // initialTotalNum_m > 2 || initialTotalNum_m == 2 || initialTotalNum_m == 1 mode

        // Update bunch and some parameters and output some info

        // Reference parameters
        double tempP2 = dot(itsBunch->P[0], itsBunch->P[0]);
        double tempGamma = sqrt(1.0 + tempP2);
        double tempBeta = sqrt(tempP2) / tempGamma;

        PathLength_m += c_mmtns * dt / 1000.0 * tempBeta; // unit: m

        t += dt;
        itsBunch->setT((t) * 1.0e-9);
        itsBunch->setLPath(PathLength_m);
        itsBunch->setLocalTrackStep((step_m + 1)); // TEMP moved this here from inside if statement below -DW
        // Here is global frame, don't do itsBunch->boundp();

        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2)
            || (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {

            // Write Phase Space and Statistics Data to h5 and dat files
            // TODO: It seems strange that ps and stat dump are not separate...-DW
            bunchDumpPhaseSpaceStatData();
        }

        if(!(step_m + 1 % 1000))
            *gmsg << "Step " << step_m + 1 << endl;

    } // end for: the integration is DONE after maxSteps_m steps!

    // Some post-integration stuff
    *gmsg << endl;
    *gmsg << "* ---------------------------- DONE TRACKING PARTICLES -------------------------------- * " << endl;

    // Calculate tunes after tracking.

    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 0) {
            double FinalMomentum2 = pow(itsBunch->P[ii](0), 2.0) + pow(itsBunch->P[ii](1), 2.0) + pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsg << "* Total phase space dump number(includes the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step (-restart " << lastDumpedStep_m << ")" << endl;
        }
    }

    Ippl::Comm->barrier();

    if(initialTotalNum_m == 2) {
	*gmsg << endl;
        *gmsg << "* **************** The result for tune calulation (NO space charge) ******************* *" << endl
              << "* Number of tracked turns: " << TturnNumber.back() << endl;
        double nur, nuz;
        getTunes(Ttime, Tdeltr, Tdeltz, TturnNumber.back(), nur, nuz);

    } else {
        // not for multibunch
        if(!(itsBunch->weHaveBins()))
          *gmsg << "*" << endl;
	  *gmsg << "* Finished during turn " << turnnumber_m << " (" << turnnumber_m - 1 << " turns completed)" << endl;
          *gmsg << "* Cave: Turn number is not correct for restart mode"<< endl;
    }

    Ippl::Comm->barrier();

    if(myNode_m == 0) outfTrackOrbit_m.close();

    if(initialTotalNum_m == 1) closeFiles();

    *gmsg << *itsBunch << endl;
}

bool ParallelCyclotronTracker::getFieldsAtPoint(const double &t, const size_t &Pindex, Vector_t &Efield, Vector_t &Bfield) {

    beamline_list::iterator sindex = FieldDimensions.begin();

    bool outOfBound = (((*sindex)->second).second)->apply(Pindex, t, Efield, Bfield);

    Bfield *= 0.10;  // kGauss --> T
    Efield *= 1.0e6; // kV/mm  --> V/m

    // For runs without space charge effects, override this step to save time
    if(itsBunch->hasFieldSolver()) {

        // Don't do for reference particle
        if(itsBunch->ID[Pindex] != 0) {

            // add external Field and self space charge field
            Efield += itsBunch->Ef[Pindex];
            Bfield += itsBunch->Bf[Pindex];
        }
    }

    return outOfBound;
}
#endif //GENERICTRACKER

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
bool ParallelCyclotronTracker::derivate(double *y, const double &t, double *yp, const size_t &Pindex) {

    Vector_t externalE, externalB, tempR;

    externalB = Vector_t(0.0, 0.0, 0.0);
    externalE = Vector_t(0.0, 0.0, 0.0);

    for(int i = 0; i < 3; i++)
        tempR(i) = y[i];

    itsBunch->R[Pindex] = tempR;

#ifdef GENERICTRACKER
    bool outOfBound = getFieldsAtPoint(t, Pindex, externalE, externalB);

#else
    beamline_list::iterator sindex = FieldDimensions.begin();

    bool outOfBound = (((*sindex)->second).second)->apply(Pindex, t, externalE, externalB);

    for(int i = 0; i < 3; i++) externalB(i) = externalB(i) * 0.10;  // kGauss --> T
    for(int i = 0; i < 3; i++) externalE(i) = externalE(i) * 1.0e6; // kV/mm  --> V/m

    // for working modes without space charge effects, override this step to save time
    if(itsBunch->hasFieldSolver()) {
        if(itsBunch->ID[Pindex] != 0) {
            // add external Field and self space charge field
            externalE += itsBunch->Ef[Pindex];
            externalB += itsBunch->Bf[Pindex];
        }
    }

#endif // GENERICTRACKER

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

    return outOfBound;
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
bool ParallelCyclotronTracker::rk4(double x[], const double &t, const double &tau, const size_t &Pindex) {
    // Forth order Runge-Kutta integrator
    // arguments:
    //   x          Current value of dependent variable
    //   t          Independent variable (usually time)
    //   tau        Step size (usually time step)
    //   Pindex     index of particel, not used yet

    IpplTimings::startTimer(IntegrationTimer_m);

    bool outOfBound = false;
    double  deriv1[PSdim];
    double  deriv2[PSdim];
    double  deriv3[PSdim];
    double  deriv4[PSdim];
    double  xtemp[PSdim];

    // Evaluate f1 = f(x,t).

    outOfBound = derivate(x, t, deriv1 , Pindex);
    if (outOfBound) {
        IpplTimings::stopTimer(IntegrationTimer_m);
        return false;
    }

    // Evaluate f2 = f( x+tau*f1/2, t+tau/2 ).
    const double half_tau = 0.5 * tau;
    const double t_half = t + half_tau;

    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + half_tau * deriv1[i];

    outOfBound = derivate(xtemp, t_half, deriv2 , Pindex);
    if (outOfBound) {
        IpplTimings::stopTimer(IntegrationTimer_m);
        return false;
    }

    // Evaluate f3 = f( x+tau*f2/2, t+tau/2 ).
    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + half_tau * deriv2[i];

    outOfBound = derivate(xtemp, t_half, deriv3 , Pindex);
    if (outOfBound) {
        IpplTimings::stopTimer(IntegrationTimer_m);
        return false;
    }

    // Evaluate f4 = f( x+tau*f3, t+tau ).
    double t_full = t + tau;
    for(int i = 0; i < PSdim; i++)
        xtemp[i] = x[i] + tau * deriv3[i];

    outOfBound = derivate(xtemp, t_full, deriv4 , Pindex);
    if (outOfBound) {
        IpplTimings::stopTimer(IntegrationTimer_m);
        return false;
    }

    // Return x(t+tau) computed from fourth-order R-K.
    for(int i = 0; i < PSdim; i++)
        x[i] += tau / 6.*(deriv1[i] + deriv4[i] + 2.*(deriv2[i] + deriv3[i]));

    IpplTimings::stopTimer(IntegrationTimer_m);
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
    rfcavity -> getMomentaKick(nomalRadius, tempP, t, dt, itsBunch->ID[Pindex], itsBunch->getM(), itsBunch->getQ()); // t : ns

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

    *gmsg << endl;
    *gmsg << "* ************************************* nuR ******************************************* *" << endl;
    *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

    int nhis_lomb = 10;
    int  stat = 0;
    // book tune class
    tune = new TUNE_class();
    stat = tune->LombAnalysis(t, r, nhis_lomb, T / lastTurn);
    if(stat != 0)
        *gmsg << "* TUNE: Lomb analysis failed" << endl;
    *gmsg << "* ************************************************************************************* *" << endl;

    delete tune;
    tune = NULL;
    // FIXME: FixMe: need to come from the inputfile
    nhis_lomb = 100;

    if(zsum != 0.0) {

        *gmsg << endl;
        *gmsg << "* ************************************* nuZ ******************************************* *" << endl;
        *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

        // book tune class
        tune = new TUNE_class();
        stat = tune->LombAnalysis(t, z, nhis_lomb, T / lastTurn);
        if(stat != 0)
            *gmsg << "* TUNE: Lomb analysis failed" << endl;
        *gmsg << "* ************************************************************************************* *" << endl;

        delete tune;
        tune = NULL;
    }
    return true;
}

void ParallelCyclotronTracker::saveOneBunch() {

    *gmsg << endl;
    *gmsg << "* ---------------- Clone Beam---------------- *" << endl;

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

    std::string fn_appendix = "-onebunch";
    itsDataSink->storeOneBunch(*itsBunch, fn_appendix);

}

/**
 *
 * @param BinID
 * @param step
 */
bool ParallelCyclotronTracker::readOneBunch(const size_t BinID) {

    *gmsg << endl;
    *gmsg << "* ---------------- Copy Beam---------------- *" << endl;

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

    static bool restartflag = true;

    if(restartflag) {
        *gmsg << "----------------Read beam from hdf5 file----------------" << endl;
        const size_t LocalNum = itsBunch->getLocalNum();

        std::string fn_appendix = "-onebunch";
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

double ParallelCyclotronTracker::getHarmonicNumber() const {
    if (opalRing_m != NULL)
        return opalRing_m->getHarmonicNumber();
    Cyclotron* elcycl = dynamic_cast<Cyclotron*>(((*FieldDimensions.begin())->second).second);
    if (elcycl != NULL)
        return elcycl->getCyclHarm();
    throw OpalException("ParallelCyclotronTracker::Tracker_MTS()",
             std::string("The first item in the FieldDimensions list does not ")
            +std::string("seem to be an OpalRing or a Cyclotron element"));
}

void ParallelCyclotronTracker::Tracker_MTS() {
	IpplTimings::startTimer(IpplTimings::getTimer("MTS"));
	IpplTimings::startTimer(IpplTimings::getTimer("MTS-Various"));

    const double harm = getHarmonicNumber();
    const double dt = itsBunch->getdT() * harm;
    if(numBunch_m > 1) {
        *gmsg << "Time interval between neighbour bunches is set to " << itsBunch->getStepsPerTurn() * dt * 1.0e9 << "[ns]" << endl;
    }

    initTrackOrbitFile();

    // get data from h5 file for restart run
    if(OpalData::getInstance()->inRestartRun()) {
        restartStep0_m = itsBunch->getLocalTrackStep();
        step_m = restartStep0_m;
	    if (numBunch_m > 1) itsBunch->resetPartBinID2(eta_m);
        *gmsg << "* Restart at integration step " << restartStep0_m << endl;
    }
    if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        lastDumpedStep_m = 0;
        itsBunch->setT(0.0);
    }
    *gmsg << "* Beginning of this run is at t= " << itsBunch->getT() * 1e9 << " [ns]" << endl;
    *gmsg << "* The time step is set to dt= " << dt * 1e9 << " [ns]" << endl;

    // for single Particle Mode, output at zero degree.
    if(initialTotalNum_m == 1) {
        openFiles(OpalData::getInstance()->getInputBasename());
    }

    initDistInGlobalFrame();
    itsBunch->R *= Vector_t(0.001); // In MTS method, we use [R] = m for calculation
    RLastTurn_m *= 0.001;
    RThisTurn_m *= 0.001;
    double const initialReferenceTheta = referenceTheta / 180.0 * pi;
    *gmsg << "* Single particle trajectory dump frequency is set to " << Options::sptDumpFreq << endl;
    *gmsg << "* Repartition frequency is set to " << Options::repartFreq << endl;
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

    int const numSubsteps = std::max(Options::mtsSubsteps, 1);
    *gmsg << "MTS: Number of substeps per step is " << numSubsteps << endl;
    double const dt_inner = dt / double(numSubsteps);
    *gmsg << "MTS: The inner time step is therefore " << dt_inner << endl;
    int SteptoLastInj = itsBunch->getSteptoLastInj();
    double oldReferenceTheta = initialReferenceTheta;
    turnnumber_m = 1;
    bool flagTransition = false; // flag to determine when to transit from single-bunch to multi-bunches mode
    int stepsNextCheck = step_m + itsBunch->getStepsPerTurn(); // step point determining the next time point of check for transition
    const double deltaTheta = pi / itsBunch->getStepsPerTurn();
    *gmsg << "* ---------------------------- Start tracking ----------------------------" << endl;
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Various"));
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-SpaceCharge"));
    if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
        evaluateSpaceChargeField();
    }
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-SpaceCharge"));
    for(; step_m < maxSteps_m; step_m++) {
    	IpplTimings::startTimer(IpplTimings::getTimer("MTS-Dump"));
        bool dumpEachTurn = false;
        if(step_m % Options::sptDumpFreq == 0) {
            itsBunch->R *= Vector_t(1000.0);
            singleParticleDump();
            itsBunch->R *= Vector_t(0.001);
        }
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Dump"));

        // First half kick from space charge force
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Kick"));
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            kick(0.5 * dt);
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Kick"));

        // Substeps for external field integration
        for(int n = 0; n < numSubsteps; ++n) {
            borisExternalFields(dt_inner);
        }

		IpplTimings::startTimer(IpplTimings::getTimer("MTS-Various"));
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
                        *gmsg << "*** Save beam distribution at turn #" << turnnumber_m << " ***" << endl;
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
                        itsBunch->resetPartBinID2(eta_m);
                    } else if(multiBunchMode_m == 2) {
                        if(OpalData::getInstance()->inRestartRun())
                            readOneBunchFromFile(BunchCount_m - 1);
                        else
                            readOneBunch(BunchCount_m - 1);

                        itsBunch->resetPartBinID2(eta_m);
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
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Various"));

	IpplTimings::startTimer(IpplTimings::getTimer("MTS-SpaceCharge"));

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
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-SpaceCharge"));

        // Second half kick from space charge force
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Kick"));
        if(itsBunch->hasFieldSolver() && initialTotalNum_m >= 1000) {
            kick(0.5 * dt);
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Kick"));

        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Various"));
        if((itsBunch->weHaveBins()) && BunchCount_m > 1) {
            if(step_m % Options::rebinFreq == 0) {
                itsBunch->resetPartBinID2(eta_m);
            }
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Various"));

        // dump some data after one push in single particle tracking
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Dump"));
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
                ++turnnumber_m;
                *gmsg << "Turn " << turnnumber_m << endl;
                dumpEachTurn = true;
                outfThetaEachTurn_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
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
                outfTheta0_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta0_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle1 - deltaTheta) && (temp_meanTheta >= azimuth_angle1 - deltaTheta)) {
                outfTheta1_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta1_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }

            if((oldReferenceTheta < azimuth_angle2 - deltaTheta) && (temp_meanTheta >= azimuth_angle2 - deltaTheta)) {
                outfTheta2_m << "#Turn number = " << turnnumber_m << ", Time = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
                outfTheta2_m << " " << sqrt(variable_m[0]*variable_m[0] + variable_m[1]*variable_m[1])
                             << " " << variable_m[3]*cos(temp_meanTheta) + variable_m[4]*sin(temp_meanTheta)
                             << " " << temp_meanTheta / pi * 180
                             << " " << -variable_m[3]*sin(temp_meanTheta) + variable_m[4]*cos(temp_meanTheta)
                             << " " << variable_m[2]
                             << " " << variable_m[5] << endl;
            }
            oldReferenceTheta = temp_meanTheta;
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Dump"));

        // check whether one turn over for multi-bunch tracking.
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Various"));
        if(Options::psDumpEachTurn && initialTotalNum_m > 2) {
            Vector_t const meanR = calcMeanR();

            // in global Cartesian frame, calculate the location in global frame of bunch
            oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

            // both for single bunch and multi-bunch
            // avoid dump at the first step
            // dumpEachTurn has not been changed in first push
            if((step_m > 10) && ((step_m + 1) % itsBunch->getStepsPerTurn()) == 0) {
                ++turnnumber_m;
                dumpEachTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch->getTotalNum() << endl;
            }
        }
        // reset Bin ID for each particle
        if((itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m % Options::rebinFreq == 0)
            itsBunch->resetPartBinID2(eta_m);
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Various"));

        // dump phase space distribution of bunch
        IpplTimings::startTimer(IpplTimings::getTimer("MTS-Dump"));

        if((((step_m + 1) % Options::psDumpFreq == 0) && initialTotalNum_m != 2) ||
           (Options::psDumpEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {

            IpplTimings::startTimer(DumpTimer_m);

            itsBunch->setSteptoLastInj(SteptoLastInj);
            itsBunch->setLocalTrackStep((step_m + 1));

            // Go from m --> mm to be compatible with bunchDumpPhaseSpaceStatData()
            // which is called also from Tracker_RK4 and Tracker_LF which are using mm
            itsBunch->R *= Vector_t(1000.0);

            // Write Phase Space and Statistics Data to h5 and dat files
	    bunchDumpPhaseSpaceStatData();

	    // Go back from mm --> m
            itsBunch->R *= Vector_t(0.001);

            // //--------------------- calculate mean coordinates of bunch -------------------------------//
            // //------------ and calculate the external field at the mass of bunch -----------------------//

            // Vector_t const meanR = calcMeanR();
            // Vector_t const meanP = calcMeanP();

            // // define longitudinal direction of the bunch
            // //double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

	    // // *gmsg << "meanR = ( " << meanR(0) * 1000.0 << " " << meanR(1) * 1000.0 << " " << meanR(2) * 1000.0 << " ) [mm] " << endl;
	    // extE_m = Vector_t(0.0, 0.0, 0.0);
            // extB_m = Vector_t(0.0, 0.0, 0.0);
            // beamline_list::iterator DumpSindex = FieldDimensions.begin();
            // (((*DumpSindex)->second).second)->apply(meanR * 1000.0, Vector_t(0.0), itsBunch->getT() * 1e9, extE_m, extB_m);
            // FDext_m[0] = extB_m * 0.1; // kgauss -> T
            // FDext_m[1] = extE_m;

            // // --------------------------- Particle dumping ------------------------------------------ //
            // // Note: Don't dump when
            // // 1. after one turn
            // // in order to sychronize the dump step for multi-bunch and single bunch for compare
            // // with each other during post-process phase.
            // double const E = itsBunch->get_meanEnergy();

            // if(!(Options::psDumpLocalFrame)) {
	    //     // --------------------------- Dump in global frame ----------------------------------- //

            //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m,E, referencePr, referenceR, referenceTheta);
            //     // itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);  // TODO: why no stat data in global frame?

            //     *gmsg << "* Phase space dump " << lastDumpedStep_m << " (global frame) at integration step "
            //           << step_m + 1 << " T = " << itsBunch->getT() * 1e9 << " [ns]" << endl;

            // } else {
	    //     // --------------------------- Dump in global frame ----------------------------------- //

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
	    //     Quaternion_t quaternionToXAxis;
            //     Vector_t const xaxis = Vector_t(1.0, 0.0, 0.0);

            //     getQuaternionTwoVectors(meanP, xaxis, quaternionToXAxis);

            //     globalToLocal(itsBunch->R, quaternionToXAxis, meanR);
            //     globalToLocal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //globalToLocal(itsBunch->R, phi, meanR);
            //     //globalToLocal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

            //     lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m, E, referencePr, referenceR, referenceTheta);
            //     itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0, E);

            //     // TEMP for testing -DW ********************************************
            //     // NEW:
            //     localToGlobal(itsBunch->R, quaternionToXAxis, meanR);
            //     localToGlobal(itsBunch->P, quaternionToXAxis, meanP);

            //     // OLD:
            //     //localToGlobal(itsBunch->R, phi, meanR);
            //     //localToGlobal(itsBunch->P, phi, meanP);
            //     // TEMP END ********************************************************

            //     *gmsg << "* Phase space dump " << lastDumpedStep_m
            //           << " (local frame) at integration step " << step_m + 1
            //           << " T = " << itsBunch->getT() * 1e9 << " [ns]" << endl;
            IpplTimings::stopTimer(DumpTimer_m);
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Dump"));
    }
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-Dump"));

    for(size_t ii = 0; ii < itsBunch->getLocalNum(); ++ii) {
        if(itsBunch->ID[ii] == 0) {
            double FinalMomentum2  = pow(itsBunch->P[ii](0), 2.0) +
                                     pow(itsBunch->P[ii](1), 2.0) +
                                     pow(itsBunch->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch->getM() * 1.0e-6;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsg << "* Total phase space dump number including the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step ( -restart " << lastDumpedStep_m << " )" << endl;
        }
    }
    Ippl::Comm->barrier();
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Dump"));
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-Various"));
    if(myNode_m == 0) outfTrackOrbit_m.close();
    if(initialTotalNum_m == 1) closeFiles();
    *gmsg << *itsBunch << endl;


    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Various"));
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS"));
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
    particleVectors -= translationToGlobal;
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
    particleVectors += translationToGlobal;
}


inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors, Quaternion_t const quaternion, Vector_t const meanR) {

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateWithQuaternion(particleVectors, quaternion);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors, Quaternion_t const quaternion, Vector_t const meanR) {

    // Reverse the quaternion by multiplying the axis components (x,y,z) with -1
    Quaternion_t reverseQuaternion = quaternion * -1.0;
    reverseQuaternion(0) *= -1.0;

    // Rotation back to original P_mean direction
    rotateWithQuaternion(particleVectors, reverseQuaternion);

    // Translation from local to global
    particleVectors += meanR;
}

inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors, double const phi, double const psi, Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, phi);

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, psi);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors, double const phi, double const psi, Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate back around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, -psi);

     // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, -phi);

    // Translation from local to global
    particleVectors += meanR;
}

inline void ParallelCyclotronTracker::rotateWithQuaternion(ParticleAttrib<Vector_t> & particleVectors, Quaternion_t const quaternion) {

    Vector_t const quaternionVectorComponent = Vector_t(quaternion(1), quaternion(2), quaternion(3));
    double const quaternionScalarComponent = quaternion(0);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        particleVectors[i] = 2.0f * dot(quaternionVectorComponent, particleVectors[i]) * quaternionVectorComponent +
                             (quaternionScalarComponent * quaternionScalarComponent -
                              dot(quaternionVectorComponent, quaternionVectorComponent)) * particleVectors[i] + 2.0f *
                              quaternionScalarComponent * cross(quaternionVectorComponent, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::normalizeQuaternion(Quaternion_t & quaternion){

    double tolerance = 1.0e-10;
    double length2 = dot(quaternion, quaternion);

    if (fabs(length2) > tolerance && fabs(length2 - 1.0f) > tolerance) {

        double length = sqrt(length2);
        quaternion /= length;
    }
}

inline void ParallelCyclotronTracker::normalizeVector(Vector_t & vector) {

    double tolerance = 1.0e-10;
    double length2 = dot(vector, vector);

    if (fabs(length2) > tolerance && fabs(length2 - 1.0f) > tolerance) {

        double length = sqrt(length2);
        vector /= length;
    }
}

inline void ParallelCyclotronTracker::rotateAroundZ(ParticleAttrib<Vector_t> & particleVectors, double const phi) {
    // Clockwise rotation of particles 'particleVectors' by 'phi' around Z axis

    Tenzor<double, 3> const rotation( cos(phi), sin(phi), 0,
                                     -sin(phi), cos(phi), 0,
                                             0,        0, 1);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::rotateAroundX(ParticleAttrib<Vector_t> & particleVectors, double const psi) {
    // Clockwise rotation of particles 'particleVectors' by 'psi' around X axis

    Tenzor<double, 3> const rotation(1,  0,          0,
				     0,  cos(psi), sin(psi),
				     0, -sin(psi), cos(psi));

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::getQuaternionTwoVectors(Vector_t u, Vector_t v, Quaternion_t & quaternion) {
    // four vector (w,x,y,z) of the quaternion of P_mean with the positive x-axis

    normalizeVector(u);
    normalizeVector(v);

    double k_cos_theta = dot(u, v);
    double k = sqrt(dot(u, u) * dot(v, v));
    double tolerance1 = 1.0e-5;
    double tolerance2 = 1.0e-8;
    Vector_t resultVectorComponent;

    if (fabs(k_cos_theta / k + 1.0) < tolerance1) {
        // u and v are almost exactly antiparallel so we need to do
        // 180 degree rotation around any vector orthogonal to u

        resultVectorComponent = cross(u, xaxis);

        // If by chance u is parallel to xaxis, use zaxis instead
	if (dot(resultVectorComponent, resultVectorComponent) < tolerance2) {

	    resultVectorComponent = cross(u, zaxis);
        }

        double halfAngle = 0.5 * pi;
        double sinHalfAngle = sin(halfAngle);

        resultVectorComponent *= sinHalfAngle;

        k = 0.0;
        k_cos_theta = cos(halfAngle);

    } else {

        resultVectorComponent = cross(u, v);

    }

    quaternion(0) = k_cos_theta + k;
    quaternion(1) = resultVectorComponent(0);
    quaternion(2) = resultVectorComponent(1);
    quaternion(3) = resultVectorComponent(2);

    normalizeQuaternion(quaternion);
}

void ParallelCyclotronTracker::push(double h) {
    IpplTimings::startTimer(IntegrationTimer_m);

    std::list<CavityCrossData> cavCrossDatas;
    for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); ++sindex) {
        if(((*sindex)->first) == "CAVITY") {
            RFCavity * cav = static_cast<RFCavity *>(((*sindex)->second).second);
            CavityCrossData ccd = {cav, cav->getSinAzimuth(), cav->getCosAzimuth(), cav->getPerpenDistance() * 0.001};
            cavCrossDatas.push_back(ccd);
        }
    }
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        Vector_t const oldR = itsBunch->R[i];
        double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
        double const c_gamma = c / gamma;
        Vector_t const v = itsBunch->P[i] * c_gamma;
        itsBunch->R[i] += h * v;
        for(std::list<CavityCrossData>::const_iterator ccd = cavCrossDatas.begin(); ccd != cavCrossDatas.end(); ++ccd) {
            double const distNew = (itsBunch->R[i][0] * ccd->sinAzimuth - itsBunch->R[i][1] * ccd->cosAzimuth) - ccd->perpenDistance;
            bool tagCrossing = false;
            double distOld;
            if(distNew <= 0.0) {
                distOld = (oldR[0] * ccd->sinAzimuth - oldR[1] * ccd->cosAzimuth) - ccd->perpenDistance;
                if(distOld > 0.0) tagCrossing = true;
            }
            if(tagCrossing) {
                double const dt1 = distOld / sqrt(dot(v, v));
                double const dt2 = h - dt1;

                // Retrack particle from the old postion to cavity gap point
                itsBunch->R[i] = oldR + dt1 * v;

                // Momentum kick
                itsBunch->R[i] *= 1000.0; // RFkick uses [itsBunch->R] == mm
                RFkick(ccd->cavity, itsBunch->getT() * 1.0e9, dt1 * 1.0e9, i);
                itsBunch->R[i] *= 0.001;

                itsBunch->R[i] += dt2 * itsBunch->P[i] * c_gamma;
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
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-PushAndRFKick"));
    push(0.5 * h);
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-PushAndRFKick"));

    // Evaluate external fields
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-EvalExternal"));
    IpplTimings::startTimer(IntegrationTimer_m);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->Ef[i] = Vector_t(0.0, 0.0, 0.0);
        itsBunch->Bf[i] = Vector_t(0.0, 0.0, 0.0);
        beamline_list::iterator sindex = FieldDimensions.begin();
        itsBunch->R[i] *= 1000.0;
        (((*sindex)->second).second)->apply(i, itsBunch->getT() * 1e9, itsBunch->Ef[i], itsBunch->Bf[i]);
        itsBunch->R[i] *= 0.001; // mm -> m
        itsBunch->Bf[i] *= 0.1;  // kgauss -> T
    }
    IpplTimings::stopTimer(IntegrationTimer_m);
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-EvalExternal"));

    // Kick particles for full step
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-Kick"));
    kick(h);
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-Kick"));

    // push particles for second half step
    IpplTimings::startTimer(IpplTimings::getTimer("MTS-PushAndRFKick"));
    push(0.5 * h);
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-PushAndRFKick"));

	IpplTimings::startTimer(IpplTimings::getTimer("MTS-PluginElements"));
    // apply the plugin elements: probe, collimator, stripper, septum
    itsBunch->R *= Vector_t(1000.0); // applyPluginElements expects [R] = mm
    applyPluginElements(h * 1e9); // expects [dt] = ns
    // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global apeture
    bool const flagNeedUpdate = deleteParticle();

    itsBunch->R *= Vector_t(0.001);
    if(itsBunch->weHaveBins() && flagNeedUpdate) itsBunch->resetPartBinID2(eta_m);
    IpplTimings::stopTimer(IpplTimings::getTimer("MTS-PluginElements"));
}

void ParallelCyclotronTracker::applyPluginElements(const double dt) {

  for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++) {
    if(((*sindex)->first) == "SEPTUM")    {
      (static_cast<Septum *>(((*sindex)->second).second))->checkSeptum(*itsBunch);
    }

    if(((*sindex)->first) == "PROBE")    {
      (static_cast<Probe *>(((*sindex)->second).second))->checkProbe(*itsBunch, turnnumber_m, itsBunch->getT() * 1e9, dt);
    }

    if(((*sindex)->first) == "STRIPPER")    {
      bool flag_stripper = (static_cast<Stripper *>(((*sindex)->second).second))
        -> checkStripper(*itsBunch, turnnumber_m, itsBunch->getT() * 1e9, dt);
      if(flag_stripper) {
        itsBunch->boundp();
        *gmsg << "* Total number of particles after stripping = " << itsBunch->getTotalNum() << endl;
      }
    }

    if(((*sindex)->first) == "CCOLLIMATOR") {
      Collimator * collim;
      collim = static_cast<Collimator *>(((*sindex)->second).second);
      collim->checkCollimator(*itsBunch, turnnumber_m, itsBunch->getT() * 1e9, dt);
    }
  }
}

bool ParallelCyclotronTracker::deleteParticle(){
  // Update immediately if any particles are lost during this step

  bool flagNeedUpdate = (min(itsBunch->Bin) < 0);
  reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

  size_t lostParticleNum = 0;

  if(flagNeedUpdate) {
      Vector_t const meanR = calcMeanR();
      Vector_t const meanP = calcMeanP();

      // Bunch (local) azimuth at meanR w.r.t. y-axis
      double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

      // Bunch (local) elevation at meanR w.r.t. xy plane
      double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

      // For statistics data, the bunch is transformed into a local coordinate system
      // with meanP in direction of y-axis -DW
      globalToLocal(itsBunch->R, phi, psi, meanR);
      globalToLocal(itsBunch->P, phi, psi, Vector_t(0.0)); // P should be rotated, but not shifted

      itsBunch->R *= Vector_t(0.001); // mm --> m

      for(unsigned int i = 0; i < itsBunch->getLocalNum(); i++) {
          if(itsBunch->Bin[i] < 0) {
              lostParticleNum++;
              itsBunch->destroy(1, i);
          }
      }

      // Now destroy particles and update pertinent parameters in local frame
      // Note that update() will be called within boundp() -DW
      itsBunch->boundp();
      //itsBunch->update();

      itsBunch->calcBeamParameters_cycl();

      itsBunch->R *= Vector_t(1000.0); // m --> mm

      localToGlobal(itsBunch->R, phi, psi, meanR);
      localToGlobal(itsBunch->P, phi, psi, Vector_t(0.0));

      reduce(lostParticleNum, lostParticleNum, OpAddAssign());
      INFOMSG("Step " << step_m << ", " << lostParticleNum << " particles lost on stripper, collimator, septum, or out of cyclotron aperture" << endl);
  }
  return flagNeedUpdate;
}

void ParallelCyclotronTracker::initTrackOrbitFile() {

    std::string f = OpalData::getInstance()->getInputBasename() + std::string("-trackOrbit.dat");

    outfTrackOrbit_m.setf(ios::scientific, ios::floatfield);
    outfTrackOrbit_m.precision(8);

    if(myNode_m == 0) {

        if(OpalData::getInstance()->inRestartRun()) {

            outfTrackOrbit_m.open(f.c_str(), ios::app);
            outfTrackOrbit_m << "# Restart at integration step " << itsBunch->getLocalTrackStep() << endl;
        } else {

            outfTrackOrbit_m.open(f.c_str());
            outfTrackOrbit_m << "# The six-dimensional phase space data in the global Cartesian coordinates" << endl;
            outfTrackOrbit_m << "# Part. ID    x [mm]       beta_x*gamma       y [mm]      beta_y*gamma        z [mm]      beta_z*gamma" << endl;
        }
    }
}

void ParallelCyclotronTracker::initDistInGlobalFrame() {
    if(!OpalData::getInstance()->inRestartRun()) {
        // Start a new run (no restart)

        double const initialReferenceTheta = referenceTheta * PIOVER180;
        PathLength_m = 0.0;

        // Force the initial phase space values of the particle with ID = 0 to zero,
        // to set it as a reference particle.
        if(initialTotalNum_m > 2) {
            for(size_t i = 0; i < initialLocalNum_m; ++i) {
                if(itsBunch->ID[i] == 0) {
                    itsBunch->R[i] = Vector_t(0.0);
                    itsBunch->P[i] = Vector_t(0.0);
                }
            }
        }

        // Changed the order of things here: Do all the initializations first and
        // Save initial phase space dump at the end (local or global) in order
        // To correctly consider the z momentum in case of axial injection -DW

        // Initialize global R
        itsBunch->R *= Vector_t(1000.0); // m --> mm

        Vector_t const initMeanR = Vector_t(referenceR * cosRefTheta_m, // [referenceR] = mm
                                            referenceR * sinRefTheta_m,
                                            referenceZ);                // [referenceZ] = mm

        localToGlobal(itsBunch->R, initialReferenceTheta, initMeanR);

        // Initialize global P (Cartesian, but input P_ref is in Pr, Ptheta, Pz,
        // so translation has to be done before the rotation this once)
        // Cave: In the local frame, the positive y-axis is the direction of movement -DW
        for(size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch->P[i](0) += referencePr;
            itsBunch->P[i](1) += referencePt;
            itsBunch->P[i](2) += referencePz;
	}

        // Out of the three coordinates of meanR (R, Theta, Z) only the angle
        // changes the momentum vector...
        localToGlobal(itsBunch->P, initialReferenceTheta);

        // Initialize the bin number of the first bunch to 0
        for(size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch->Bin[i] = 0;
        }

        // Backup initial distribution if multi bunch mode
        if((initialTotalNum_m > 2) && (numBunch_m > 1) && (multiBunchMode_m == 1)) {

            initialR_m = new Vector_t[initialLocalNum_m];
            initialP_m = new Vector_t[initialLocalNum_m];

            for(size_t i = 0; i < initialLocalNum_m; ++i) {

                initialR_m[i] = itsBunch->R[i];
                initialP_m[i] = itsBunch->P[i];
            }
        }
    } else {
        // Restart from the distribution in the h5 file

        PathLength_m = itsBunch->getLPath();
        itsBunch->R *= Vector_t(1000.0); // m --> mm

        /*
        // Omit restarting in local fram in future releases! -DW
        if((Options::psDumpLocalFrame)) {
            // Restart from the distribution in the local frame
            // Multi-bunch must not be done in the local frame!
            // If you restart from OPAL-T h5 file, you must dump in the local frame!
            // TODO: In the future, this should be decided by the mode of the
            // H5 file containing restart data, not from the mode of the current
            // run! -DW

            *gmsg<<"* Restart in the local frame" <<endl;

            double const initialReferenceTheta = referenceTheta / 180.0 * pi;
            PathLength_m = 0.0;

            // Initialize global R
            itsBunch->R *= Vector_t(1000.0); // m --> mm

            Vector_t const initMeanR = Vector_t(referenceR * cosRefTheta_m, // [referenceR] == mm
                                                referenceR * sinRefTheta_m,
                                                referenceZ);                // [referenceZ] == mm

            localToGlobal(itsBunch->R, initialReferenceTheta, initMeanR);

            // Initialize global P (Cartesian, but input P_ref is in Pr, Ptheta, Pz,
            // so translation has to be done before the rotation this once)
            // Cave: In the local frame, the positive y-axis is the direction of movement -DW
            for(size_t i = 0; i < initialLocalNum_m; ++i) {
                itsBunch->P[i](0) += referencePr;
                itsBunch->P[i](1) += referencePt;
                itsBunch->P[i](2) += referencePz;
            }

            localToGlobal(itsBunch->P, initialReferenceTheta);

            // Initialize the bin number of the first bunch to 0
            for(size_t i = 0; i < initialLocalNum_m; ++i) {
                itsBunch->Bin[i] = 0;
            }

        } else {
            // Restart from the distribution in the global frame

            *gmsg<<"* Restart in the global frame" <<endl;

            PathLength_m = itsBunch->getLPath();
            itsBunch->R *= Vector_t(1000.0); // m --> mm
        }
        */
        //END omission of local restart -DW
    }

    // ------------ Get some Values ---------------------------------------------------------- //
    Vector_t const meanR = calcMeanR();
    Vector_t const meanP = calcMeanP();

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

    // AUTO mode
    if(multiBunchMode_m == 2) {

        RLastTurn_m = sqrt(meanR[0] * meanR[0] + meanR[1] * meanR[1]);
        RThisTurn_m = RLastTurn_m;

        if(OpalData::getInstance()->inRestartRun()) {

            *gmsg << "Radial position at restart position = ";

        } else {

            *gmsg << "Initial radial position = ";
        }

        *gmsg << RThisTurn_m << " mm" << endl;
    }

    // Do boundp and repartition in the local frame at beginning of this run
    globalToLocal(itsBunch->R, phi, psi, meanR);
    globalToLocal(itsBunch->P, phi, psi, Vector_t(0.0)); // P should be rotated, but not shifted

    itsBunch->R *= Vector_t(0.001); // mm --> m

    itsBunch->boundp();

    checkNumPart(std::string("* Before repartition: "));
    repartition();
    checkNumPart(std::string("* After repartition: "));

    itsBunch->calcBeamParameters_cycl();

    itsBunch->R *= Vector_t(1000.0); // m --> mm

    localToGlobal(itsBunch->R, phi, psi, meanR);
    localToGlobal(itsBunch->P, phi, psi, Vector_t(0.0));

    // Save initial distribution if not a restart
    if(!OpalData::getInstance()->inRestartRun()) {

        step_m -= 1;

        bunchDumpPhaseSpaceStatData();

        step_m += 1;
    }

    // Print out the Bunch information at beginning of the run
    *gmsg << *itsBunch << endl;
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

                if(rmsg == 0)
                    ERRORMSG("Could not receive from client nodes in main." << endl);

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

            for(int i = 0; i < counter; i++) {

                smsg->put(itsBunch->ID[found[i]]);

                for(int j = 0; j < 3; j++) {

                    smsg->put(itsBunch->R[found[i]](j));
                    smsg->put(itsBunch->P[found[i]](j));
                }
            }

            if(!Ippl::Comm->send(smsg, 0, tag))
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
        }

        Ippl::Comm->barrier();

    } else {

        for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {

            if(itsBunch->ID[i] == 0 || itsBunch->ID[i] == 1) {

                outfTrackOrbit_m << "ID" << itsBunch->ID[i] << " ";
                outfTrackOrbit_m << itsBunch->R[i](0) << " " << itsBunch->P[i](0) << " ";
                outfTrackOrbit_m << itsBunch->R[i](1) << " " << itsBunch->P[i](1) << " ";
                outfTrackOrbit_m << itsBunch->R[i](2) << " " << itsBunch->P[i](2) << endl;
            }
        }
    }

    IpplTimings::stopTimer(DumpTimer_m);
}

void ParallelCyclotronTracker::bunchDumpPhaseSpaceStatData() {
    // --------------------------- Particle dumping ------------------------------------------ //
    // Note: Don't dump when
    // 1. after one turn
    // in order to sychronize the dump step for multi-bunch and single bunch for comparison
    // with each other during post-process phase.
    // --------------------------------------------------------------------------------------- //
    IpplTimings::startTimer(DumpTimer_m);
    // ------------ Get some Values ---------------------------------------------------------- //
    double const temp_t = itsBunch->getT() * 1e9; // s -> ns
    Vector_t const meanR = calcMeanR();
    Vector_t const meanP = calcMeanP();

    double const betagamma_temp = sqrt(dot(meanP, meanP));

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double const theta = atan2(meanR(1), meanR(0));

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

    // Re-calculate reference values in format of input values
    // Position:
    referenceR = sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1));
    referenceTheta = theta * 180.0 / pi;
    referenceZ = meanR(2);
    // Momentum in Theta-hat, R-hat, Z-hat coordinates at position meanR:
    referencePtot = betagamma_temp;
    referencePz = meanP(2);
    referencePr = meanP(0) * cos(theta) + meanP(1) * sin(theta);
    referencePt = sqrt(referencePtot * referencePtot - referencePz * referencePz - referencePr * referencePr);

    // -----------  Calculate the external fields at the center of the bunch ------------------//
    beamline_list::iterator DumpSindex = FieldDimensions.begin();

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    (((*DumpSindex)->second).second)->apply(meanR, Vector_t(0.0), temp_t, extE_m, extB_m);
    FDext_m[0] = extB_m / 10.0; // kgauss --> T
    FDext_m[1] = extE_m;

    // if(!(Options::psDumpLocalFrame)) {
    // --------------------- Always dump whole bunch in global frame ------------------------ //

    // For statistics data, the bunch is transformed into a local coordinate system
    // with meanP in direction of y-axis -DW
    globalToLocal(itsBunch->R, phi, psi, meanR);  // meanR mm --> m
    globalToLocal(itsBunch->P, phi, psi, Vector_t(0.0)); // P should be rotated, but not shifted

    itsBunch->R *= Vector_t(0.001); // mm --> m

    // Calculate the beam parameters in a local frame but with the magnitude of the momentum
    // vector unchanged
    //itsBunch->calcBeamParameters_cycl();
    itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0, 0.0, 0.0);

    localToGlobal(itsBunch->R, phi, psi, meanR / Vector_t(1000.0)); // Right now: R (m), meanR (mm)
    localToGlobal(itsBunch->P, phi, psi, Vector_t(0.0));

    double const E = itsBunch->get_meanEnergy();

    lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, // Global and in (m)
                                                         FDext_m, E,
                                                         referencePr,
                                                         referencePt,
                                                         referencePz,
                                                         referenceR,
                                                         referenceTheta,
                                                         referenceZ,
                                                         phi * 180.0 / pi,  // P_mean azimuth at ref. R/Th/Z
                                                         psi * 180.0 / pi); // P_mean elevation at ref. R/Th/Z

    itsBunch->R *= Vector_t(1000.0); // m --> mm

    // Tell user we are dumping data
    *gmsg << endl;
    *gmsg << "* Phase space dump " << lastDumpedStep_m
          << " at integration step " << step_m + 1 << endl;

	//}
    /* // No more dumping in local frame! -DW
    else {
        // --------------------------- Dump in local frame ---------------------------------- //

        globalToLocal(itsBunch->R, phi, psi, meanR);
        globalToLocal(itsBunch->P, phi, psi, meanP);

	itsBunch->R /= Vector_t(1000.0); // mm --> m

	lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch,
                                                             FDext_m, E,
                                                             referencePr,
                                                             referencePz,
                                                             referenceR,
                                                             referenceTheta,
                                                             referenceZ);

	itsDataSink->writeStatData(*itsBunch, FDext_m , 0.0, 0.0, 0.0, E);

	itsBunch->R *= Vector_t(1000.0); // m --> mm

        localToGlobal(itsBunch->R, phi, psi, meanR);
        localToGlobal(itsBunch->P, phi, psi, meanP);

        // Tell user in which mode we are dumping
	*gmsg << "* Phase space dump " << lastDumpedStep_m
              << " (local frame) at integration step " << step_m + 1 << endl;

    }
    */
    // END omission of local frame dumping

    // Print dump information on screen
    *gmsg << "* T = " << temp_t << " ns"
          << ", Live Particles: " << itsBunch->getTotalNum() << endl;
    *gmsg << "* E = " << E << " MeV"
          << ", beta * gamma = " << betagamma_temp << endl;
    *gmsg << "* Bunch position: R =  " << referenceR << " mm"
          << ", Theta = " << referenceTheta << " Deg"
          << ", Z = " << referenceZ << " mm" << endl;
    *gmsg << "* Local Azimuth = " << phi / pi * 180.0 << " Deg"
          << ", Local Elevation = " << psi / pi * 180.0 << " Deg" << endl;

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
