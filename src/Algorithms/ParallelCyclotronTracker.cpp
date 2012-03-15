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
// $Author: adelmann jang $
//
// ------------------------------------------------------------------------

/*
  - Have RK4 and LF the very same functions implemented?
  


*/
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
#include "H5Part.h"

class Beamline;
class PartData;
using Physics::c;
using Physics::m_p; // GeV
using Physics::PMASS;
using Physics::PCHARGE;
using Physics::pi;

const double c_mmtns = c * 1.0e-6; // m/s --> mm/ns

#define PSdim 6

extern Inform *gmsg;

typedef FVector<double,PSdim> Vector;
typedef FMatrix<double,PSdim,PSdim> Matrix;
typedef FTps<double,PSdim> Series;
typedef FVps<double,PSdim> Map, VSeries;
typedef FMatrix<FTps<double,PSdim>,PSdim,PSdim> MxSeries;

static Vector implicitIntStep(const Vector &zin, const VSeries &f, double s, double ds, int nx = 20);
static Vector implicitInt2(const Vector &zin, const VSeries &f, double s, double ds, int nx = 20, int cx = 4);
static Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx = 20, int cx = 4);
static Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx = 50);
static Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx = 50, int cx = 4);

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
	flagDoTune_m(false),
	myNode_m(Ippl::myNode()),
	initialLocalNum_m(0),
	initialTotalNum_m(0)
{
    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
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
	flagDoTune_m(false),
	myNode_m(Ippl::myNode()),
	initialLocalNum_m(bunch.getLocalNum()),
	initialTotalNum_m(bunch.getTotalNum())
{
    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    //  scaleFactor_m = itsBunch->getdT() * c;
    scaleFactor_m = 1;
    multiBunchMode_m = 0;

    // default value for proton
    ratioCh_M_m = 1.0;
    chtmc_m = PCHARGE/(PMASS*c);
    chtm_m = PCHARGE/PMASS;

    IntegrationTimer_m = IpplTimings::getTimer("Integration");
    TransformTimer_m   = IpplTimings::getTimer("Frametransform");
    DumpTimer_m        = IpplTimings::getTimer("Dump");
    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");
}

/** 
 * Destructor ParallelCyclotronTracker
 * 
 */
ParallelCyclotronTracker::~ParallelCyclotronTracker()
{
    for (list<Component*>::iterator compindex = myElements.begin(); compindex != myElements.end(); compindex++)
		{
			delete (*compindex);
		}
    for (beamline_list::iterator fdindex = FieldDimensions.begin(); fdindex != FieldDimensions.end(); fdindex++)
		{
			delete (*fdindex);
		}
    delete itsBeamline;
}

/** 
 * 
 * 
 * @param fn Base file name 
 */
void ParallelCyclotronTracker::openFiles(string SfileName) {

	string  SfileName2 = SfileName+string("-Angle0.dat");

	outfTheta0_m.precision(8);
	outfTheta0_m.setf(ios::scientific, ios::floatfield);

	outfTheta0_m.open(SfileName2.c_str());
	outfTheta0_m<<"#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
				  << endl;
			
	SfileName2 = SfileName+string("-Angle1.dat");
		
	outfTheta1_m.precision(8);
	outfTheta1_m.setf(ios::scientific, ios::floatfield);
		
	outfTheta1_m.open(SfileName2.c_str());
	outfTheta1_m<<"#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
				  << endl;

	SfileName2 = SfileName+string("-Angle2.dat");
	
	outfTheta2_m.precision(8);
	outfTheta2_m.setf(ios::scientific, ios::floatfield);

	outfTheta2_m.open(SfileName2.c_str());
	outfTheta2_m<<"#  r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"
				  << endl;

	// for single Particle Mode, output after each turn, to define matched initial phase ellipse.

	SfileName2 = SfileName+string("-afterEachTurn.dat");

	outfThetaEachTurn_m.precision(8);
	outfThetaEachTurn_m.setf(ios::scientific, ios::floatfield);
	
	outfThetaEachTurn_m.open(SfileName2.c_str());
	outfThetaEachTurn_m<<"# r [mm]          p_r[rad]       theta [mm]          p_theta[rad]        z [mm]          p_z[rad]"<<endl;

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
void ParallelCyclotronTracker::visitCyclotron(const Cyclotron &cycl)
{
    Inform msg("visitCyclotron ");

    myElements.push_back(dynamic_cast<Cyclotron*>(cycl.clone()));
    Component *elptr = *(--myElements.end());

    double ri = elptr->getRinit();
    msg << "RINIT= " << ri << " [mm]" << endl;
    referenceR = ri;

    double pri = elptr->getPRinit();
    //msg << "PRINIT= " << pri << " [CU]" << endl;
    referencePr = pri;

    double phii = elptr->getPHIinit();
    msg << "PHIINIT= " << phii << " [deg]" << endl;
    referenceTheta = phii;
    if ( referenceTheta <= -180.0 || referenceTheta > 180.0 ){
		throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron","PHIINIT is out of [-180, 180)!");
    }

    referencePz = 0.0;
    referencePtot =  itsReference.getGamma() * itsReference.getBeta();
    referencePt = sqrt(referencePtot*referencePtot - referencePr*referencePr);
    if ( referencePtot < 0.0) referencePt *= -1.0;

    sinRefTheta = sin(phii/180.0*pi);
    cosRefTheta = cos(phii/180.0*pi);

    msg<<"Initial gamma = "<<itsReference.getGamma()<<endl;

    msg<<"Initial beta = "<<itsReference.getBeta()<<endl;

    msg<<"Total reference momentum   = "<<referencePtot*1000.0<<" [MCU]"<<endl;

    msg<<"Reference azimuthal momentum  = "<<referencePt*1000.0<<" [MCU]" <<endl;

    msg<<"Reference radial momentum     = "<<referencePr*1000.0<<" [MCU]"<<endl;

    double sym = elptr->getSymmetry();
    msg << "Field symmetry= " << sym << " " << endl;

    double rff = elptr->getRfFrequ();
    msg << "Rf frequency= " << rff << " [MHz]" << endl;

    string fmfn = elptr->getFieldMapFN();
    msg << "Field map file name= " << fmfn << " " << endl;

    string type = elptr->getType();
    msg << "Type of cyclotron= " << type << " " << endl;

    double h = elptr->getCyclHarm();
    msg << "Harmonic number h= " << h << " " << endl;

    // FIXME: this is total wrong ..... needs e ficed:
	// a) doc below has to go to the cyclotron cmd
	// b) type can be pasased to the init routine

	// fieldflag is a flag defining which field function is used to readin the magnetic field.its value is 
	// defined by the option "TYPE" of the element  "CYCLOTRON"
	// fieldflag = 1, readin PSI format measured field file (default)
	// fieldflag = 2, readin carbon cyclotron field file created by Jianjun Yang, TYPE=CARBONCYCL
	// fieldflag = 3, readin ANSYS format file for CYCIAE-100 created by Jianjun Yang, TYPE=CYCIAE
	// fieldflag = 4, readin AVFEQ format file for Riken cyclotrons

	// you can add the field flag of yourself here
	int  fieldflag;
	if ( type == string("CARBONCYCL")){
		fieldflag = 2;
	} else if ( type == string("CYCIAE")){
		fieldflag = 3;
	} else if ( type == string("AVFEQ")){
		fieldflag = 4;
	} else if ( type == string("FFAG")){
		fieldflag = 5;
	} else 
		fieldflag = 1;
  
	// read field map on the  middle plane of cyclotron.
	// currently scalefactor is set to 1.0
	elptr->initialise( itsBunch, fieldflag, 1.0 );

	double BcParameter[8];
	for (int i =0; i<8;i++) BcParameter[i]=0.0;
	string ElementType = "CYCLOTRON";
	BcParameter[0] = elptr->getRmin();
	BcParameter[1] = elptr->getRmax();

	// store inner radius and outer radius of cyclotron field map in the list
	buildupFieldList(BcParameter,ElementType, elptr);

}

/** 
 * Not implemented and most probable never used
 * 
 */
void ParallelCyclotronTracker::visitBeamBeam(const BeamBeam &)
{
    *gmsg << "In BeamBeam tracker is missing "<< endl;
}

/** 
 * 
 * 
 * @param coll 
 */
void ParallelCyclotronTracker::visitCollimator(const Collimator &coll)
{
    *gmsg << "In Collimator; L= " << coll.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Collimator*>(coll.clone()));
    //   applyDrift(flip_s * coll.getElementLength());
}

/** 
 * 
 * 
 * @param corr 
 */
void ParallelCyclotronTracker::visitCorrector(const Corrector &corr)
{
    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Corrector*>(corr.clone()));
}

/** 
 * 
 * 
 * @param diag 
 */
void ParallelCyclotronTracker::visitDiagnostic(const Diagnostic &diag)
{
    *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Diagnostic*>(diag.clone()));
}

/** 
 * 
 * 
 * @param drift 
 */
void ParallelCyclotronTracker::visitDrift(const Drift &drift)
{
    *gmsg << "In drift L= " << drift.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Drift*>(drift.clone()));
}

/** 
 * 
 * 
 * @param lamb 
 */
void ParallelCyclotronTracker::visitLambertson(const Lambertson &lamb)
{
    *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Lambertson*>(lamb.clone()));
}

/** 
 * 
 * 
 * @param marker 
 */
void ParallelCyclotronTracker::visitMarker(const Marker &marker)
{
    //   *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Marker*>(marker.clone()));
    // Do nothing.
}

/** 
 * 
 * 
 * @param corr 
 */
void ParallelCyclotronTracker::visitMonitor(const Monitor &corr)
{
    //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Monitor*>(corr.clone()));
    //   applyDrift(flip_s * corr.getElementLength());
}

/** 
 * 
 * 
 * @param mult 
 */
void ParallelCyclotronTracker::visitMultipole(const Multipole &mult)
{
    *gmsg << "In Multipole; L= " << mult.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Multipole*>(mult.clone()));
}

/** 
 * 
 * 
 * @param prob 
 */
void ParallelCyclotronTracker::visitProbe(const Probe &prob)
{
	Inform msg("visitProbe ");

	myElements.push_back(dynamic_cast<Probe*>(prob.clone()));
	Component *elptr = *(--myElements.end());

	double xstart = (dynamic_cast<Probe*>(elptr))->getXstart();
	msg << "XStart= " << xstart << " [mm]" << endl;

	double xend = (dynamic_cast<Probe*>(elptr)) ->getXend();
	msg << "XEnd= " << xend << " [mm]" << endl;

	double ystart =  (dynamic_cast<Probe*>(elptr))->getYstart();
	msg << "YStart= " << ystart<< " [mm]" << endl;

	double yend = (dynamic_cast<Probe*>(elptr)) ->getYend();
	msg << "YEnd= " << yend<< " [mm]" << endl;

	double width = (dynamic_cast<Probe*>(elptr)) ->getWidth();
	msg << "Width= " << width<< " [mm]" << endl;

   
	// initialise, do nothing
	elptr->initialise(itsBunch,1.0);

	double BcParameter[8];
	for (int i =0; i<8;i++) 
		BcParameter[i]=0.0;
	string ElementType = "PROBE";
	BcParameter[0] =xstart ;
	BcParameter[1] =xend;
	BcParameter[2] =ystart ;
	BcParameter[3] =yend;
	BcParameter[4] =width ;
   
	// store probe parameters in the list
	buildupFieldList(BcParameter,ElementType, elptr);
}

/** 
 * 
 * 
 * @param bend 
 */
void ParallelCyclotronTracker::visitRBend(const RBend &bend)
{
    *gmsg << "In RBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RBend*>(bend.clone()));
}

/** 
 * 
 * 
 * @param as 
 */
void ParallelCyclotronTracker::visitRFCavity(const RFCavity &as)
{

    Inform msg("visitRFCavity ");
    msg << "---------------------------------------" << endl;
    myElements.push_back(dynamic_cast<RFCavity*>(as.clone()));
    Component *elptr = *(--myElements.end());

    if ( (elptr->getComponentType() !="SINGLEGAP")&&(elptr->getComponentType() !="DOUBLEGAP") ) {
	    msg<<(elptr->getComponentType())<<endl;
	    throw OpalException("ParallelCyclotronTracker::visitRFCavity",
							"The ParallelCyclotronTracker can only play with cyclotron type RF system currently ...");
	}

    double rmin = (dynamic_cast<RFCavity*>(elptr))->getRmin();
    msg << "Minimal radius of cavity = " << rmin << " [mm]" << endl;

    double rmax = (dynamic_cast<RFCavity*>(elptr))->getRmax();
    msg << "Maximal radius of cavity = " << rmax << " [mm]" << endl;

    double rff = (dynamic_cast<RFCavity*>(elptr))->getCycFrequency();
    msg << "RF frequency (2*pi*f) = " << rff << " [rad./s]" << endl;

    string fmfn = (dynamic_cast<RFCavity*>(elptr))->getFieldMapFN();
    msg << "RF Field map file name = " << fmfn <<endl;

    double angle = (dynamic_cast<RFCavity*>(elptr))->getAzimuth();
    msg << "Cavity azimuth position = " << angle << " [deg.] " << endl;

    double gap = (dynamic_cast<RFCavity*>(elptr))->getGapWidth();
    msg << "Cavity gap width = " << gap << " [mm] " << endl;

    double pdis = (dynamic_cast<RFCavity*>(elptr))->getPerpenDistance();
    msg << "Cavity Shift distance = " << pdis << " [mm] "<< endl;


    double phi0 = (dynamic_cast<RFCavity*>(elptr))->getPhi0();
    msg << "Initial RF phase (t=0) of cavity   = " << phi0 << "[deg.] " << endl;

    // read cavity voltage profile data from file.
    elptr->initialise(itsBunch,1.0);

    double BcParameter[8];
    for (int i =0; i<8;i++) 
		BcParameter[i]=0.0;
    string ElementType = "CAVITY";
    BcParameter[0] = rmin;
    BcParameter[1] = rmax;
    BcParameter[2] = pdis;
    BcParameter[3] = angle;

    buildupFieldList(BcParameter,ElementType, elptr);
}

/** 
 * 
 * 
 * @param rfq 
 */
void ParallelCyclotronTracker::visitRFQuadrupole(const RFQuadrupole &rfq)
{
    *gmsg << "In RFQuadrupole; L= " << rfq.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RFQuadrupole*>(rfq.clone()));
}

/** 
 * 
 * 
 * @param bend 
 */
void ParallelCyclotronTracker::visitSBend(const SBend &bend)
{
    *gmsg << "In SBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<SBend*>(bend.clone()));
}

/** 
 * 
 * 
 * @param sep 
 */
void ParallelCyclotronTracker::visitSeparator(const Separator &sep)
{
    *gmsg << "In Seapator L= " << sep.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Separator*>(sep.clone()));
}

/** 
 * 
 * 
 * @param sept 
 */
void ParallelCyclotronTracker::visitSeptum(const Septum &sept)
{
    Inform msg("visitSeptum ");

    myElements.push_back(dynamic_cast<Septum*>(sept.clone()));
    Component *elptr = *(--myElements.end());

    double xstart = (dynamic_cast<Septum*>(elptr))->getXstart();
    msg << "XStart= " << xstart << " [mm]" << endl;

    double xend = (dynamic_cast<Septum*>(elptr)) ->getXend();
    msg << "XEnd= " << xend << " [mm]" << endl;

    double ystart =  (dynamic_cast<Septum*>(elptr))->getYstart();
    msg << "YStart= " << ystart<< " [mm]" << endl;

    double yend = (dynamic_cast<Septum*>(elptr)) ->getYend();
    msg << "YEnd= " << yend<< " [mm]" << endl;

    double width = (dynamic_cast<Septum*>(elptr)) ->getWidth();
    msg << "Width= " << width<< " [mm]" << endl;

   
    // initialise, do nothing
    elptr->initialise(itsBunch,1.0);

    double BcParameter[8];
    for (int i =0; i<8;i++) 
		BcParameter[i]=0.0;
    string ElementType = "SEPTUM";
    BcParameter[0] =xstart ;
    BcParameter[1] =xend;
    BcParameter[2] =ystart ;
    BcParameter[3] =yend;
    BcParameter[4] =width ;

    // store septum parameters in the list
    buildupFieldList(BcParameter,ElementType, elptr);
}

/** 
 * 
 * 
 * @param solenoid 
 */
void ParallelCyclotronTracker::visitSolenoid(const Solenoid &solenoid)
{
    myElements.push_back(dynamic_cast<Solenoid*>(solenoid.clone()));
    Component *elptr = *(--myElements.end());
    if (!elptr->hasAttribute("ELEMEDGE")) {
	    *gmsg << "Solenoid: no position of the element given!" << endl;
	    return;
	}
	
    double startField = elptr->getAttribute("ELEMEDGE");
    double endField;

    //buildupFieldList(startField, endField, elptr);

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
                                                   const BMultipoleField &field, double scale)
{

}

/** 
 * not used
 * 
 * @param angle 
 * @param curve 
 * @param field 
 * @param scale 
 */
void ParallelCyclotronTracker::applyExitFringe(double angle, double curve,
                                               const BMultipoleField &field, double scale)
{

}


/** 
 * 
 * 
 * @param BcParameter 
 * @param ElementType 
 * @param elptr 
 */
void ParallelCyclotronTracker::buildupFieldList(double BcParameter[], string ElementType, Component* elptr)
{
    beamline_list::iterator sindex;

    type_pair* localpair = new type_pair();
    localpair->first = ElementType;

    for(int i=0; i<8; i++ ) 
		*(((localpair->second).first)+i)=*(BcParameter+i);

    (localpair->second).second = elptr;
	
    // always put cyclotron as the first element in the list.
    if (ElementType == "CYCLOTRON"){
		sindex = FieldDimensions.begin();
    } 
	else {
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
void ParallelCyclotronTracker::visitBeamline(const Beamline & bl){
    itsBeamline->iterate(*dynamic_cast<BeamlineVisitor*>(this),false);
}

/** 
 * 
 * 
 * @param zin 
 * @param f 
 * @param gradf 
 * @param ds 
 * @param nx 
 * 
 * @return 
 */
static Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx)
{

    Vector zf;
    return zf;
}

/** 
 * 
 * 
 * @param zin 
 * @param f 
 * @param s 
 * @param ds 
 * @param nx 
 * @param cx 
 * 
 * @return 
 */
static Vector implicitInt2(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{

    Vector zf=zin;
    return zf;
}

/** 
 * 
 * 
 * @param zin 
 * @param f 
 * @param s 
 * @param ds 
 * @param nx 
 * @param cx 
 * 
 * @return 
 */
static Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
    Vector zf=zin;
    return zf;
}

/** 
 * 
 * 
 * @param zin 
 * @param f 
 * @param ds 
 * @param nx 
 * 
 * @return 
 */
static Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx)
{
    Vector zf;
    return zf;
}

/** 
 * 
 * 
 * @param zin 
 * @param f 
 * @param s 
 * @param ds 
 * @param nx 
 * @param cx 
 * 
 * @return 
 */
static Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx)
{
    Vector zf=zin;
    return zf;
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
    minnlp=0;
    maxnlp=111111;
    reduce(nlp,minnlp,OpMinAssign());
    reduce(nlp,maxnlp,OpMaxAssign());
}

/** 
 * 
 * 
 */
void ParallelCyclotronTracker::execute(){
  
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
	*gmsg<<"-----------------------------"<<endl;
	*gmsg<<"The selected Beam line elements are :"<<endl;
	for ( beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++)
		*gmsg<<((*sindex)->first)<<endl;
	*gmsg<<"-----------------------------"<<endl;
	
	// external field arrays for dumping
	for(int k=0; k<2; k++) 
		FDext_m[k] = Vector_t(0.0,0.0,0.0);
	extE_m = Vector_t(0.0, 0.0, 0.0);
	extB_m = Vector_t(0.0, 0.0, 0.0);
	
	// mark the dumpstep to inject new bunch from here for AUTO mode of restart run of multibunch
	backupDumpStep_m = -1;

	if (timeIntegrator_m == 0) {
		*gmsg<<"* 4th order Runge-Kuta integrator"<<endl;
		Tracker_RK4();
	}
	else if(timeIntegrator_m == 1) {
		*gmsg<<"* 2nd order Leap-Frog integrator"<<endl;
		Tracker_LF();
	}
	else {
		*gmsg<<"ERROR: Invalid name of TIMEINTEGRATOR in Track command"<<endl;
		exit(1);
    }

	*gmsg<<"-----------------------------"<<endl;
	*gmsg<<"Finalize i.e. write data and close files :"<<endl;
	for ( beamline_list::iterator sindex = FieldDimensions.begin(); sindex != FieldDimensions.end(); sindex++) {
	    ( ((*sindex)->second).second )->finalise();
	}
	*gmsg<<"-----------------------------"<<endl;
}

/**
   In general the two tracker have much code in common.
   This is a great source of erros.
   Need to avoid this

*/



/** 
 * 
 * 
 */
void ParallelCyclotronTracker::Tracker_LF(){

	Inform *gmsgAll;
	gmsgAll = new  Inform("CycTracker LF",INFORM_ALL_NODES);

	BorisPusher pusher(itsReference);
 
	// time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();

    // determine how many sample particles used for tuning with SC
    const int NoSamplePart = 500;
    // determine when to record tune data;  
    const int StartTurnnumber = 30;

	beamline_list::iterator sindex = FieldDimensions.begin();
	const double harm = (((*sindex)->second).second )-> getCyclHarm();

	// load time
	double t  = itsBunch->getT()*1.0e9;
	const double dt = itsBunch->getdT()*1.0e9*harm; //[s]-->[ns]

	// find the injection time interval
	if (numBunch_m > 1){
		const double RfFreq = (((*sindex)->second).second )-> getRfFrequ(); // [MHz]
		*gmsg << "Time interval between neighbour bunches is set to "<<stepsPerTurn*dt<<"[ns]"<<endl;
	}
	
	/* 
	   transform initial coordinates and momenta from local beam frame (relative value to the bunchcenter )
	   to global Cartesian frame (absolute value ).
	   Do some other initialisation work.
	   for multi-particle Mode and single Particle Mode, output particles of ID = 0 and 1 for each  dumpfreq steps
	*/
	
	string SfileName = OPAL.getInputFn();
	int pdot = SfileName.find(string("."),0);
	SfileName.erase(pdot,SfileName.size()-pdot);

	string  SfileName1 = SfileName+string("-trackOrbit.dat");
	ofstream outf;
	outf.setf(ios::scientific, ios::floatfield );
	outf.precision(8);

	if (initialTotalNum_m > 2)
		if ( myNode_m == 0 &&  flagDoTune_m){		
			// FixMe: this can go into a function
			ifstream inf;
			char skipChar[100];
			double skipNum;
			Vector_t tempSEO;
			string  SfileName5 = string("S")+SfileName+string("-trackOrbit.dat");
			inf.open(SfileName5.c_str(),ios::in);
			if(inf.fail()) {
				*gmsg<<"Cannot open file "<<SfileName5<<" for tune calculation, please check if it really exists."<<endl;
				exit(1);
			}

			if(!inf.eof()){
				for(int i=0; i<14; i++) {
					inf>>skipChar;
				}
			}
			
			while(!inf.eof()){
				inf>>skipChar>> tempSEO(1) >> skipNum >> tempSEO(2) >> skipNum >> skipNum >> skipNum;
				if(inf.eof()) break;
				tempSEO(0) = calculateAngle(tempSEO(1),tempSEO(2));
				variable_SEO_m.push_back(tempSEO);
			}
			inf.close();
			*gmsg<<"Finish reading in SEO file for tune calculation"<<endl;
		}

	// FixMe: make new element striper 

	// flag and position of the striper for multi-bunch extraction
	
	const bool flagExtraction= false;
	// const double R_stripper=400.0;  //R mm
	// const double angle_stripper=45.0/180.0*pi; //theta
	//  const double R_stripper=1852.30;  //R mm
	// const double angle_stripper=63.41247/180.0*pi; //theta
	const double R_stripper=1875.5;  //R mm
	const double angle_stripper=59.63/180.0*pi; //theta
	
	// parameter for reset bin in multi-bunch run, todo: readin from inputfile
	const  double eta = 0.01; 
  
	int SteptoLastInj = itsBunch->getSteptoLastInj();
	int restart_readstep;
  
	// get data from h5 file for restart run
	if(OPAL.inRestartRun() && numBunch_m > 1){
		if(flagExtraction) {
			FILE *f=NULL;
			string  temp_fileName = SfileName+string("-readinfrom.in");
			if( (f=fopen(temp_fileName.c_str(),"r")) == NULL ) {
				*gmsg<<" Error in ParallCyclotronTracker cannot open file "<< temp_fileName <<" please check."<<endl;
				exit(1);
			}
			assert(fscanf(f, "%d", &restart_readstep));
			*gmsg <<"read in bunch phase space from step "<<restart_readstep<<endl;
			fclose(f);
		}
		
		restartStep0_m = itsBunch->getTrackStep();
		step_m = restartStep0_m;
		itsBunch->resetPartBinID2(eta);
		*gmsg <<"Restart at integration step "<<restartStep0_m <<endl;
	}
	
	if (OPAL.hasBunchAllocated() && Options::scan ) {
		lastDumpedStep_m = 0;
		t = 0.0;
	}
  
	*gmsg<< "* Beginning of this run is at t= "<< t << " [ns]"<<endl;
	*gmsg<< "* The time step is set to dt= "<<dt << " [ns]"<<endl;

	// add header in the data dump file for particle ID equal 0 and 1
	if(myNode_m == 0) {
		if(OPAL.inRestartRun()){
			outf.open(SfileName1.c_str(),ios::app);
			outf      << "# Restart at "<<step_m<<" dumping step"<<endl;
		}
		else{
			outf.open(SfileName1.c_str());
			outf      <<"# ID   x [mm]          px [rad]       y [mm]          py [rad]        z [mm]          pz [rad]"
					  << endl;
		}
	}
	
	// for single Particle Mode, output at zero degree.
	if (initialTotalNum_m == 1)
		openFiles(SfileName);


	if(!OPAL.inRestartRun()) {

		PathLength_m = 0.0;
		// Force the initial phase space values of the particle with ID=0 to zero, to set it as a reference particle.
		if (initialTotalNum_m > 2) { // only for mulit-particle mode
			for(int ii=0; ii<initialLocalNum_m; ii++) {
				if (itsBunch->ID[ii] == 0) {
					itsBunch->R[ii] = Vector_t(0.0);
					itsBunch->P[ii] = Vector_t(0.0);
				}
				// for tuning calculation
				if(flagDoTune_m){
					if (itsBunch->ID[ii] == 1) {
						itsBunch->R[ii] = Vector_t(0.003, 0.0, 0.003);
						itsBunch->P[ii] = Vector_t(0.0);
					}
				}
			}
		}
		
		if ( Options::psDumpLocalFrame ) {

			// dump the initial distribution
			lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
			*gmsg <<"phase space dumping in local frame at initial position "<<endl;
			*gmsg << "Dump step = "<< lastDumpedStep_m <<endl;
		}


		itsBunch->R *= Vector_t(1000.0); // m --> mm
		for (int ii = 0; ii < initialLocalNum_m; ii++) {

			// rotate the local beam frame for angle of minus theta
			double tempDx = itsBunch->R[ii](0) * cosRefTheta - itsBunch->R[ii](1) * sinRefTheta;
			double tempDy = itsBunch->R[ii](0) * sinRefTheta + itsBunch->R[ii](1) * cosRefTheta;
			// move the origin point of global frame
			itsBunch->R[ii](0) = tempDx + referenceR * cosRefTheta;
			itsBunch->R[ii](1) = tempDy + referenceR * sinRefTheta;

			double tempPx = ( itsBunch->P[ii](0) + referencePr) * cosRefTheta -  ( itsBunch->P[ii](1) + referencePt) * sinRefTheta;
			double tempPy = ( itsBunch->P[ii](0) + referencePr) * sinRefTheta +  ( itsBunch->P[ii](1) + referencePt) * cosRefTheta;

			itsBunch->P[ii](0)=tempPx;
			itsBunch->P[ii](1)=tempPy;
			// initialize the bin number of the first bunch to 0
			itsBunch->Bin[ii]=0;

		}

		if ( !(Options::psDumpLocalFrame) ) {
			
			// dump the initial distribution
			itsBunch->R /= Vector_t(1000.0); // mm --> m
			lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
			itsBunch->R *= Vector_t(1000.0); // m --> mm
			*gmsg<<"meanR=( "<<referenceR * cosRefTheta<<" "<<referenceR * sinRefTheta<<" "<<0.0<<" ) [mm] "<<endl;
			*gmsg <<"phase space dumping in glabol frame at initial position "<<endl;
			*gmsg << "Dump step = "<< lastDumpedStep_m <<endl;
		}

		// AUTO mode
		if (multiBunchMode_m == 2 ) {
			Vector_t Rmean = itsBunch->get_centroid()*1000.0; // m->mm
			RLastTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );
			RThisTurn_m = RLastTurn_m;
			
			*gmsg<<"initial Radial position = "<<RThisTurn_m<<" [mm]"<<endl;
		}

		/*
		//todo: more generic initial azimuthal position, not finished yet
		double referencePx = referencePr*cosRefTheta - referencePt*sinRefTheta;
		double referencePy = referencePr*sinRefTheta + referencePt*cosRefTheta;
		double referencePhi = calculateAngle(referencePx,referencePy);
		
		referencePhi -= pi/2.0;
		// end todo
		*/
		
		if ((initialTotalNum_m > 2 ) && ( numBunch_m > 1) && (multiBunchMode_m == 1) ){

			// backup initial distribution
			assert( initialR_m = new Vector_t[initialLocalNum_m]);
			assert( initialP_m = new Vector_t[initialLocalNum_m]);

			for (int ii=0; ii<initialLocalNum_m; ii++){
				initialR_m[ii] = itsBunch->R[ii];
				initialP_m[ii] = itsBunch->P[ii];
			}
		}
		
	}else{
		PathLength_m = itsBunch->getLPath();
		
		// AUTO mode
		if (multiBunchMode_m == 2 ) {
			itsBunch->calcBeamParameters_cycl();
			
			Vector_t Rmean = itsBunch->get_centroid()*1000.0; // m->mm
			RLastTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );
			RThisTurn_m = RLastTurn_m;
			
			*gmsg<<"Radial position at restart position = "<<RThisTurn_m<<" [mm]"<<endl;
		}
		
		itsBunch->R *= Vector_t(1000.0); // m --> mm
		
	}

	//  read in some control parameters
	const int SinglePartDumpFreq = Options::sptDumpFreq;
	const int RepartitionFreq = Options::repartFreq;
	const int resetBinFreq = Options::rebinFreq;
	const int scSolveFreq = Options::scSolveFreq;
	const bool doDumpAfterEachTurn =Options::psDumpEachTurn;
    
	
	int boundpDestroyFreq = 10; // todo: is it better treat as a control parameter
	
	// prepare for dump after each turn
	const double initialReferenceTheta = referenceTheta/180.0*pi;
	double oldReferenceTheta = initialReferenceTheta;

	*gmsg<<"single particle trajectory dump frequency is set to "<<SinglePartDumpFreq<<endl;
	*gmsg<<"particles repartition frequency is set to "<<RepartitionFreq<<endl;
	if(numBunch_m > 1)
		*gmsg<<"particles energy bin ID reset frequency is set to "<<resetBinFreq<<endl;

	// if initialTotalNum_m = 2, trigger SEO mode
	// prepare for transverse tuning calculation
	vector<double> Ttime,Tdeltr,Tdeltz;
	// if flagDoTune == tune for multi-particles mode
	// prepare for transverse tuning calculation
	vector<double> TdeltrSC[NoSamplePart],TdeltzSC[NoSamplePart];
	vector<int> TturnNumber;
	int turnnumber = 1;
	int lastTurn = 1;


	bool flagNoDeletion;

	// flag to determine when to transit from single-bunch to multi-bunches mode
	bool flagTransition = false;
	// flag to omit the dump step after backup dump step for multi-bunches mode
	bool OneStepAfterBackUp = false;
	// step point determining the next time point of check for transition
	int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();

	const  double deltaTheta = pi/(stepsPerTurn);
	// record at which angle the space charge are solved
	double angleSpaceChargeSolve;
	// local particle number
	int nlp;
	double changeunite = 1000.0*c*dt*1.0e-9;  
  
	if(initialTotalNum_m == 1){
		*gmsg << "***---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** "<<endl;
		*gmsg << "Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," <<endl
			  << " The initial distribution file must be specified which should contain only one line for the single particle "<<endl
			  << "***------------WARNING: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** "<<endl;
		if (Ippl::getNodes() != 1)
			throw OpalException("Error in ParallelCyclotronTracker::execute","SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");

	}else if (initialTotalNum_m == 2){
		*gmsg << "***------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** "<<endl;
		*gmsg << "Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." <<endl
			  << "This mode does NOT include any RF cavities. The initial distribution file must be specified"<<endl
			  << "In the file the first line is for reference particle and the second line is for offcenter particle."<<endl
			  << "The tuning is calculated by FFT routines based on these two particles. "<<endl
			  << "***------------WARNING: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** "<<endl;
		if (Ippl::getNodes() != 1)
			throw OpalException("Error in ParallelCyclotronTracker::execute","SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
	}else{
  
    	*gmsg << "***---------------------------- MULTI-PARTICLES MODE------ ----------------------------*** "<<endl;
		nlp=itsBunch->getLocalNum();
		/**
		   Need itsBunch->boundpNoRep();
		   otherwise dead lock!
		*/
		// itsBunch->boundpNoRep();
	
		nlp=itsBunch->getLocalNum();
		checkNumPart(string("Before boundp "), nlp);

		itsBunch->boundp();
		nlp=itsBunch->getLocalNum();

		checkNumPart(string("After boundp"),nlp);
	}

    // do repartation for restart run
	if(OPAL.inRestartRun() && numBunch_m > 1){

		Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
		Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);

		for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
			for (int jj=0; jj<3; jj++){
				temp_meanR(jj) += itsBunch->R[ii](jj);
			}
		}
    
		// sum all particles in the same bin using reduce() for parallel enviroment,
		// and change N to initialTotalNum_m.
		reduce(temp_meanR,temp_meanR,OpAddAssign());
		temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
    
		for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
			for (int jj=0; jj<3; jj++){
				temp_meanP(jj) += itsBunch->P[ii](jj);
			}
		}
    
		reduce(temp_meanP,temp_meanP,OpAddAssign());
		temp_meanP /= Vector_t(itsBunch->getTotalNum());
    
		// in global Cartesian frame, calculate the direction of longitudinal angle of bunch
    
		double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));
    
		// the angle from positive radial direction to propagational  direction
		// the angle from positive azimuthal direction to propagational direction
		temp_meanPhi -= pi/2.0;
		double angleSpaceChargeSolve = temp_meanPhi;
    
		double cosTemp_meanPhi = cos(temp_meanPhi);
		double sinTemp_meanPhi = sin(temp_meanPhi);
    
		double temp_meanPLongitudinal2 = pow(temp_meanP(0),2.0) + pow(temp_meanP(1),2.0);
    
		double temp_meangamma = sqrt(1.0 +temp_meanPLongitudinal2);
    
		// remove mean coordinates
		itsBunch->R -= temp_meanR;
    
		//scale coordinates
		itsBunch->R /= Vector_t(1000.0); // mm --> m
    
		// rotate from global frame to local reference frame(transverse horizontal, longitudinal, transverse vertical )
		// For single bin, rotate the frame for temp_meanPhi degree
		for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
			double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  + itsBunch->R[ii](1)*sinTemp_meanPhi;
			double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
      
			itsBunch->R[ii](0) = temp_RHorizontal;
			itsBunch->R[ii](1) = temp_RLongitudinal;
		}
    
		itsBunch->boundp();

		*gmsg<< "Do binaryRepartation for restart run"<<endl;
		nlp=itsBunch->getLocalNum();
		checkNumPart(string("Before repartation: "), nlp);
		IpplTimings::startTimer(BinRepartTimer_m);
		itsBunch->do_binaryRepart();
		IpplTimings::stopTimer(BinRepartTimer_m);
		Ippl::Comm->barrier();
		nlp=itsBunch->getLocalNum();
		checkNumPart(string("After repartation: "), nlp);
		*gmsg<< "Done binaryRepartation for restart run"<<endl;
        
		// HERE transform particles coordinates back to global frame (rotate and shift)
		// rotate back from local reference frame to global frame
		for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
			double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  - itsBunch->R[ii](1)*sinTemp_meanPhi;
			double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
      
			itsBunch->R[ii](0) = temp_RHorizontal;
			itsBunch->R[ii](1) = temp_RLongitudinal;
		}
    
		//scale coordinates back
		itsBunch->R *= Vector_t(1000.0); // m --> mm
    
		// retrieve mean coordinates
		itsBunch->R += temp_meanR;
	}
    
	// *****************II***************
	// main integration loop
	// *****************II***************
	*gmsg << "***---------------------------- Start tracking   ------------------------------------------*** "<<endl;
	for(step_m; step_m < maxSteps_m; step_m++){
		bool dumpEachTurn = false;
		bool flagNeedUpdate = false;
		// *****************II-1***************
		// dump The particle with ID =0 & 1
		// *****************II-1***************
		if ( step_m%SinglePartDumpFreq == 0 ){
			IpplTimings::startTimer(DumpTimer_m);

			double x;
			int  id;
			vector<double> tmpr;
			vector<int> tmpi;

			int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
      
			// for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
			int found[2] = {-1,-1};
			int counter = 0;

			for(int ii=0; ii<( itsBunch->getLocalNum()); ii++) {
				if (itsBunch->ID[ii] == 0) {
					found[counter] = ii;
					counter++;
				}
				if (itsBunch->ID[ii] == 1) {
					found[counter] = ii;
					counter++;
				}
			} 
      
			if(myNode_m == 0) {
				// for root node
				int notReceived =  Ippl::getNodes() - 1;
				int numberOfPart =0;

				while (notReceived > 0) {
					int node = COMM_ANY_NODE;
					Message* rmsg =  Ippl::Comm->receive_block(node, tag);
					if (rmsg == 0)
						ERRORMSG("Could not receive from client nodes in main." << endl);
					notReceived--;
					rmsg->get(&numberOfPart);
					for (int ii=0; ii < numberOfPart; ii++) {
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
				for (int ii=0; ii < counter; ii++) {
					tmpi.push_back(itsBunch->ID[found[ii]]);
					for (int jj=0;jj<3;jj++) {
						tmpr.push_back(itsBunch->R[found[ii]](jj));
						tmpr.push_back(itsBunch->P[found[ii]](jj));
					}
				}
				vector<double>::iterator itParameter = tmpr.begin();
				vector<int>::iterator  itId = tmpi.begin();

				for ( itId = tmpi.begin(); itId != tmpi.end(); itId++ ){
					outf<<"ID"<<*itId;
					for (int ii = 0;ii < 6; ii++ ){
						outf<<" "<<*itParameter;
						itParameter++;
					}
					outf <<endl;
				}
         
			}
			else {
				// for other nodes
				Message* smsg = new Message();
				smsg->put(counter);
				for (int ii=0; ii < counter; ii++) {
					smsg->put(itsBunch->ID[found[ii]]);
					for (int jj=0;jj<3;jj++) {
						smsg->put(itsBunch->R[found[ii]](jj));
						smsg->put(itsBunch->P[found[ii]](jj));
					}
				}
				bool res = Ippl::Comm->send(smsg, 0, tag);
				if (!res)
					ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
			}
      
			IpplTimings::stopTimer(DumpTimer_m);
		}
		//end dump
		Ippl::Comm->barrier();
    
		IpplTimings::startTimer(IntegrationTimer_m);
		// push particles for first half setp
		for (int i = 0; i < itsBunch->getLocalNum(); ++i) {  //f1
		        rold_m = itsBunch->R[i];
			pold_m = itsBunch->P[i];
				  
			itsBunch->R[i] /= changeunite ; //mm --> dimensionless, [cdt]
			pusher.push(itsBunch->R[i], itsBunch->P[i], dt*1.0e-9);
			itsBunch->R[i] *= changeunite; //  dimensionless, [cdt] --> mm
				  
			rnew_m = itsBunch->R[i];
    
			/* check whether gap crossing happens during the first half push step, 
			   restore the old data of this particle before the first half push step,
			   and (1) drift for (dt1),(2) do momenta kicking (3) drift step for (dt/2-dt1)
			*/
			for ( beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++){
				// here check gap cross in the list, if do , set tag_crossing to TRUE
				bool tag_crossing = false;
				double DistOld=0.0; //mm
				if( ((*sindex)->first) == "CAVITY")
					tag_crossing = checkGapCross(rold_m, rnew_m, ((*sindex)->second).second, DistOld);
	
				if(tag_crossing){
					// gap crossing happens
					double oldMomentum2 = dot(pold_m,pold_m); 
					double oldBetgam = sqrt(oldMomentum2);
					double oldGamma = sqrt(1.0 +oldMomentum2);
					double oldBeta = oldBetgam/oldGamma;
					double dt1 = DistOld/( c* oldBeta * 1.0e-6);  // ns
					double dt2 = dt/2 - dt1;
					double tempP[3];
					const double temp_RscaleFactor1 = c*dt1*1.0e-6/oldGamma; 
					const double temp_RscaleFactor2 = c*dt2*1.0e-6/oldGamma; 
	  
					// retrack particle from the old postion to cavity gap point
					// restore the old coordinates and momenta

					itsBunch->R[i] = rold_m;
					itsBunch->P[i] = pold_m;

					// drift for dt1
					itsBunch->R[i] += itsBunch->P[i]*temp_RscaleFactor1;
	  
					//momentum kick
					double radius = sqrt( pow(itsBunch->R[i](0), 2.0) + pow(itsBunch->R[i](1), 2.0) - 
										  pow( (dynamic_cast<RFCavity*>(((*sindex)->second).second))->getPerpenDistance() , 2.0) );
					double rmin = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmin());
					double rmax = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmax());
					double nomalRadius = (radius - rmin)/(rmax - rmin);
					if (nomalRadius < 0 || nomalRadius > 1) 
						*gmsg<<"Current radial location of particle is less than minimal radius of cavity, "
							 << "or larger then maximal radius. "<<endl;
					else {
						for (int j =0;j<3;j++) 
							tempP[j] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
							  
						// here evaluate voltage and conduct momenta kicking note t in ns
						(dynamic_cast<RFCavity*>(((*sindex)->second).second)) -> getMomentaKick(nomalRadius, tempP, t, dt1, itsBunch->ID[i],itsReference.getM());
							  
						for (int j =0;j<3;j++) 
							itsBunch->P[i](j) = tempP[j];
					}
					// drift for dt2
					itsBunch->R[i] += itsBunch->P[i]*temp_RscaleFactor2;	
				}
			}
		} //f1 
	
		IpplTimings::stopTimer(IntegrationTimer_m);
		t += dt/2; //update for the half push
		itsBunch->setT( (t)*1.0e-9 );
		// path length update
		double temp_P2_0= pow(itsBunch->P[0](0),2.0) + pow(itsBunch->P[0](1),2.0) + pow(itsBunch->P[0](2),2.0);
		PathLength_m += c_mmtns*dt/2*(sqrt(temp_P2_0)/sqrt( 1.0 + temp_P2_0))/1000.0; // unit: m
		itsBunch->setLPath(PathLength_m);

		// bunch injection
		if (numBunch_m > 1 ){
			  
			size_t injectLocalNUM =initialTotalNum_m;
			  
			if( (BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)){
				if (step_m == stepsNextCheck ){
					// under 3 conditions, following code will be execute
					// to check the distance between two neighborring bunches
					// 1.multi-bunch mode, AUTO sub-mode
					// 2.After each revolution
					// 3.only one bunch exists

					*gmsg<<"checking for automatically injecting new bunch ..." << endl;

					itsBunch->R /= Vector_t(1000.0); // mm --> m
					itsBunch->calcBeamParameters_cycl();
					itsBunch->R *= Vector_t(1000.0); // m --> mm

					Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm

					RThisTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );

					Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

					double XYrms =  sqrt( pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0) );

					// if the distance between two neighbour bunch is less than 5 times of its 2D rms size
					// start multi-bunch simulation, fill current phase space to initialR and initialP arrays

					if ( (RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms ) {
						// since next turn, start multi-bunches
					
						itsBunch->R /= Vector_t(1000.0); // mm --> m
						lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
						itsBunch->R *= Vector_t(1000.0); // m --> mm

						backupDumpStep_m = lastDumpedStep_m;
						flagTransition = true;
						OneStepAfterBackUp =true;

						*gmsg << "In order to inject a new bunch after one revolution, Dump step = " 
							  << lastDumpedStep_m <<", time ="<<t<<" [ns]"<<endl;
						*gmsg << "**************** Another revolution later, Multi-Bunch "
							  << "Mode will be invorked automatically ****************"<<endl;
					}

					stepsNextCheck += stepsPerTurn;
	  
					*gmsg<<"RLastTurn = "<<RLastTurn_m<<" [mm]"<<endl;
					*gmsg<<"RThisTurn = "<<RThisTurn_m<<" [mm]"<<endl;
					*gmsg<<"    XYrms = "<<XYrms    <<" [mm]"<<endl;
				
					RLastTurn_m = RThisTurn_m;
				}
			} 
			else if ( SteptoLastInj == stepsPerTurn-1 ) {
				if (((!flagExtraction) && BunchCount_m < numBunch_m) || flagExtraction ){
				
					// under 4 conditions, following code will be execute
					// to read new bunch from hdf5 format file for FORCE or AUTO mode
					// 1.multi-bunch mode
					// 2.after each revolution
					// 3.existing bunches is less than the specified bunches
					// 4.FORCE mode, or AUTO mode with flagTransition = true
					// Note: restart from 1 < BunchCount < numBunch_m must be avoided.

					*gmsg<<"time = "<<t<<" [ns], Now inject a new bunch by reading data from h5 file ... ... ..."<<endl;
					BunchCount_m++;

					// read initial distribution from h5 file
					if(multiBunchMode_m == 1 ){
						readOneBunch(BunchCount_m-1, 0);
						itsBunch->resetPartBinID2(eta);
					} else if (multiBunchMode_m == 2 ){
						if (flagExtraction && OPAL.inRestartRun())
							readOneBunch(BunchCount_m-1,restart_readstep);
						else
							readOneBunch(BunchCount_m-1, backupDumpStep_m);	
						itsBunch->resetPartBinID2(eta);
					}
				
					SteptoLastInj = 0;
				
					itsBunch->setNumBunch(BunchCount_m);
				
					stepsNextCheck += stepsPerTurn;
				
					// update  after injection
					itsBunch->boundp();

					Ippl::Comm->barrier();
					*gmsg<<"read initial distribution from h5 file"<<endl;
					*gmsg<<BunchCount_m<<"'th bunch injected, total particle number = "<<itsBunch->getTotalNum()<<endl;
				}
			}
			else if ((!flagExtraction) && BunchCount_m == numBunch_m ){
				// After this, numBunch_m is wrong but useless
				numBunch_m--;
			
			}
			else{
				SteptoLastInj++;
			}
		}
    
		// FixMe: 1000 as parameter !?
		// calculate self fields Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
		if(itsBunch->hasFieldSolver() && initialTotalNum_m >=1000) {
			if(step_m % scSolveFreq == 0){
				*gmsg << "Calculate space charge at step " << step_m<<endl;
				// Firstly reset E and B to zero before fill new space charge field data for each track step
				itsBunch->Bf = Vector_t(0.0);
				itsBunch->Ef = Vector_t(0.0);
			    
				//HERE transform particles coordinates to local frame (rotate and shift)
				//
				Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
				Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);
			
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					for (int jj=0; jj<3; jj++){
						temp_meanR(jj) += itsBunch->R[ii](jj);
					}
				}
			
				// sum all particles in the same bin using reduce() for parallel enviroment,
				// and change N to initialTotalNum_m.
				reduce(temp_meanR,temp_meanR,OpAddAssign());
				temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
				if ( (itsBunch->weHaveBins()) && BunchCount_m > 1 ) {

					IpplTimings::startTimer(TransformTimer_m);

					double temp_binsMeanPhi = itsBunch->calcMeanPhi();
					// the angle from positive radial direction to propagational  direction
					// the angle from positive azimuthal direction to propagational direction
					temp_binsMeanPhi -= pi/2.0;
					angleSpaceChargeSolve = temp_binsMeanPhi;
								    
					double cosTemp_binsMeanPhi = cos(temp_binsMeanPhi);
					double sinTemp_binsMeanPhi = sin(temp_binsMeanPhi);


					// remove mean coordinates
					itsBunch->R -= temp_meanR;

					//scale coordinates
					itsBunch->R /= Vector_t(1000.0); // mm --> m

					// rotate from global frame to local frame(transverse horizontal, longitudinal, transverse vertical )
					// For multi-bin, rotate the frame for temp_binsMeanPhi degree
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_binsMeanPhi  + 
							itsBunch->R[ii](1)*sinTemp_binsMeanPhi;
						double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_binsMeanPhi  + 
							itsBunch->R[ii](1)*cosTemp_binsMeanPhi;
					
						itsBunch->R[ii](0) = temp_RHorizontal;
						itsBunch->R[ii](1) = temp_RLongitudinal;
					}

					if ( (step_m+1) % boundpDestroyFreq == 0 )
						itsBunch->boundp_destroy();
					else
						itsBunch->boundp();
				
					IpplTimings::stopTimer(TransformTimer_m);
	      
					// calcualte gamma for each energy bin
					itsBunch->calcGammas_cycl();
	      
					// repartation
					if ((step_m%RepartitionFreq) == 0 ) {
						nlp=itsBunch->getLocalNum();
						checkNumPart(string("Before repartation: "), nlp);
						IpplTimings::startTimer(BinRepartTimer_m);
						itsBunch->do_binaryRepart();
						IpplTimings::stopTimer(BinRepartTimer_m);
						Ippl::Comm->barrier();
						nlp=itsBunch->getLocalNum();
						checkNumPart(string("After repartation: "), nlp);
					}
					// calculate space charge field for each energy bin
					for (int b=0; b<itsBunch->getLastemittedBin() ; b++) {
						if( itsBunch->getNumPartInBin(b) >= 500){
							itsBunch->setBinCharge(b,itsBunch->getChargePerParticle());
							//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
							itsBunch->computeSelfFields_cycl(b);
							//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
						}
						else{
							INFOMSG("Note: Bin "<< b <<": less than 1000 particles, omit space charge fields"<<endl);
						}
					}
	      
					itsBunch->Q = itsBunch->getChargePerParticle();
				
					IpplTimings::startTimer(TransformTimer_m);
					// HERE transform particles coordinates back to global frame (rotate and shift)
					// rotate back from local reference frame to global frame
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_binsMeanPhi  - 
							itsBunch->R[ii](1)*sinTemp_binsMeanPhi;
						double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_binsMeanPhi  + 
							itsBunch->R[ii](1)*cosTemp_binsMeanPhi;
						itsBunch->R[ii](0) = temp_RHorizontal;
						itsBunch->R[ii](1) = temp_RLongitudinal;
					}

					//scale coordinates back
					itsBunch->R *= Vector_t(1000.0); // m --> mm

					// retrieve mean coordinates
					itsBunch->R += temp_meanR;
				
					// HERE transform self field back to global frame (rotate)
				
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_binsMeanPhi  - itsBunch->Ef[ii](1)*sinTemp_binsMeanPhi;
						double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_binsMeanPhi  + itsBunch->Ef[ii](1)*cosTemp_binsMeanPhi;
					
						itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
						itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
					}

					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_binsMeanPhi  - itsBunch->Bf[ii](1)*sinTemp_binsMeanPhi;
						double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_binsMeanPhi  + itsBunch->Bf[ii](1)*cosTemp_binsMeanPhi;
					
						itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
						itsBunch->Bf[ii](1) = temp_E2;  // By,T
					}
				
				}
				else{
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						for (int jj=0; jj<3; jj++){
							temp_meanP(jj) += itsBunch->P[ii](jj);
						}
					}
				
					reduce(temp_meanP,temp_meanP,OpAddAssign());
					temp_meanP /= Vector_t(itsBunch->getTotalNum());
				
					// in global Cartesian frame, calculate the direction of longitudinal angle of bunch

					double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));
				
					// the angle from positive radial direction to propagational  direction
					//                                   ==>
					// the angle from positive azimuthal direction to propagational direction
					temp_meanPhi -= pi/2.0;
					angleSpaceChargeSolve = temp_meanPhi;
								
					double cosTemp_meanPhi = cos(temp_meanPhi);
					double sinTemp_meanPhi = sin(temp_meanPhi);
				
					double temp_meanPLongitudinal2 = pow(temp_meanP(0),2.0) + pow(temp_meanP(1),2.0);

					double temp_meangamma = sqrt(1.0 +temp_meanPLongitudinal2);
	    
					// remove mean coordinates
					itsBunch->R -= temp_meanR;

					//scale coordinates
					itsBunch->R /= Vector_t(1000.0); // mm --> m

					// rotate from global frame to local frame(transverse horizontal, longitudinal, transverse vertical )
					// For single bin, rotate the frame for temp_meanPhi degree
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  + itsBunch->R[ii](1)*sinTemp_meanPhi;
						double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
					
						itsBunch->R[ii](0) = temp_RHorizontal;
						itsBunch->R[ii](1) = temp_RLongitudinal;
					}

					if ( (step_m+1) % boundpDestroyFreq == 0 )
						itsBunch->boundp_destroy();
					else
						itsBunch->boundp();

					IpplTimings::stopTimer(TransformTimer_m);

					if ((step_m%RepartitionFreq) == 0 ) {
						nlp=itsBunch->getLocalNum();
						checkNumPart(string("Before repartation: "), nlp);
						IpplTimings::startTimer(BinRepartTimer_m);
						itsBunch->do_binaryRepart();
						IpplTimings::stopTimer(BinRepartTimer_m);
						Ippl::Comm->barrier();
						nlp=itsBunch->getLocalNum();
						checkNumPart(string("After repartation: "), nlp);
					}
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
					itsBunch->computeSelfFields_cycl( temp_meangamma );
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
				
					IpplTimings::startTimer(TransformTimer_m);
					// HERE transform particles coordinates back to global frame (rotate and shift)
					// rotate back from local reference frame to global frame
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  - itsBunch->R[ii](1)*sinTemp_meanPhi;
						double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
						itsBunch->R[ii](0) = temp_RHorizontal;
						itsBunch->R[ii](1) = temp_RLongitudinal;
					}

					//scale coordinates back
					itsBunch->R *= Vector_t(1000.0); // m --> mm

					// retrieve mean coordinates
					itsBunch->R += temp_meanR;

					// HERE transform self field back to global frame (rotate)

					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_meanPhi  - itsBunch->Ef[ii](1)*sinTemp_meanPhi;
						double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_meanPhi  + itsBunch->Ef[ii](1)*cosTemp_meanPhi;
						itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
						itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
					}
				
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_meanPhi  - itsBunch->Bf[ii](1)*sinTemp_meanPhi;
						double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_meanPhi  + itsBunch->Bf[ii](1)*cosTemp_meanPhi;
						itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
						itsBunch->Bf[ii](1) = temp_E2;  // By,T
					}
				}
				if (step_m%SinglePartDumpFreq == 0) {

				}
				IpplTimings::stopTimer(TransformTimer_m);
			}
			else {
				Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);

				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					for (int jj=0; jj<3; jj++){
						temp_meanP(jj) += itsBunch->P[ii](jj);
					}
				}

				reduce(temp_meanP,temp_meanP,OpAddAssign());
				temp_meanP /= Vector_t(itsBunch->getTotalNum());
			    
				// in global Cartesian frame, calculate the direction of longitudinal angle of bunch
			    
				double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));
	  
				// the angle from positive radial direction to propagational  direction
				//                                   ==>
				// the angle from positive azimuthal direction to propagational direction
				temp_meanPhi -= pi/2.0;
				double temp_deltaPhi = temp_meanPhi-angleSpaceChargeSolve;
			    
				double cosTemp_deltaPhi = cos(temp_deltaPhi);
				double sinTemp_deltaPhi = sin(temp_deltaPhi);
			    			    
	  
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_deltaPhi  - itsBunch->Ef[ii](1)*sinTemp_deltaPhi;
					double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_deltaPhi  + itsBunch->Ef[ii](1)*cosTemp_deltaPhi;
					itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
					itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
				}
	  
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_deltaPhi  - itsBunch->Bf[ii](1)*sinTemp_deltaPhi;
					double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_deltaPhi  + itsBunch->Bf[ii](1)*cosTemp_deltaPhi;
					itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
					itsBunch->Bf[ii](1) = temp_E2;  // By,T
				}

			}
		}
		else {
			// if field solver is not available , only update bunch, to transfer particles between nodes if needed,
			// reset parameters such as LocalNum, initialTotalNum_m.
			// INFOMSG("No space charge Effects are included!"<<endl;);
			if ((step_m%RepartitionFreq*100) == 0 && initialTotalNum_m >=1000) {
				Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
				Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);

				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					for (int jj=0; jj<3; jj++){
						temp_meanR(jj) += itsBunch->R[ii](jj);
					}
				}
	  
				// sum all particles in the same bin using reduce() for parallel enviroment,
				// and change N to initialTotalNum_m.
				reduce(temp_meanR,temp_meanR,OpAddAssign());
				temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
	  
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					for (int jj=0; jj<3; jj++){
						temp_meanP(jj) += itsBunch->P[ii](jj);
					}
				}
			
				reduce(temp_meanP,temp_meanP,OpAddAssign());
				temp_meanP /= Vector_t(itsBunch->getTotalNum());
				
				// in global Cartesian frame, calculate the direction of longitudinal angle of bunch

				double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));
				
				// the angle from positive radial direction to propagational  direction
				//                                   ==>
				// the angle from positive azimuthal direction to propagational direction
				temp_meanPhi -= pi/2.0;
				angleSpaceChargeSolve = temp_meanPhi;
							
				double cosTemp_meanPhi = cos(temp_meanPhi);
				double sinTemp_meanPhi = sin(temp_meanPhi);

				double temp_meanPLongitudinal2 = pow(temp_meanP(0),2.0) + pow(temp_meanP(1),2.0);

				double temp_meangamma = sqrt(1.0 +temp_meanPLongitudinal2);
	    
				// remove mean coordinates
				itsBunch->R -= temp_meanR;

				//scale coordinates
				itsBunch->R /= Vector_t(1000.0); // mm --> m

				// rotate from global frame to local frame(transverse horizontal, longitudinal, transverse vertical )
				// For single bin, rotate the frame for temp_meanPhi degree
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  + itsBunch->R[ii](1)*sinTemp_meanPhi;
					double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;

					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;
				}

				if ( (step_m+1) % boundpDestroyFreq == 0 )
					itsBunch->boundp_destroy();
				else
					itsBunch->boundp();
			
				if ((step_m%RepartitionFreq) == 0 ) {
					nlp=itsBunch->getLocalNum();
					checkNumPart(string("Before repartation: "), nlp);
					IpplTimings::startTimer(BinRepartTimer_m);
					itsBunch->do_binaryRepart();
					IpplTimings::stopTimer(BinRepartTimer_m);
					Ippl::Comm->barrier();
					nlp=itsBunch->getLocalNum();
					checkNumPart(string("After repartation: "), nlp);
				}
	    	    
				// HERE transform particles coordinates back to global frame (rotate and shift)
				// rotate back from local reference frame to global frame
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  - itsBunch->R[ii](1)*sinTemp_meanPhi;
					double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;
				}
	    
				//scale coordinates back
				itsBunch->R *= Vector_t(1000.0); // m --> mm

				// retrieve mean coordinates
				itsBunch->R += temp_meanR;
	  
			}
		}

		IpplTimings::startTimer(IntegrationTimer_m);
		//  kick particles for one step
		for (int i = 0; i < itsBunch->getLocalNum(); ++i) {
			Vector_t externalE, externalB;
			double partR;

			externalB = Vector_t(0.0, 0.0, 0.0);
			externalE = Vector_t(0.0, 0.0, 0.0);
      
			partR = sqrt(dot(itsBunch->R[i],itsBunch->R[i]));
      
			beamline_list::iterator sindex = FieldDimensions.begin();
           
			if ( ((((*sindex)->second).first)[0] <= partR ) &&  ( (((*sindex)->second).first)[1] >= partR ) ) 
				( ((*sindex)->second).second )->apply(itsBunch->R[i], Vector_t(0.0), t, externalE, externalB);
            
			externalB= externalB/10.0; // kgauss -> T
      
			if(itsBunch->hasFieldSolver()){
				externalE += itsBunch->Ef[i];
				externalB += itsBunch->Bf[i];
			}
			pusher.kick(itsBunch->R[i], itsBunch->P[i],externalE , externalB, dt*1.0e-9);
		}
    
		// push particles for second half setp - f1
		for (int i = 0; i < itsBunch->getLocalNum(); ++i) {  
			rold_m =  itsBunch->R[i]; // used for gap cross checking
			pold_m =  itsBunch->P[i];
	
			itsBunch->R[i] /= changeunite ; //mm --> dimensionless, [cdt]
			pusher.push(itsBunch->R[i], itsBunch->P[i], dt*1.0e-9);
			itsBunch->R[i] *= changeunite; //  dimensionless, [cdt] --> mm
      		
			rnew_m = itsBunch->R[i];  //[x,y,z]  units: [mm]
    
			// check whether gap crossing happens during the first half push step, 
			// restore the old data of this particle before the first half push step,
			// and (1) drift for (dt1),(2) do momenta kicking (3) drift step for (dt/2-dt1). 
			for (beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++){
			        // here check gap cross in the list, if do , set tag_crossing to TRUE
				bool tag_crossing = false;
				double DistOld=0.0; //mm
				if( ((*sindex)->first) == "CAVITY")
					tag_crossing = checkGapCross(rold_m, rnew_m, ((*sindex)->second).second, DistOld);
	
				if(tag_crossing){
					// gap crossing happens 
					double oldMomentum2  = dot(pold_m,pold_m);
					double oldBetgam = sqrt(oldMomentum2);
					double oldGamma = sqrt(1.0 +oldMomentum2);
					double oldBeta = oldBetgam/oldGamma;
					double dt1 = DistOld/( c* oldBeta * 1.0e-6);  // ns
					double dt2 = dt/2 - dt1;
					double tempP[3];
					const double temp_RscaleFactor1 = c*dt1*1.0e-6/oldGamma; 
					const double temp_RscaleFactor2 = c*dt2*1.0e-6/oldGamma; 
	  
					// retrack particle from the old postion to cavity gap point
					// restore the old coordinates and momenta
					itsBunch->R[i] = rold_m;
					itsBunch->P[i] = pold_m;
					  
					// drift for dt1
					itsBunch->R[i] += itsBunch->P[i]*temp_RscaleFactor1;
					  
					//momentum kick
					double radius = sqrt( pow(itsBunch->R[i](0), 2.0) + pow(itsBunch->R[i](1), 2.0) - pow( (dynamic_cast<RFCavity*>(((*sindex)->second).second))->getPerpenDistance() , 2.0) );
					double rmin = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmin());
					double rmax = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmax());
					double nomalRadius = (radius - rmin)/(rmax - rmin);
					if (nomalRadius < 0 || nomalRadius > 1) {
						*gmsg<<"Current radial location of particle is less than minimal radius of cavity, or larger then maximal radius. " <<endl;
					} 
					else {
						for (int j =0;j<3;j++) tempP[j] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
						// here evaluate voltage and conduct momenta kicking;
						(dynamic_cast<RFCavity*>(((*sindex)->second).second)) -> getMomentaKick(nomalRadius, tempP, t, dt1, itsBunch->ID[i],itsReference.getM()); // t : ns
						for (int j =0;j<3;j++) itsBunch->P[i](j) = tempP[j];
					}
					// drift for dt2
					itsBunch->R[i] += itsBunch->P[i]*temp_RscaleFactor2;	
				}
			} 
		}//f1 
		IpplTimings::stopTimer(IntegrationTimer_m);
    
		t += dt/2; //update for the half push
		itsBunch->setT( (t)*1.0e-9 );
		// path length update
		temp_P2_0= pow(itsBunch->P[0](0),2.0) + pow(itsBunch->P[0](1),2.0) + pow(itsBunch->P[0](2),2.0);
		PathLength_m += c_mmtns*dt/2*(sqrt(temp_P2_0)/sqrt( 1.0 + temp_P2_0))/1000.0; // unit: m
		itsBunch->setLPath(PathLength_m);
    
		// extract particles which hit particles
		if (flagExtraction){
			// extraction particles
			Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
			for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				for (int jj=0; jj<3; jj++){
					temp_meanR(jj) += itsBunch->R[ii](jj);
				}
			}
			// sum all particles in the same bin using reduce() for parallel enviroment,
			reduce(temp_meanR,temp_meanR,OpAddAssign());
			temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
			// [-pi,pi)
			double temp_meantheta = calculateAngle2(temp_meanR(0), temp_meanR(1));
		
			if (abs(temp_meantheta - angle_stripper) < 200*deltaTheta) {
				bool isextracted=itsBunch->boundp_extraction(R_stripper, angle_stripper, deltaTheta, t);
				// update bin stuff  
				if (isextracted && itsBunch->weHaveBins())
					itsBunch->resetPartBinID2(eta);
			}
		}
        
		// recalculate bingamma and reset the BinID for each particles according to its current gamma 
		if ( (itsBunch->weHaveBins()) && BunchCount_m > 1 ){
			if (step_m%resetBinFreq == 0 ){
				itsBunch->resetPartBinID2(eta);
			}
		}
		// dump some data after one push in single particle tracking
		if (initialTotalNum_m == 1) {
			int i = 0; flagNoDeletion = true;
      
			// change phase space parameters from local reference frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
			for (int j =0;j<3;j++){ 
				variable_m[j]   = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
				variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
			}

			if ( (step_m%SinglePartDumpFreq == 0) ){
				outf<<"ID"<<( itsBunch->ID[i] );
				outf<<" "<<variable_m[0]<<" "<<variable_m[3]<<" "<<variable_m[1]<<" "
					<<variable_m[4]<<" "<<variable_m[2]<<" "<<variable_m[5]<<endl;
			}
		
			double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]);//[ -pi ~ pi ]
			if((oldReferenceTheta < initialReferenceTheta-deltaTheta) && 
			   (temp_meanTheta >= initialReferenceTheta-deltaTheta )){
				*gmsg<< turnnumber <<" turn finished!"<<endl;
				turnnumber++; 
				dumpEachTurn = true;
				outfThetaEachTurn_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfThetaEachTurn_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
								 <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
								 <<" "<<temp_meanTheta/pi*180
								 <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
								 <<" "<<variable_m[2]
								 <<" "<<variable_m[5]<<endl;
			}
			// FixMe: should be defined elesewhere ! 
			// define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
			const double azimuth_angle0=0.0;
			const double azimuth_angle1=22.5/180.0*pi;
			const double azimuth_angle2=45.0/180.0*pi;
			if((oldReferenceTheta < azimuth_angle0-deltaTheta) && (temp_meanTheta >= azimuth_angle0-deltaTheta )){
				outfTheta0_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta0_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}
		
			if((oldReferenceTheta < azimuth_angle1-deltaTheta) && (temp_meanTheta >= azimuth_angle1-deltaTheta )){
				outfTheta1_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta1_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}

			if((oldReferenceTheta < azimuth_angle2-deltaTheta) && (temp_meanTheta >= azimuth_angle2-deltaTheta )){
				outfTheta2_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta2_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}
			oldReferenceTheta = temp_meanTheta;
		}


		// check whether one turn over for multi-bunch tracking.       
		if (doDumpAfterEachTurn && initialTotalNum_m > 2) {
			Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
			Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);
			for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				for (int jj=0; jj<3; jj++){
					temp_meanR(jj) += itsBunch->R[ii](jj);
				}
			}

			reduce(temp_meanR,temp_meanR,OpAddAssign());
			temp_meanR /= Vector_t(itsBunch->getTotalNum()); 
     
			// in global Cartesian frame, calculate the location in global frame of bunch
			double temp_meanTheta = calculateAngle2(temp_meanR(0), temp_meanR(1));
            
			// both for single bunch and multi-bunch 
			// avoid dump at the first step
			// dumpEachTurn has not been changed in first push
			if ( step_m > 10 )
				if((oldReferenceTheta < initialReferenceTheta-deltaTheta) && 
				   (temp_meanTheta >= initialReferenceTheta-deltaTheta )){
					*gmsg<< "one turn finished!"<<endl;
					turnnumber++;
					*gmsg<<"turn "<< turnnumber <<", total particle number = "<<itsBunch->getTotalNum()<<endl;
					dumpEachTurn = true;
				}
			oldReferenceTheta = temp_meanTheta;
		}    
		// reset Bin ID for each particle
		if ( (itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m%resetBinFreq == 0 )
			itsBunch->resetPartBinID2(eta);    

		// dump phase space distribution of bunch
		if ( ( ( (step_m+1) % Options::psDumpFreq == 0 ) && initialTotalNum_m != 2) || 
			 (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2) ) {
			  
			IpplTimings::startTimer(DumpTimer_m);
	
			itsBunch->setSteptoLastInj(SteptoLastInj);
	
			itsBunch->setTrackStep((step_m+1));

			extE_m = Vector_t(0.0, 0.0, 0.0);
			extB_m = Vector_t(0.0, 0.0, 0.0);

			//--------------------- calculate mean coordinates  of bunch -------------------------------//
			//------------  and calculate the external field at the mass of bunch-----------------------//

			Vector_t meanR=Vector_t(0.0,0.0,0.0);
			Vector_t meanP=Vector_t(0.0,0.0,0.0);

			for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				for (int j=0; j<3; j++){
					meanR(j) += itsBunch->R[ii](j);
					meanP(j) += itsBunch->P[ii](j);
				}
			}

			reduce(meanR, meanR,OpAddAssign());
			meanR /= Vector_t(( itsBunch->getTotalNum()));

			reduce(meanP, meanP,OpAddAssign());
			meanP /= Vector_t(( itsBunch->getTotalNum()));

			// FixMe: does not make sense to me
			double temp_Phi = calculateAngle(meanP(0), meanP(1));
			// the angle from positive radial direction to propagational  direction
			//                                   ==>
			// the angle from positive azimuthal direction to propagational direction
			temp_Phi -= pi/2.0;

			double temp_cosPhi = cos(temp_Phi);
			double temp_sinPhi = sin(temp_Phi);
			  
			double temp_Theta = calculateAngle(meanR(0), meanR(1));
			double temp_cosTheta = cos(temp_Theta);
			double temp_sinTheta = sin(temp_Theta);

			*gmsg<<"meanR=( "<<meanR(0)<<" "<<meanR(1)<<" "<<meanR(2)<<" ) [mm] "<<endl;

			double meanRadius = sqrt( meanR(0)*meanR(0) +  meanR(1)*meanR(1) );
			  
			beamline_list::iterator DumpSindex = FieldDimensions.begin();

			if ( ((((*DumpSindex)->second).first)[0] <= meanRadius ) &&  
				 ((((*DumpSindex)->second).first)[1] >= meanRadius ) ) {
				( ((*DumpSindex)->second).second )->apply(meanR, Vector_t(0.0), t, extE_m, extB_m);
			}
		    
			FDext_m[0] = extB_m/10.0; // kgauss -> T
			FDext_m[1] = extE_m;

			//----------------------------dump in global frame-------------------------------------//
			// Note: Don't dump when
			// 1. after one turn
			// 2.for multi-bunch mode, just after backup.
			// in order to sychronize the dump step for multi-bunch and single bunch for compare 
			// with each other during post-process phase.
			if ( !(Options::psDumpLocalFrame) ) {
				if(!OneStepAfterBackUp ){
					itsBunch->R /= Vector_t(1000.0); // mm --> m
					lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m );
					itsBunch->R *= Vector_t(1000.0); // m --> mm
					  
					*gmsg <<" Phase space dumping in global frame after step "<<step_m+1<<endl;
					*gmsg <<" Dump step = "<<lastDumpedStep_m<<", time = " << t << " [ns]" << endl;
				} 
				else {
					// only omit the dump step just after backup
					*gmsg << "Omit this Dump step = "<<lastDumpedStep_m<<endl;
					OneStepAfterBackUp = false;
				}
				//----------------------------dump in local frame-------------------------------------//
			}
			else {
				itsBunch->R -= meanR;
				itsBunch->P -= meanP;

				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_RHorizontal     =   itsBunch->R[ii](0)*temp_cosTheta  + itsBunch->R[ii](1)*temp_sinTheta;
					double temp_RLongitudinal   =  -itsBunch->R[ii](0)*temp_sinTheta  + itsBunch->R[ii](1)*temp_cosTheta;

					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;

					double temp_PHorizontal     =   itsBunch->P[ii](0)*temp_cosTheta  + itsBunch->P[ii](1)*temp_sinTheta;
					double temp_PLongitudinal   =  -itsBunch->P[ii](0)*temp_sinTheta  + itsBunch->P[ii](1)*temp_cosTheta;

					itsBunch->P[ii](0) = temp_PHorizontal;
					itsBunch->P[ii](1) = temp_PLongitudinal;
				}

				itsBunch->R /= Vector_t(1000.0); // mm --> m
				lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m );
				itsBunch->R *= Vector_t(1000.0); // m --> mm
	
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					// FixMe: make this const
					double temp_Rx   =  -itsBunch->R[ii](1)*temp_sinTheta  + itsBunch->R[ii](0)*temp_cosTheta;
					double temp_Ry   =   itsBunch->R[ii](1)*temp_cosTheta  + itsBunch->R[ii](0)*temp_sinTheta;

					itsBunch->R[ii](0) = temp_Rx;
					itsBunch->R[ii](1) = temp_Ry;

					double temp_Px   =  -itsBunch->P[ii](1)*temp_sinTheta  + itsBunch->P[ii](0)*temp_cosTheta;
					double temp_Py   =   itsBunch->P[ii](1)*temp_cosTheta  + itsBunch->P[ii](0)*temp_sinTheta;
					  
					itsBunch->P[ii](0) = temp_Px;
					itsBunch->P[ii](1) = temp_Py;
				}

				itsBunch->R += meanR;
				itsBunch->P += meanP;
				  
				*gmsg <<" Phase space dumping in local frame after step "<<step_m+1<<endl;
				*gmsg <<" Dump step = "<<lastDumpedStep_m<<", time = " << t << " [ns]" << endl;
			}
			IpplTimings::stopTimer(DumpTimer_m);
		}
	}

	for(int ii=0; ii<( itsBunch->getLocalNum()); ii++) {
		if (itsBunch->ID[ii] == 0) {
			// FixMe: FinalMomentum2  = dot(itsBunch->P[ii],itsBunch->P[ii]);
			double FinalMomentum2  = pow(itsBunch->P[ii](0),2.0) + 
				pow(itsBunch->P[ii](1),2.0) + 
				pow(itsBunch->P[ii](2),2.0);
			double FinalEnergy = (sqrt(1.0 +FinalMomentum2)-1.0)*itsReference.getM()*1.0e-6;
			*gmsgAll<<"Final energy of reference particle = "<< FinalEnergy<<" [MeV]"<<endl;
			*gmsgAll<<"Total phase space dump number(includes the initial distribution) = "<< lastDumpedStep_m+1 <<endl;
			*gmsgAll<<"One can restart simulation from the last dump step ( -restart "<< lastDumpedStep_m <<" )"<<endl;
		}
	}
	
	Ippl::Comm->barrier();
}
  
void ParallelCyclotronTracker::Tracker_RK4(){

    Inform *gmsgAll;
    gmsgAll = new  Inform("CycTracker RK4",INFORM_ALL_NODES);
	
    // prepare the tracker for doing its job
    // prepare the elements of the line to be tracked, i.e. load fieldmaps and so on

    // time steps interval between bunches for multi-bunch simulation.
    const int stepsPerTurn = itsBunch->getStepsPerTurn();

    // record how many bunches has already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch->getNumBunch();

    // decide how many energy bins. ONLY FOR MPM
    // For the time being, we set bin number equal to bunch number.
    BinCount_m = BunchCount_m;

    // determine how many sample particles used for tuning with SC
    const int NoSamplePart = 500;
    // determine when to record tune data;  
    const int StartTurnnumber = 30;
  
    // temporal 3 real space variables of particle ID=0~NoSamplePart, [x,y,z]. for tune with SC.  Unit: mm
    // todo: only processor 0 need this matrix. 
    double variable_tuneSamParts_X[NoSamplePart];
    double variable_tuneSamParts_Y[NoSamplePart];
    double variable_tuneSamParts_Z[NoSamplePart];

    beamline_list::iterator sindex = FieldDimensions.begin();
    const double harm = (((*sindex)->second).second )-> getCyclHarm();

    // load time
    double t  = itsBunch->getT()*1.0e9;
    const double dt = itsBunch->getdT()*1.0e9*harm; //[s]-->[ns]

    // find the injection time interval
	if (numBunch_m > 1){
		const double RfFreq = (((*sindex)->second).second )-> getRfFrequ(); // [MHz]
		*gmsg << "Time interval between neighbour bunches is set to "<<stepsPerTurn*dt<<"[ns]"<<endl;
	}

	/*
	  transform initial coordinates and momenta from local beam frame (relative value to the bunchcenter )
	  to global Cartesian frame (absolute value ).
	  Do some other initialisation work.
	  for multi-particle Mode and single Particle Mode, output particles of ID = 0 and 1 for each  dumpfreq steps
	*/
	
	string SfileName = OPAL.getInputFn();
	int pdot = SfileName.find(string("."),0);
	SfileName.erase(pdot,SfileName.size()-pdot);
	
	string  SfileName1 = SfileName+string("-trackOrbit.dat");
	ofstream outf;
	outf.setf(ios::scientific, ios::floatfield );
	outf.precision(8);
	
	if (initialTotalNum_m > 2)
		if ( myNode_m == 0 &&  flagDoTune_m ){
			ifstream inf;
			char skipChar[100];
			double skipNum;
			Vector_t tempSEO;
			string  SfileName5 = string("S")+SfileName+string("-trackOrbit.dat");
			inf.open(SfileName5.c_str(),ios::in);
			if(inf.fail()) {
				*gmsg<<"Cannot open file "<<SfileName5<<" for tuning calculation, please check if it really exists."<<endl;
				exit(1);
			}
			if(!inf.eof()){
				for(int i=0; i<14; i++) {
					inf>>skipChar;
				}
			}
			
			while(!inf.eof()){
				inf>>skipChar>> tempSEO(1) >> skipNum >> tempSEO(2) >> skipNum >> skipNum >> skipNum;
				if(inf.eof()) break;
				tempSEO(0) = calculateAngle(tempSEO(1),tempSEO(2));
				variable_SEO_m.push_back(tempSEO);
			}
			inf.close();
			*gmsg<<"Finish reading in SEO file for tuning calculation"<<endl;
		}

	// get data from h5 file for restart run
	if(OPAL.inRestartRun()){
		restartStep0_m = itsBunch->getTrackStep();
		step_m = restartStep0_m;
		*gmsg <<"Restart at integration step "<<restartStep0_m<<endl;
	}
	
	if (OPAL.hasBunchAllocated() && Options::scan ) {
		lastDumpedStep_m = 0;
		t = 0.0;
	}

	*gmsg<< "* Beginning of this run is at t= "<< t << " [ns]"<<endl;
	*gmsg<< "* The time step is set to dt= "<<dt << " [ns]"<<endl;

	// add header in the data dump file for particle ID equal 0 and 1
	if(myNode_m == 0) {
		if(OPAL.inRestartRun()){
			outf.open(SfileName1.c_str(),ios::app);
			outf      << "# Restart at "<<step_m<<" dumping step"<<endl;
		}
		else{
			outf.open(SfileName1.c_str());
			outf      <<"# ID   x [mm]          px [rad]       y [mm]          py [rad]        z [mm]          pz [rad]"
					  << endl;
		}
	}

	// for single Particle Mode, output at zero degree.
	if (initialTotalNum_m == 1)
		openFiles(SfileName);

	if(!OPAL.inRestartRun()) {
		PathLength_m = 0.0;
		// Force the initial phase space values of the particle with ID=0 to zero, to set it as a reference particle.
		if (initialTotalNum_m > 2) { // only for mulit-particle mode
			for(int ii=0; ii<initialLocalNum_m; ii++) {
				if (itsBunch->ID[ii] == 0) {
					itsBunch->R[ii] = Vector_t(0.0);
					itsBunch->P[ii] = Vector_t(0.0);
				}
				// for tuning calculation
				if(flagDoTune_m){
					if (itsBunch->ID[ii] == 1) {
						itsBunch->R[ii] = Vector_t(0.003, 0.0, 0.003);
						itsBunch->P[ii] = Vector_t(0.0);
					}
				}
			}
		}

		if ( Options::psDumpLocalFrame ) {
			// dump the initial distribution
			lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
			*gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the local frame "<<endl;
		}

		itsBunch->R *= Vector_t(1000.0); // m --> mm
		for (int ii = 0; ii < initialLocalNum_m; ii++) {
			// rotate the local beam frame for angle of minus theta
			double tempDx = itsBunch->R[ii](0) * cosRefTheta - itsBunch->R[ii](1) * sinRefTheta;
			double tempDy = itsBunch->R[ii](0) * sinRefTheta + itsBunch->R[ii](1) * cosRefTheta;

			// move the origin point of global frame
			itsBunch->R[ii](0) = tempDx + referenceR * cosRefTheta;
			itsBunch->R[ii](1) = tempDy + referenceR * sinRefTheta;

			double tempPx = ( itsBunch->P[ii](0) + referencePr) * cosRefTheta 
				-  ( itsBunch->P[ii](1) + referencePt) * sinRefTheta;
			double tempPy = ( itsBunch->P[ii](0) + referencePr) * sinRefTheta 
				+  ( itsBunch->P[ii](1) + referencePt) * cosRefTheta;

			itsBunch->P[ii](0)=tempPx;
			itsBunch->P[ii](1)=tempPy;
			// initialize the bin number of the first bunch to 0
			itsBunch->Bin[ii]=0;
		}
		if ( !(Options::psDumpLocalFrame) ) {
			// dump the initial distribution
			itsBunch->R /= Vector_t(1000.0); // mm --> m
			lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
			itsBunch->R *= Vector_t(1000.0); // m --> mm
			*gmsg << "* meanR=( "<<referenceR * cosRefTheta<<" "<<referenceR * sinRefTheta<<" "<<0.0<<" ) [mm] "<<endl;
			*gmsg << "* Phase space dump at step " << lastDumpedStep_m << " in the global frame "<<endl;
		}
		
		// AUTO mode
		if (multiBunchMode_m == 2 ) {
			Vector_t Rmean = itsBunch->get_centroid()*1000.0; // m->mm
			RLastTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );
			RThisTurn_m = RLastTurn_m;	
			*gmsg<<"* Initial Radial position = "<<RThisTurn_m<<" [mm]"<<endl;
		}

		/*
		  FixMe: more generic initial azimuthal position, not finished yet
		  double referencePx = referencePr*cosRefTheta - referencePt*sinRefTheta;
		  double referencePy = referencePr*sinRefTheta + referencePt*cosRefTheta;
		
		  double referencePhi = calculateAngle(referencePx,referencePy);
		
		  referencePhi -= pi/2.0;
		*/

		if ((initialTotalNum_m > 2 ) && ( numBunch_m > 1) && (multiBunchMode_m == 1) ){
			// backup initial distribution
			assert( initialR_m = new Vector_t[initialLocalNum_m]);
			assert( initialP_m = new Vector_t[initialLocalNum_m]);
			
			for (int ii=0; ii<initialLocalNum_m; ii++){
				initialR_m[ii] = itsBunch->R[ii];
				initialP_m[ii] = itsBunch->P[ii];
			}
		}
	}
	else{
		PathLength_m = itsBunch->getLPath();
		// AUTO mode
		if (multiBunchMode_m == 2 ) {
			itsBunch->calcBeamParameters_cycl();
			Vector_t Rmean = itsBunch->get_centroid()*1000.0; // m->mm
			RLastTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );
			RThisTurn_m = RLastTurn_m;		
			*gmsg<<"* Radial position at restart position = "<<RThisTurn_m<<" [mm]"<<endl;
		}
		itsBunch->R *= Vector_t(1000.0); // m --> mm
	}
	
	//  read in some control parameters
	const int SinglePartDumpFreq = Options::sptDumpFreq;
	const int RepartitionFreq = Options::repartFreq;
	const int resetBinFreq = Options::rebinFreq;
	const int scSolveFreq = Options::scSolveFreq;
	const bool doDumpAfterEachTurn =Options::psDumpEachTurn;

	int boundpDestroyFreq = 10; // todo: is it better treat as a control parameter

	// prepare for dump after each turn
	const double initialReferenceTheta = referenceTheta/180.0*pi;
	double oldReferenceTheta = initialReferenceTheta;

	*gmsg<<"* Single particle trajectory dump frequency is set to "<<SinglePartDumpFreq<<endl;
	*gmsg<<"* The frequency to solve space charge fields is set to "<<scSolveFreq<<endl;
	*gmsg<<"* The repartition frequency is set to "<<RepartitionFreq<<endl;

	if(numBunch_m > 1)
		*gmsg<<"* The particles energy bin reset frequency is set to "<<resetBinFreq<<endl;


	// if initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
	vector<double> Ttime,Tdeltr,Tdeltz;
	// if flagDoTune == tune for multi-particles mode and prepare for transverse tuning calculation
	vector<double> TdeltrSC[NoSamplePart],TdeltzSC[NoSamplePart];
	vector<int> TturnNumber;
	int turnnumber = 1;
	int lastTurn = 1;
    
	bool flagNoDeletion;

	// flag to determine when to transit from single-bunch to multi-bunches mode
	bool flagTransition = false;
	// flag to omit the dump step after backup dump step for multi-bunches mode
	bool OneStepAfterBackUp = false;
	// step point determining the next time point of check for transition
	int stepsNextCheck = step_m + itsBunch->getStepsPerTurn();
    
	const double deltaTheta = pi/(stepsPerTurn); // half of the average angle per step
	// record at which angle the space charge are solved
	double angleSpaceChargeSolve;
	// local particle number
	int nlp;
    
	if(initialTotalNum_m == 1){
		*gmsg << "* ---------------------------- SINGLE PARTICLE MODE------ ----------------------------*** "<<endl;
		*gmsg << "* Instruction: when the total particle number equal to 1, single particle mode is triggered automatically," <<endl
			  << "* The initial distribution file must be specified which should contain only one line for the single particle "<<endl
			  << "* ------------WARNING: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** "<<endl;
		if (Ippl::getNodes() != 1)
			throw OpalException("Error in ParallelCyclotronTracker::execute","SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
		
	}
	else if (initialTotalNum_m == 2){
		*gmsg << "* ------------------------ STATIC EQUILIBRIUM ORBIT MODE ----------------------------*** "<<endl;
		*gmsg << "* Instruction: when the total particle number equal to 2, SEO mode is triggered automatically." <<endl
			  << "* This mode does NOT include any RF cavities. The initial distribution file must be specified"<<endl
			  << "* In the file the first line is for reference particle and the second line is for offcenter particle."<<endl
			  << "* The tuning is calculated by FFT routines based on these two particles. "<<endl
			  << "* ------------WARNING: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------*** "<<endl;
		if (Ippl::getNodes() != 1)
			throw OpalException("Error in ParallelCyclotronTracker::execute","SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
	}
	else{
		*gmsg << "* ---------------------------- MULTI-PARTICLES MODE------ ----------------------------*** "<<endl;
		nlp=itsBunch->getLocalNum();

		/**
		   FixMe: nedd clarification on that
		   Need itsBunch->boundpNoRep() otherwise dead lock!
			 
		   itsBunch->boundpNoRep();
		*/
		
	
		nlp=itsBunch->getLocalNum();
		checkNumPart(string("Before boundp "), nlp);

		itsBunch->boundp();
		nlp=itsBunch->getLocalNum();
		
		checkNumPart(string("After boundp"),nlp);
		  
	}
	
	// main integration loop
	*gmsg << "***---------------------------- Start tracking   ------------------------------------------*** "<<endl;
	for(step_m; step_m < maxSteps_m; step_m++) {
		bool dumpEachTurn = false;
		bool flagNeedUpdate = false;
		
		/* 
		   execute one step for initialTotalNum_m > 2 mode , initialTotalNum_m = 2 mode or initialTotalNum_m =1 mode
		   
		*/
		if (initialTotalNum_m > 2) {
			// singel particle dumping
			if ( step_m%SinglePartDumpFreq == 0 ){ // dump
				IpplTimings::startTimer(DumpTimer_m);

				double x;
				int  id;
				vector<double> tmpr;
				vector<int> tmpi;
				
				int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
				
				if (!flagDoTune_m) {
					// for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
					int found[2] = {-1,-1};
					int counter = 0;

					for(int ii=0; ii<( itsBunch->getLocalNum()); ii++) {
						if (itsBunch->ID[ii] == 0) {
							found[counter] = ii;
							counter++;
						}
						if (itsBunch->ID[ii] == 1) {
							found[counter] = ii;
							counter++;
						}
					}
					// for the regular modes only the space data of particles with ID = 0 and 1 need be transfored 
					if(myNode_m == 0) {
						// for root node
						int notReceived =  Ippl::getNodes() - 1;
						int numberOfPart =0;

						while (notReceived > 0) {
							int node = COMM_ANY_NODE;
							Message* rmsg =  Ippl::Comm->receive_block(node, tag);
							if (rmsg == 0)
								ERRORMSG("Could not receive from client nodes in main." << endl);
							notReceived--;
							rmsg->get(&numberOfPart);
							for (int ii=0; ii < numberOfPart; ii++) {
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
						for (int ii=0; ii < counter; ii++) {
							tmpi.push_back(itsBunch->ID[found[ii]]);
							for (int jj=0;jj<3;jj++) {
								tmpr.push_back(itsBunch->R[found[ii]](jj));
								tmpr.push_back(itsBunch->P[found[ii]](jj));
							}
						}
						vector<double>::iterator itParameter = tmpr.begin();
						vector<int>::iterator  itId = tmpi.begin();

						for ( itId = tmpi.begin(); itId != tmpi.end(); itId++ ){
							outf<<"ID"<<*itId;
							for (int ii = 0;ii < 6; ii++ ){
								outf<<" "<<*itParameter;
								itParameter++;
							}
							outf <<endl;
						}
						// sample frequency = SinglePartDumpFreq
					}else {
						// for other nodes
						Message* smsg = new Message();
						smsg->put(counter);
						for (int ii=0; ii < counter; ii++) {
							smsg->put(itsBunch->ID[found[ii]]);
							for (int jj=0;jj<3;jj++) {
								smsg->put(itsBunch->R[found[ii]](jj));
								smsg->put(itsBunch->P[found[ii]](jj));
							}
						}
						bool res = Ippl::Comm->send(smsg, 0, tag);
						if (!res)
							ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
						// FixMe: why it block at here if I delete smsg? delete smsg;
						 
					}
				}
				else {
					if ( turnnumber <= StartTurnnumber )
						StartTime_m = t ; // [ns]
					else {
						// for all nodes, find the location of particle with ID = 0 & 1~NoSamplePart in bunch containers
						int found[NoSamplePart];
						int counter = 0;	 
						int counterAll = 0;
						for(int ii=0; ii<NoSamplePart; ii++) 
							found[ii] = -1;
						for(int ii=0; ii<( itsBunch->getLocalNum()); ii++) {
							if (itsBunch->ID[ii] < NoSamplePart){
								found[counter] = ii;
								counter++;
							}
						}
						reduce(counter, counterAll, OpAddAssign());
						*gmsg<<"Found "<< counterAll<<" sample particles."<<endl;
	   	  	    
						if(myNode_m == 0) {
							// for root node
							int notReceived =  Ippl::getNodes()-1;
							int numberOfPart = counter;
							
							while (notReceived > 0) {
								int node = COMM_ANY_NODE;
								Message* rmsg =  Ippl::Comm->receive_block(node, tag);
								if (rmsg == 0)
									ERRORMSG("Could not receive from client nodes in main." << endl);
								notReceived--;
								rmsg->get(&numberOfPart);
								for (int ii=0; ii < numberOfPart; ii++) {
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
	   
							// check particles at local memory 
							for (int ii=0; ii < counter; ii++) {
								tmpi.push_back(itsBunch->ID[found[ii]]);
								for (int jj=0;jj<3;jj++) {
									tmpr.push_back(itsBunch->R[found[ii]](jj));
									tmpr.push_back(itsBunch->P[found[ii]](jj));
								}
							}
							vector<double>::iterator itParameter = tmpr.begin();
							vector<int>::iterator  itId = tmpi.begin();
	    
							// dump into file 'INPUT_FILE_NAME-trackorbit.dat' 
							for ( itId = tmpi.begin(); itId != tmpi.end(); itId++ ){
								if ( *itId < 10){
									outf<<"ID"<<*itId;
									for (int ii = 0;ii < 6; ii++ ){
										outf<<" "<<*itParameter;
										itParameter++;
									}
									outf <<endl;
								}
							}
							*gmsg<<"Finish receiving!"<<endl;

							double r_tuning[2], z_tuning[2] ;
							itParameter = tmpr.begin();
							itId = tmpi.begin();

							for ( itId = tmpi.begin(); itId != tmpi.end(); itId++ ){
								variable_tuneSamParts_X[*itId]=*itParameter;
								itParameter = itParameter+2;
								variable_tuneSamParts_Y[*itId]=*itParameter;
								itParameter = itParameter+2;
								variable_tuneSamParts_Z[*itId]=*itParameter;
								itParameter = itParameter+2;
							}
							
							for (int ii=0; ii<NoSamplePart; ii++ ){
								double tempTheta = 0.0;
								double minDTheta = 100.0;
								vector<Vector_t>::iterator itSEO, itSEO2;
								double tempx, tempy;
								double tempDTheta;
								variable_tune0_m[0]=variable_tuneSamParts_X[ii];
								variable_tune0_m[1]=variable_tuneSamParts_Y[ii];
								variable_tune0_m[2]=variable_tuneSamParts_Z[ii];
								tempTheta = calculateAngle(variable_tune0_m[0], variable_tune0_m[1]);
								
								for ( itSEO = variable_SEO_m.begin(); itSEO != variable_SEO_m.end(); itSEO++ ){
									tempDTheta = abs(tempTheta-((*itSEO)(0)));
									if (tempDTheta < minDTheta){
										minDTheta = tempDTheta;
										itSEO2 = itSEO;
									}
								}
								
								tempx = (*itSEO2)(1);
								tempy = (*itSEO2)(2);
								double angle = calculateAngle(tempx, tempy);

								r_tuning[0] = tempx*cos(tempTheta)+tempy*sin(tempTheta);
								
								r_tuning[1] = variable_tuneSamParts_X[ii]*cos(tempTheta)+variable_tuneSamParts_Y[ii]*sin(tempTheta);
								z_tuning[1] = variable_tuneSamParts_Z[ii];
								
								TdeltrSC[ii].push_back(r_tuning[1]  - r_tuning[0]); 
								TdeltzSC[ii].push_back(z_tuning[1]);
							}
							
							Ttime.push_back((t-StartTime_m)*1.0e-9);
							TturnNumber.push_back(turnnumber-StartTurnnumber);
							*gmsg<<"Time = "<<t<<", turn number = "<<turnnumber<<endl;
							*gmsg<<"Acount time = "<<t-StartTime_m <<", turn number = "<< turnnumber-StartTurnnumber<<endl; 
						}
						else {
							// for other nodes
							Message* smsg = new Message();
							smsg->put(counter);
							for (int ii=0; ii < counter; ii++) {
								smsg->put(itsBunch->ID[found[ii]]);
								for (int jj=0;jj<3;jj++) {
									smsg->put(itsBunch->R[found[ii]](jj));
									smsg->put(itsBunch->P[found[ii]](jj));
								}
							}
							bool res = Ippl::Comm->send(smsg, 0, tag);
							if (!res)
								ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
						}
					}
				}
				IpplTimings::stopTimer(DumpTimer_m);
			}
			//end dump
			Ippl::Comm->barrier();
		    
			// bunch injection
			if (numBunch_m > 1 ){
				size_t injectLocalNUM =initialTotalNum_m;
				if( (BunchCount_m == 1) && (multiBunchMode_m == 2) && (!flagTransition)) {
					if (step_m == stepsNextCheck ){
						// under 3 conditions, following code will be execute
						// to check the distance between two neighborring bunches
						// 1.multi-bunch mode, AUTO sub-mode
						// 2.After each revolution
						// 3.only one bunch exists
						
						*gmsg<<"checking for automatically injecting new bunch ..." << endl;
						
						itsBunch->R /= Vector_t(1000.0); // mm --> m
						itsBunch->calcBeamParameters_cycl();
						itsBunch->R *= Vector_t(1000.0); // m --> mm

						Vector_t Rmean = itsBunch->get_centroid() * 1000.0; // m --> mm
						
						RThisTurn_m = sqrt( pow(Rmean[0], 2.0) + pow(Rmean[1], 2.0) );
						
						Vector_t Rrms = itsBunch->get_rrms() * 1000.0; // m --> mm

						double XYrms =  sqrt( pow(Rrms[0], 2.0) + pow(Rrms[1], 2.0) );

						// if the distance between two nieghbour bunch is less than 5 times of its 2D rms size
						// start multi-bunch simulation, fill current phase space to initialR and initialP arrays
						
						if ((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
							
							// since next turn, start multi-bunches

							itsBunch->R /= Vector_t(1000.0); // mm --> m
							lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m);
							itsBunch->R *= Vector_t(1000.0); // m --> mm
													  
							backupDumpStep_m = lastDumpedStep_m;
							flagTransition = true;
							OneStepAfterBackUp =true;

							*gmsg << "In order to inject a new bunch after one revolution, Dump step = "
								  << lastDumpedStep_m <<", time ="<<t<<" [ns]"<<endl;
							*gmsg<<"**************** Another revolution later, Multi-Bunch Mode will be invorked automatically ****************"<<endl;

						}
						stepsNextCheck += stepsPerTurn;
						
						*gmsg<<"RLastTurn = "<<RLastTurn_m<<" [mm]"<<endl;
						*gmsg<<"RThisTurn = "<<RThisTurn_m<<" [mm]"<<endl;
						*gmsg<<"    XYrms = "<<XYrms    <<" [mm]"<<endl;
						
						RLastTurn_m = RThisTurn_m;
					}

				}
				else if ( (BunchCount_m < numBunch_m) && (step_m == stepsNextCheck) ) {
									  
					// under 4 conditions, following code will be execute
					// to read new bunch from hdf5 format file for FORCE or AUTO mode
					// 1.multi-bunch mode
					// 2.after each revolution
					// 3.existing bunches is less than the specified bunches
					// 4.FORCE mode, or AUTO mode with flagTransition = true
					// Note: restart from 1 < BunchCount < numBunch_m must be avoided.
					*gmsg<<"time = "<<t<<" [ns], Now inject a new bunch by reading data from h5 file ... ... ..."<<endl;
					
					BunchCount_m++;

					// read initial distribution from h5 file
					if(multiBunchMode_m == 1 )
						readOneBunch(BunchCount_m-1, 0);
					else if (multiBunchMode_m == 2 )
						readOneBunch(BunchCount_m-1, backupDumpStep_m);
					
					itsBunch->setNumBunch(BunchCount_m);

					stepsNextCheck += stepsPerTurn;

					// update  after injection
					itsBunch->boundp();

					Ippl::Comm->barrier();
					*gmsg<<"read initial distribution from h5 file"<<endl;
					*gmsg<<BunchCount_m<<"th bunch injected, total particle number = "<<itsBunch->getTotalNum()<<endl;

				}
				else if ( BunchCount_m == numBunch_m ) {
					// After this, numBunch_m is wrong but useless
					numBunch_m--;
				}
			}
		    
			// Calculate SC field before each time step and keep constant during integration.
			// Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
			if(itsBunch->hasFieldSolver() && initialTotalNum_m >=1000) {
				if(step_m % scSolveFreq == 0){
					// Firstly reset E and B to zero before fill new space charge field data for each track step
					itsBunch->Bf = Vector_t(0.0);
					itsBunch->Ef = Vector_t(0.0);
			    
					IpplTimings::startTimer(TransformTimer_m);
			    
					//HERE transform particles coordinates to local frame (rotate and shift)
					Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
					Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);

					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						for (int jj=0; jj<3; jj++){
							temp_meanR(jj) += itsBunch->R[ii](jj);
						}
					}

					// sum all particles in the same bin using reduce() for parallel enviroment,
					// and change N to initialTotalNum_m.
					reduce(temp_meanR,temp_meanR,OpAddAssign());
					temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX

					// in global Cartesian frame, calculate the location in global frame of bunch
					double temp_meanTheta = calculateAngle2(temp_meanR(0), temp_meanR(1));

					// both for single bunch and multi-bunch 
					if ( step_m > 10 ) // avoid dump at the first step
						if((oldReferenceTheta < initialReferenceTheta-deltaTheta) && 
						   (temp_meanTheta >= initialReferenceTheta-deltaTheta )){

							*gmsg<< "one turn finished!"<<endl;
							turnnumber++;
							dumpEachTurn = true;
						}

					oldReferenceTheta = temp_meanTheta;

					if ( (itsBunch->weHaveBins()) && BunchCount_m > 1 ) {
						double temp_binsMeanPhi = itsBunch->calcMeanPhi();

						// the angle from positive radial direction to propagational  direction
						//                                   ==>
						// the angle from positive azimuthal direction to propagational direction
						temp_binsMeanPhi -= pi/2.0;
						angleSpaceChargeSolve = temp_binsMeanPhi;
				    
						double cosTemp_binsMeanPhi = cos(temp_binsMeanPhi);
						double sinTemp_binsMeanPhi = sin(temp_binsMeanPhi);


						// remove mean coordinates
						itsBunch->R -= temp_meanR;
						
						//scale coordinates
						itsBunch->R /= Vector_t(1000.0); // mm --> m

						// rotate from global frame to local frame(transverse horizontal,longitudinal,transverse vertical)
						// For multi-bin, rotate the frame for temp_binsMeanPhi degree
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_binsMeanPhi  
								+ itsBunch->R[ii](1)*sinTemp_binsMeanPhi;
							double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_binsMeanPhi  
								+ itsBunch->R[ii](1)*cosTemp_binsMeanPhi;

							itsBunch->R[ii](0) = temp_RHorizontal;
							itsBunch->R[ii](1) = temp_RLongitudinal;
						}

						if ( (step_m+1) % boundpDestroyFreq == 0 )
							itsBunch->boundp_destroy();
						else
							itsBunch->boundp();
						
						IpplTimings::stopTimer(TransformTimer_m);
						
						// FixMe: check this
						//   IpplTimings::startTimer(BinRepartTimer_m);
						
						// todo: overuse binaryRepart() can cause buffer overflow during large scale job.
						// as a temporal solution, reduce its using frequency can solve such kind of problem.
						//   if( (step_m+1)%RepartitionFreq == 0 )
						// 	  {
						// 	      NDIndex<3> ldom = itsBunch->getFieldLayout().getLocalNDIndex();
						// 	      double hx, hy, hz;
						// 	      int IndexMaxX,IndexMaxY,IndexMaxZ;
						// 	      int IndexMinX,IndexMinY,IndexMinZ;
						// 	      hx = itsBunch->get_hr()(0);
						// 	      hy = itsBunch->get_hr()(1);
						// 	      hz = itsBunch->get_hr()(2);
						// 	      IndexMaxX =  ldom[0].max();
						// 	      IndexMinX =  ldom[0].min();
						// 	      IndexMaxY =  ldom[1].max();
						// 	      IndexMinY =  ldom[1].min();
						// 	      IndexMaxZ =  ldom[2].max();
						// 	      IndexMinZ =  ldom[2].min();
						// 	      *gmsg<<"hx = "<<hx<<", hy = "<<hy<<", hz = "<<hz<<endl;
						// 	      Ippl::Comm->barrier();
						// 	      *gmsgAll <<"min/max X:"<< IndexMinX<<"--"<<IndexMaxX
					    //                 << ",Y: "<< IndexMinY<<"--"<<IndexMaxY
						//                 << ", Z: "<< IndexMinZ<<"--"<<IndexMaxZ<<endl ;
						// 	      Ippl::Comm->barrier();
						// 	      *gmsgAll <<"After boundp, particle  "<<itsBunch->getLocalNum()<<endl ;
						// 	      Ippl::Comm->barrier();
						// 	  }
						//           IpplTimings::stopTimer(BinRepartTimer_m);
						//           Ippl::Comm->barrier();
						
						
						// calcualte gamma for each energy bin
						itsBunch->calcGammas_cycl();
											  
						if ((step_m%RepartitionFreq) == 0 ) {
							nlp=itsBunch->getLocalNum();
							checkNumPart(string("Before repartation: "), nlp);
							IpplTimings::startTimer(BinRepartTimer_m);
							itsBunch->do_binaryRepart();
							IpplTimings::stopTimer(BinRepartTimer_m);
							Ippl::Comm->barrier();
							nlp=itsBunch->getLocalNum();
							checkNumPart(string("After repartation: "), nlp);
						}
						// calculate space charge field for each energy bin
						for (int b=0; b<itsBunch->getLastemittedBin(); b++) {
				      
							itsBunch->setBinCharge(b,itsBunch->getChargePerParticle());
							itsBunch->computeSelfFields_cycl(b);
						}

						itsBunch->Q = itsBunch->getChargePerParticle();

						IpplTimings::startTimer(TransformTimer_m);
						
						// HERE transform particles coordinates back to global frame (rotate and shift)
						// rotate back from local reference frame to global frame
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_binsMeanPhi  
								- itsBunch->R[ii](1)*sinTemp_binsMeanPhi;
							double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_binsMeanPhi  
								+ itsBunch->R[ii](1)*cosTemp_binsMeanPhi;

							itsBunch->R[ii](0) = temp_RHorizontal;
							itsBunch->R[ii](1) = temp_RLongitudinal;
						}

						//scale coordinates back
						itsBunch->R *= Vector_t(1000.0); // m --> mm

						// retrieve mean coordinates
						itsBunch->R += temp_meanR;

						// HERE transform self field back to global frame (rotate)

						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_binsMeanPhi  
								- itsBunch->Ef[ii](1)*sinTemp_binsMeanPhi;
							double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_binsMeanPhi  
								+ itsBunch->Ef[ii](1)*cosTemp_binsMeanPhi;
							
							itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
							itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
						}

						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_binsMeanPhi  
								- itsBunch->Bf[ii](1)*sinTemp_binsMeanPhi;
							double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_binsMeanPhi  
								+ itsBunch->Bf[ii](1)*cosTemp_binsMeanPhi;

							itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
							itsBunch->Bf[ii](1) = temp_E2;  // By,T
						}
					}
					else{
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							for (int jj=0; jj<3; jj++){
								temp_meanP(jj) += itsBunch->P[ii](jj);
							}
						}

						reduce(temp_meanP,temp_meanP,OpAddAssign());
						temp_meanP /= Vector_t(itsBunch->getTotalNum());
						
						// in global Cartesian frame, calculate the direction of longitudinal angle of bunch

						double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));
				
						// the angle from positive radial direction to propagational  direction
						//                                   ==>
						// the angle from positive azimuthal direction to propagational direction
						temp_meanPhi -= pi/2.0;
						angleSpaceChargeSolve = temp_meanPhi;
				
						double cosTemp_meanPhi = cos(temp_meanPhi);
						double sinTemp_meanPhi = sin(temp_meanPhi);

						double temp_meanPLongitudinal2 = pow(temp_meanP(0),2.0) + pow(temp_meanP(1),2.0);
						double temp_meangamma = sqrt(1.0 +temp_meanPLongitudinal2);

						// remove mean coordinates
						itsBunch->R -= temp_meanR;

						//scale coordinates
						itsBunch->R /= Vector_t(1000.0); // mm --> m

						// rotate from global frame to local frame(transverse horizontal,longitudinal, transverse vertical)
						// For single bin, rotate the frame for temp_meanPhi degree
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  
								+ itsBunch->R[ii](1)*sinTemp_meanPhi;
							double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_meanPhi  
								+ itsBunch->R[ii](1)*cosTemp_meanPhi;

							itsBunch->R[ii](0) = temp_RHorizontal;
							itsBunch->R[ii](1) = temp_RLongitudinal;
						}
						
						if ( (step_m+1) % boundpDestroyFreq == 0 )
							itsBunch->boundp_destroy();
						else
							itsBunch->boundp();

						IpplTimings::stopTimer(TransformTimer_m);

						// FixMe:
						//    IpplTimings::startTimer(BinRepartTimer_m);

						// todo: overuse binaryRepart() can cause baffer overflow during large scale job.
						// as a temporal solution, reduce its using frequency can solve such kind of problem.

						//  if( (step_m+1)%RepartitionFreq == 0 )
						// 	  {
						// 	      NDIndex<3> ldom = itsBunch->getFieldLayout().getLocalNDIndex();
						// 	      double hx, hy, hz;
						// 	      int IndexMaxX,IndexMaxY,IndexMaxZ;
						// 	      int IndexMinX,IndexMinY,IndexMinZ;
						// 	      hy = itsBunch->get_hr()(1);
						// 	      hz = itsBunch->get_hr()(2);
						// 	      IndexMaxX =  ldom[0].max();
						// 	      IndexMinX =  ldom[0].min();
						// 	      IndexMaxY =  ldom[1].max();
						// 	      IndexMinY =  ldom[1].min();
						// 	      IndexMaxZ =  ldom[2].max();
						// 	      IndexMinZ =  ldom[2].min();
						// 	      *gmsg<<"hx = "<<hx<<", hy = "<<hy<<", hz = "<<hz<<endl;
						// 	      Ippl::Comm->barrier();
						// 	      *gmsgAll<<"min/max X:"<< IndexMinX<<"--"<<IndexMaxX<<", Y: "
						//                << IndexMinY<<"--"<<IndexMaxY<<", Z: "<< IndexMinZ<<"--"<<IndexMaxZ<<endl ;
						// 	      Ippl::Comm->barrier();
						// 	      *gmsgAll <<"After boundp, particle  "<<itsBunch->getLocalNum()<<endl ;
						// 	      Ippl::Comm->barrier();
						// 	      *gmsg << "***---------------------------- DO REPART   ----------------------*** "<<endl;
						// 	      itsBunch->do_binaryRepart();
						// 	      *gmsg << "***---------------------------- REPART DONE ----------------------*** "<<endl;
						// 	  }
						//        IpplTimings::stopTimer(BinRepartTimer_m);
						//           Ippl::Comm->barrier();
				
						if ((step_m%RepartitionFreq) == 0 ) {
							nlp=itsBunch->getLocalNum();
							checkNumPart(string("Before repartation: "), nlp);
							IpplTimings::startTimer(BinRepartTimer_m);
							itsBunch->do_binaryRepart();
							IpplTimings::stopTimer(BinRepartTimer_m);
							Ippl::Comm->barrier();
							nlp=itsBunch->getLocalNum();
							checkNumPart(string("After repartation: "), nlp);
						}
						itsBunch->computeSelfFields_cycl( temp_meangamma );

						IpplTimings::startTimer(TransformTimer_m);
						// HERE transform particles coordinates back to global frame (rotate and shift)
						// rotate back from local reference frame to global frame
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  
								- itsBunch->R[ii](1)*sinTemp_meanPhi;
							double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_meanPhi  
								+ itsBunch->R[ii](1)*cosTemp_meanPhi;

							itsBunch->R[ii](0) = temp_RHorizontal;
							itsBunch->R[ii](1) = temp_RLongitudinal;
						}

						//scale coordinates back
						itsBunch->R *= Vector_t(1000.0); // m --> mm
						
						// retrieve mean coordinates
						itsBunch->R += temp_meanR;

						// HERE transform self field back to global frame (rotate)
						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_meanPhi  - itsBunch->Ef[ii](1)*sinTemp_meanPhi;
							double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_meanPhi  + itsBunch->Ef[ii](1)*cosTemp_meanPhi;
							itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
							itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
						}

						for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
							double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_meanPhi  - itsBunch->Bf[ii](1)*sinTemp_meanPhi;
							double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_meanPhi  + itsBunch->Bf[ii](1)*cosTemp_meanPhi;
							itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
							itsBunch->Bf[ii](1) = temp_E2;  // By,T
						}
					}

					IpplTimings::stopTimer(TransformTimer_m);
			    
				}
				else {
					//HERE transform particles coordinates to local frame (rotate and shift)
					Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
					Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);

					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						for (int jj=0; jj<3; jj++){
							temp_meanR(jj) += itsBunch->R[ii](jj);
						}
					}
					// sum all particles in the same bin using reduce() for parallel enviroment,
					// and change N to initialTotalNum_m.
					reduce(temp_meanR,temp_meanR,OpAddAssign());
					temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
			    
					// in global Cartesian frame, calculate the location in global frame of bunch
					double temp_meanTheta = calculateAngle2(temp_meanR(0), temp_meanR(1));
			    
					// both for single bunch and multi-bunch
					if ( step_m > 10 ) // avoid dump at the first step
						if((oldReferenceTheta < initialReferenceTheta-deltaTheta) 
						   && (temp_meanTheta >= initialReferenceTheta-deltaTheta )){

							*gmsg<< "one turn finished!"<<endl;
							turnnumber++;
							dumpEachTurn = true;
						}
					
					oldReferenceTheta = temp_meanTheta;
					
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						for (int jj=0; jj<3; jj++){
							temp_meanP(jj) += itsBunch->P[ii](jj);
						}
					}

					reduce(temp_meanP,temp_meanP,OpAddAssign());
					temp_meanP /= Vector_t(itsBunch->getTotalNum());
					
					// in global Cartesian frame, calculate the direction of longitudinal angle of bunch
			    
					double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));

					// the angle from positive radial direction to propagational  direction
					//                                   ==>
					// the angle from positive azimuthal direction to propagational direction
					temp_meanPhi -= pi/2.0;
					double temp_deltaPhi = temp_meanPhi-angleSpaceChargeSolve;
			    
					double cosTemp_deltaPhi = cos(temp_deltaPhi);
					double sinTemp_deltaPhi = sin(temp_deltaPhi);
			    			    
   
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Ef[ii](0)*cosTemp_deltaPhi  - itsBunch->Ef[ii](1)*sinTemp_deltaPhi;
						double temp_E2 =  itsBunch->Ef[ii](0)*sinTemp_deltaPhi  + itsBunch->Ef[ii](1)*cosTemp_deltaPhi;
						itsBunch->Ef[ii](0) = temp_E1;  // Ex,V/m
						itsBunch->Ef[ii](1) = temp_E2;  // Ey,V/m
					}
					
					for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
						double temp_E1 =  itsBunch->Bf[ii](0)*cosTemp_deltaPhi  - itsBunch->Bf[ii](1)*sinTemp_deltaPhi;
						double temp_E2 =  itsBunch->Bf[ii](0)*sinTemp_deltaPhi  + itsBunch->Bf[ii](1)*cosTemp_deltaPhi;
						itsBunch->Bf[ii](0) = temp_E1;  // Bx,T
						itsBunch->Bf[ii](1) = temp_E2;  // By,T
					}
				}
				
			}
			else{
				// if field solver is not available , only update bunch, to transfer particles between nodes if needed,
				// reset parameters such as LocalNum, initialTotalNum_m.
				
				Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
				Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					for (int jj=0; jj<3; jj++){
						temp_meanR(jj) += itsBunch->R[ii](jj);
						temp_meanP(jj) += itsBunch->P[ii](jj);
					}
				}

				// sum all particles in the same bin using reduce() for parallel enviroment,
				// and change N to itsBunch->getTotalNum().
				reduce(temp_meanR,temp_meanR,OpAddAssign());
				reduce(temp_meanP,temp_meanP,OpAddAssign());
				temp_meanR /= Vector_t(itsBunch->getTotalNum()); // in global cartesian frame :meanXX
				temp_meanP /= Vector_t(itsBunch->getTotalNum());
				
				// in global Cartesian frame, calculate the location in global frame of bunch
				double temp_meanTheta = calculateAngle2(temp_meanR(0), temp_meanR(1));
							
				// both for single bunch and multi-bunch
				if ( step_m > 10 ) // avoid dump at the first step
					if((oldReferenceTheta < initialReferenceTheta-deltaTheta) 
					   && (temp_meanTheta >= initialReferenceTheta-deltaTheta )) {
						*gmsg<< "one turn finished!"<<endl;
						turnnumber++;
						dumpEachTurn = true;
					}
				
				oldReferenceTheta = temp_meanTheta;
				// if no space charge effects are included, don't need to call update() in local frame
			}
			// track all particles for one step
			IpplTimings::startTimer(IntegrationTimer_m);
			
			size_t *countLost;
			if(itsBunch->weHaveBins()){
				const int tempN = itsBunch->pbin_m->getLastemittedBin();
				countLost = new size_t[tempN];
				for (int ii=0; ii< tempN; ii++) countLost[ii] = 0;
			}

			for (int i = 0; i < ( itsBunch->getLocalNum()); i++) {
				flagNoDeletion = true;
				// change phase space parameters from localframe of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
				for (int j =0;j<3;j++) {
					variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
					variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
					rold_m[j] = variable_m[j]; // used for gap cross checking
					pold_m[j] = variable_m[j+3]; // used for gap cross
				}

				// integrate for one step in the lab Cartesian frame (absulate value ).
				// IpplTimings::startTimer(IntegrationTimer_m);
				flagNoDeletion = rk4( variable_m, t, dt, i );
				if (!flagNoDeletion){
					// put particle onto deletion list .
					itsBunch->destroy( 1, i, true);
					
					//update bin parameter
					if(itsBunch->weHaveBins()) countLost[itsBunch->Bin[i]] += 1 ;

					*gmsgAll <<"PARTICLE DELETION(out of region): ID = "<<itsBunch->ID[i]<<", STEP = "<< step_m <<endl;
					
					flagNeedUpdate = true;
					continue;
				}
				
				// IpplTimings::stopTimer(IntegrationTimer_m);
				for (int j =0;j<3;j++) {
					itsBunch->R[i](j) = variable_m[j];    //[x,y,z]  units: [mm]
					itsBunch->P[i](j) = variable_m[j+3];  //[px,py,pz]  units: dimensionless, beta*gama
				}
				
				//If gap crossing happens, do momenta kicking
				for (beamline_list::iterator sindex = ++(FieldDimensions.begin()); 
					 sindex != FieldDimensions.end(); sindex++) {
					
					bool tag_crossing = false;
					bool tag_septum = false;
					double DistOld=0.0; //mm

					if( ((*sindex)->first) == "CAVITY") {
						// here check gap cross in the list, if do , set tag_crossing to TRUE
						for (int j =0;j<3;j++)
							rnew_m[j] = variable_m[j];
						tag_crossing = checkGapCross(rold_m, rnew_m, ((*sindex)->second).second, DistOld);
					}
					if(tag_crossing){
						double oldMomentum2  = dot(pold_m,pold_m); 
						double oldBetgam = sqrt(oldMomentum2);
						double oldGamma = sqrt(1.0 +oldMomentum2);
						double oldBeta = oldBetgam/oldGamma;
						double dt1 = DistOld/( c* oldBeta * 1.0e-6);  // ns
						double dt2 = dt -dt1;
						double tempP[3];

						// retrack particle from the old postion to cavity gap point
						// restore the old coordinates and momenta
						for (int j =0;j<3;j++) { 
							variable_m[j]=rold_m[j];
							variable_m[j+3]=pold_m[j];
						}
						
						if (dt/dt1 < 1000.0) rk4( variable_m, t, dt1, i );
						
						for (int j =0;j<3;j++) {
							itsBunch->R[i](j) = variable_m[j] ;    //[x,y,z]  units: [mm]
							itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units:[] beta*gama
						}
						//momentum kick
						double radius = sqrt( pow(variable_m[0], 2.0) + pow(variable_m[1], 2.0) 
											  - pow((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getPerpenDistance() , 2.0) );
						double rmin = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmin());
						double rmax = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmax());
						double nomalRadius = (radius - rmin)/(rmax - rmin);
						if (nomalRadius < 0 || nomalRadius > 1) {
							// *gmsg<<"Current radial location of particle is out of region NO energy gain" << endl;
						}
						else{
							for (int j =0;j<3;j++) 
								tempP[j] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
							
							// here evaluate voltage and conduct momenta kicking;
							(dynamic_cast<RFCavity*>(((*sindex)->second).second)) -> 
								getMomentaKick(nomalRadius, tempP, t, dt1, itsBunch->ID[i],itsReference.getM()); // t : ns

							for (int j =0;j<3;j++) 
								itsBunch->P[i](j) = tempP[j];
						}
						
						// retrack particle  from cavity gap point for the left time to finish the entire timestep
						
						for (int j =0;j<3;j++) {
							variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
							variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
						}
						
						if (dt/dt2 < 1000.0) rk4( variable_m, t, dt2, i );

						for (int j =0;j<3;j++){ 
							itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
							itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [] beta*gama
						}
					} // end if: gap-crossing monentum kicking at certain cavity
				} //end for: finish checking for all cavities
			} //end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode
			
			for (beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++){
				if( ((*sindex)->first) == "SEPTUM")    {
					bool flag_septum= (dynamic_cast<Septum*>(((*sindex)->second).second)) -> checkSeptum(*itsBunch);
				}
				
				if( ((*sindex)->first) == "PROBE")    {
					bool flag_probe= (dynamic_cast<Probe*>(((*sindex)->second).second)) 
						-> checkProbe(*itsBunch,turnnumber,deltaTheta);
				}
			}
			bool globPartOutOfBounds = (min(itsBunch->Bin) < 0);
			size_t ne = 0;
			if(globPartOutOfBounds) {
				ne = itsBunch->boundp_destroyT();
				flagNeedUpdate = true;
			}
			reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

			// update immediately if some particle are lost during this step
			// if have Field Solver, deletion will be done during boundp() before next track step, 
			// cause here is in global frame
			if(flagNeedUpdate){
				Vector_t temp_meanR = Vector_t(0.0,0.0,0.0);
				Vector_t temp_meanP = Vector_t(0.0,0.0,0.0);
				
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++)
					for (int jj=0; jj<3; jj++)
						temp_meanP(jj) += itsBunch->P[ii](jj);
				
				reduce(temp_meanP,temp_meanP,OpAddAssign());
				temp_meanP /= Vector_t(itsBunch->getTotalNum());

				// in global Cartesian frame, calculate the direction of longitudinal angle of bunch

				double temp_meanPhi = calculateAngle(temp_meanP(0), temp_meanP(1));

				temp_meanPhi -= pi/2.0;

				double cosTemp_meanPhi = cos(temp_meanPhi);
				double sinTemp_meanPhi = sin(temp_meanPhi);

				// remove mean coordinates
				itsBunch->R -= temp_meanR;

				//scale coordinates
				itsBunch->R /= Vector_t(1000.0); // mm --> m

				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  + itsBunch->R[ii](1)*sinTemp_meanPhi;
					double temp_RLongitudinal = -itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
					
					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;
				}

				// now destroy particles and update pertinent parameters in local frame
				itsBunch->boundp();

				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
					double temp_RHorizontal   =  itsBunch->R[ii](0)*cosTemp_meanPhi  - itsBunch->R[ii](1)*sinTemp_meanPhi;
					double temp_RLongitudinal =  itsBunch->R[ii](0)*sinTemp_meanPhi  + itsBunch->R[ii](1)*cosTemp_meanPhi;
					
					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;
				}

				//scale coordinates back
				itsBunch->R *= Vector_t(1000.0); // m --> mm

				// retrieve mean coordinates
				itsBunch->R += temp_meanR;

				*gmsg <<"Current total particles number is : "<< itsBunch->getTotalNum()<<endl;
			}

			if(itsBunch->weHaveBins() && flagNeedUpdate )
				itsBunch->pbin_m->updatePartInBin(countLost);
		    
			if ( (itsBunch->weHaveBins()) && BunchCount_m > 1 && step_m%resetBinFreq == 0 )
				itsBunch->resetPartBinID();

			if (itsBunch->weHaveBins() && countLost!=NULL) delete[] countLost;

			IpplTimings::stopTimer(IntegrationTimer_m);
			Ippl::Comm->barrier();

		} 
		else if (initialTotalNum_m == 2) {
			// initialTotalNum_m == 2
			// trigger SEO mode (swith off cavity) and calculate betatron osciliation tuning.
			double r_tuning[2], z_tuning[2] ;

			for (int i = 0; i < ( itsBunch->getLocalNum()); i++) {
				
				// change phase space parameters from local frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
				for (int j =0;j<3;j++)variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
				for (int j =0;j<3;j++)variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: [] 
	
				if ( (step_m%SinglePartDumpFreq == 0) ){
					outf<<"ID"<<( itsBunch->ID[i] );
					outf<<" "<<variable_m[0]<<" "<<variable_m[3]<<" "<<variable_m[1]<<" "
						<<variable_m[4]<<" "<<variable_m[2]<<" "<<variable_m[5]<<endl;
				}
				
				double OldTheta=0.0;

				OldTheta = calculateAngle(variable_m[0], variable_m[1]);
				r_tuning[i] = variable_m[0]*cos(OldTheta)+variable_m[1]*sin(OldTheta);
				z_tuning[i] = variable_m[2];
				turnnumber=lastTurn;

				// integrate for one step in the lab Cartesian frame (absulate value ).
				
				rk4( variable_m, t, dt, i );

				double NewTheta=0.0;

				NewTheta = calculateAngle(variable_m[0], variable_m[1]);
				if ( (i==0) && (NewTheta < OldTheta)) lastTurn++;
				
				for (int j =0;j<3;j++) itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
				for (int j =0;j<3;j++) itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: dimensionless, beta*gama

			}//end for: finish one step tracking for all particles for initialTotalNum_m != 2 mode
			
			// store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
			if ( step_m%SinglePartDumpFreq == 0 ){ 
				Ttime.push_back(t*1.0e-9);
				Tdeltz.push_back(z_tuning [1]);
				Tdeltr.push_back(r_tuning[1]  - r_tuning[0]);
				TturnNumber.push_back(turnnumber);
			}
		}
		else if (initialTotalNum_m == 1) {
			// initialTotalNum_m == 1 trigger single particle mode
			int i = 0; flagNoDeletion = true;
			
			// change phase space parameters from local frame of bunch (dr,dtheta,dz) to global Cartesian frame (X,Y,Z)
			for (int j =0;j<3;j++) {
				variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
				variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: dimensionless
				rold_m[j] = variable_m[j]; // used for gap cross checking
				pold_m[j] = variable_m[j+3]; // used for gap cross
			}
			
			if ( (step_m%SinglePartDumpFreq == 0) ){
				outf<<"ID"<<( itsBunch->ID[i] );
				outf<<" "<<variable_m[0]<<" "<<variable_m[3]<<" "<<variable_m[1]<<" "
					<<variable_m[4]<<" "<<variable_m[2]<<" "<<variable_m[5]<<endl;
			}

			double temp_meanTheta = calculateAngle2(variable_m[0], variable_m[1]);//[ -pi ~ pi ]

			if((oldReferenceTheta < initialReferenceTheta-deltaTheta) 
			   && (temp_meanTheta >= initialReferenceTheta-deltaTheta )){
				turnnumber++;
				dumpEachTurn = true; 
				outfThetaEachTurn_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfThetaEachTurn_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
								 <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
								 <<" "<<temp_meanTheta/pi*180
								 <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
								 <<" "<<variable_m[2]
								 <<" "<<variable_m[5]<<endl;
			}
			
			//define 3 special azimuthal angles where dump particle's six parameters  at each turn into 3 ASCII files.
			const double azimuth_angle0=0.0;
			const double azimuth_angle1=22.5/180.0*pi;
			const double azimuth_angle2=45.0/180.0*pi;

			if((oldReferenceTheta < azimuth_angle0-deltaTheta) && (temp_meanTheta >= azimuth_angle0-deltaTheta )){
				outfTheta0_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta0_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}
			
			if((oldReferenceTheta < azimuth_angle1-deltaTheta) && (temp_meanTheta >= azimuth_angle1-deltaTheta )){
				outfTheta1_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta1_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}

			if((oldReferenceTheta < azimuth_angle2-deltaTheta) && (temp_meanTheta >= azimuth_angle2-deltaTheta )){
				outfTheta2_m<< "#Turn number = "<<turnnumber<<", Time = "<<t<<" [ns]"<<endl;
				outfTheta2_m<<" "<<sqrt(variable_m[0]*variable_m[0]+variable_m[1]*variable_m[1])
						  <<" "<<variable_m[3]*cos(temp_meanTheta)+variable_m[4]*sin(temp_meanTheta)
						  <<" "<<temp_meanTheta/pi*180
						  <<" "<<-variable_m[3]*sin(temp_meanTheta)+variable_m[4]*cos(temp_meanTheta)
						  <<" "<<variable_m[2]
						  <<" "<<variable_m[5]<<endl;
			}

			oldReferenceTheta = temp_meanTheta;

			// integrate for one step in the lab Cartesian frame (absulate value ).
			// IpplTimings::startTimer(IntegrationTimer_m);
			flagNoDeletion = rk4( variable_m, t, dt, i );

			if (!flagNoDeletion){
				*gmsg <<"particle"<<"is lost at "<< step_m <<"st step!"<<endl;
				throw OpalException("ParallelCyclotronTracker", "the particle is out of the region of interest.");
			}
			// IpplTimings::stopTimer(IntegrationTimer_m);
			for (int j =0;j<3;j++) itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
			for (int j =0;j<3;j++) itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [] beta*gama

			//If gap crossing happens, do momenta kicking

			for ( beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); sindex++){
				bool tag_crossing = false;
				double DistOld=0.0; //mm

				if( ((*sindex)->first) == "CAVITY") {
					// here check gap cross in the list, if do , set tag_crossing to TRUE
					for (int j =0;j<3;j++)
						rnew_m[j] = variable_m[j];
					
					tag_crossing = checkGapCross(rold_m, rnew_m, ((*sindex)->second).second, DistOld);
				}

				if(tag_crossing){
					*gmsg << "***------------------------ gap crossing happen 2 ----------------------------*** "<<endl;
					
					double oldMomentum2  = dot(pold_m,pold_m); 
					double oldBetgam = sqrt(oldMomentum2);
					double oldGamma = sqrt(1.0 +oldMomentum2);
					double oldBeta = oldBetgam/oldGamma;
					double dt1 = DistOld/( c* oldBeta * 1.0e-6);  // ns
					double dt2 = dt -dt1;
					double tempP[3];
					
					// retrack particle from the old postion to cavity gap point
					
					// restore the old coordinates and momenta
					for (int j =0;j<3;j++) { 
						variable_m[j]=rold_m[j];
						variable_m[j+3]=pold_m[j];
					}
					
					if (dt/dt1 < 1000.0) rk4( variable_m, t, dt1, i );

					for (int j =0;j<3;j++) { 
						itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
						itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [] beta*gama
					}

					//momentum kick

					double radius=sqrt(pow(variable_m[0], 2.0) + 
									   pow(variable_m[1], 2.0) -
									   pow((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getPerpenDistance(),2.0));
					double rmin = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmin());
					double rmax = ((dynamic_cast<RFCavity*>(((*sindex)->second).second))->getRmax());
					double nomalRadius = (radius - rmin)/(rmax - rmin);
					if (nomalRadius < 0 || nomalRadius > 1) 
						*gmsg<<"Current radial location of particle is out of range." << endl;
					else {
						for (int j =0;j<3;j++) 
							tempP[j] = itsBunch->P[i](j);  //[px,py,pz]  units: [] 

						// here evaluate voltage and conduct momenta kicking;
						(dynamic_cast<RFCavity*>(((*sindex)->second).second)) -> 
							getMomentaKick(nomalRadius, tempP, t, dt1, itsBunch->ID[i],itsReference.getM()); // t : ns

						for (int j =0;j<3;j++) 
							itsBunch->P[i](j) = tempP[j];
						
					}

					// retrack particle  from cavity gap point for the left time to finish the entire timestep

					for (int j =0;j<3;j++) {
						variable_m[j] = itsBunch->R[i](j);  //[x,y,z]  units: [mm]
						variable_m[j+3] = itsBunch->P[i](j);  //[px,py,pz]  units: []
					}
					
					// FixME: ??
					if (dt/dt2 < 1000.0) rk4( variable_m, t, dt2, i );
							  
					for (int j =0;j<3;j++) {
						itsBunch->R[i](j) = variable_m[j] ;  //[x,y,z]  units: [mm]
						itsBunch->P[i](j) = variable_m[j+3] ;  //[px,py,pz]  units: [], beta*gama
					}
				}// end if: gap-crossing monentum kicking at certain cavity
			}//end for: finish checking for all cavities
		}//end if: finish one step tracking either for initialTotalNum_m==2 || initialTotalNum_m==2 || initialTotalNum_m==1 mode
		
		// update bunch and some parameters and output some info. after one time step.
		// reference parameters
		// pow(itsBunch->P[0](0),2.0) + pow(itsBunch->P[0](1),2.0) + pow(itsBunch->P[0](2),2.0);
		double tempP2= dot(itsBunch->P[0],itsBunch->P[0]); 
		double tempGamma = sqrt( 1.0 + tempP2 );
		double tempBeta = sqrt(tempP2)/tempGamma;

		PathLength_m += c_mmtns*dt/1000.0*tempBeta; // unit: m
	    
		t += dt;
		itsBunch->setT((t)*1.0e-9 );
		itsBunch->setLPath(PathLength_m);
		// Here is global frame, don't do itsBunch->boundp();
	    
		// FixMe:
		// dump phase space distribution of bunch
		// if ((((step_m+1) % Options::psDumpFreq == 0 ) && initialTotalNum_m > 2) || (dumpEachTurn && initialTotalNum_m > 2) )
		// if(step_m>(1185*(Options::rDump))) {
		//  itsDataSink->writeStatData(*itsBunch, FDext_m ,0.0,0.0,0.0);
		//}
		if ((((step_m+1) % Options::psDumpFreq==0) && initialTotalNum_m != 2) 
			|| (doDumpAfterEachTurn && dumpEachTurn && initialTotalNum_m != 2)) {
			IpplTimings::startTimer(DumpTimer_m);

			itsBunch->setTrackStep((step_m+1));

			extE_m = Vector_t(0.0, 0.0, 0.0);
			extB_m = Vector_t(0.0, 0.0, 0.0);

			//--------------------- calculate mean coordinates  of bunch -------------------------------//
			//------------  and calculate the external field at the mass of bunch-----------------------//
			
			Vector_t meanR=Vector_t(0.0,0.0,0.0);
			Vector_t meanP=Vector_t(0.0,0.0,0.0);
			
			for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				for (int j=0; j<3; j++){
					meanR(j) += itsBunch->R[ii](j);
					meanP(j) += itsBunch->P[ii](j);
				}
			}
			
			reduce(meanR, meanR,OpAddAssign());
			meanR /= Vector_t(( itsBunch->getTotalNum()));

			reduce(meanP, meanP,OpAddAssign());
			meanP /= Vector_t(( itsBunch->getTotalNum()));

			double temp_Phi = calculateAngle(meanP(0), meanP(1));
			// the angle from positive radial direction to propagational  direction
			//                                   ==>
			// the angle from positive azimuthal direction to propagational direction
			temp_Phi -= pi/2.0;

			double temp_cosPhi = cos(temp_Phi);
			double temp_sinPhi = sin(temp_Phi);
			
			double temp_Theta = calculateAngle(meanR(0), meanR(1));
			double temp_cosTheta = cos(temp_Theta);
			double temp_sinTheta = sin(temp_Theta);

			double meanRadius = sqrt( meanR(0)*meanR(0) +  meanR(1)*meanR(1) );

			beamline_list::iterator DumpSindex = FieldDimensions.begin();

			if (((((*DumpSindex)->second).first)[0] <= meanRadius ) 
				&& ((((*DumpSindex)->second).first)[1] >= meanRadius)) {
				( ((*DumpSindex)->second).second )->apply(meanR, Vector_t(0.0), t, extE_m, extB_m);
			}
		    
			FDext_m[0] = extB_m/10.0; // kgauss -> T
			FDext_m[1] = extE_m;

			//----------------------------dump in global frame-------------------------------------//
			// Note: Don't dump when
			// 1. after one turn
			// 2.for multi-bunch mode, just after backup.
			// in order to sychronize the dump step for multi-bunch and single 
			// bunch for compare with each other during post-process phase.
			if ( !(Options::psDumpLocalFrame) ) {
				if(!OneStepAfterBackUp ){
					
					itsBunch->R /= Vector_t(1000.0); // mm --> m
					lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m );
					itsBunch->R *= Vector_t(1000.0); // m --> mm

					*gmsg <<" Phase space dumping in global frame after step "<<step_m+1<<endl;
					*gmsg <<" Dump step = "<<lastDumpedStep_m<<", time = " << t << " [ns]" << endl;
				}
				else {
					// only omit the dump step just after backup
					*gmsg << "Omit this Dump step = "<<lastDumpedStep_m<<endl;
					OneStepAfterBackUp = false;
				}
				//----------------------------dump in local frame-------------------------------------//
			} 
			else {
				itsBunch->R -= meanR;
				itsBunch->P -= meanP;
				/*
				  for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				  
				  double temp_RHorizontal     =   itsBunch->R[ii](0)*temp_cosPhi  + itsBunch->R[ii](1)*temp_sinPhi;
				  double temp_RLongitudinal   =  -itsBunch->R[ii](0)*temp_sinPhi  + itsBunch->R[ii](1)*temp_cosPhi;
				  
				  itsBunch->R[ii](0) = temp_RLongitudinal;
				  itsBunch->R[ii](1) = temp_RHorizontal;
				  
				  double temp_PHorizontal     =   itsBunch->P[ii](0)*temp_cosPhi  + itsBunch->P[ii](1)*temp_sinPhi;
				  double temp_PLongitudinal   =  -itsBunch->P[ii](0)*temp_sinPhi  + itsBunch->P[ii](1)*temp_cosPhi;
				  
				  itsBunch->P[ii](0) = temp_PLongitudinal;
				  itsBunch->P[ii](1) = temp_PHorizontal;
				  }
				*/
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){

					double temp_RHorizontal     =   itsBunch->R[ii](0)*temp_cosTheta  + itsBunch->R[ii](1)*temp_sinTheta;
					double temp_RLongitudinal   =  -itsBunch->R[ii](0)*temp_sinTheta  + itsBunch->R[ii](1)*temp_cosTheta;

					itsBunch->R[ii](0) = temp_RHorizontal;
					itsBunch->R[ii](1) = temp_RLongitudinal;

					double temp_PHorizontal     =   itsBunch->P[ii](0)*temp_cosTheta  + itsBunch->P[ii](1)*temp_sinTheta;
					double temp_PLongitudinal   =  -itsBunch->P[ii](0)*temp_sinTheta  + itsBunch->P[ii](1)*temp_cosTheta;

					itsBunch->P[ii](0) = temp_PHorizontal;
					itsBunch->P[ii](1) = temp_PLongitudinal;
				}

				itsBunch->R /= Vector_t(1000.0); // mm --> m
				lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(*itsBunch, FDext_m );
				itsBunch->R *= Vector_t(1000.0); // m --> mm
	
				/*
				  for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){
				  
				  double temp_Rx   =  -itsBunch->R[ii](0)*temp_sinPhi  + itsBunch->R[ii](1)*temp_cosPhi;
				  double temp_Ry   =   itsBunch->R[ii](0)*temp_cosPhi  + itsBunch->R[ii](1)*temp_sinPhi;

				  itsBunch->R[ii](0) = temp_Rx;
				  itsBunch->R[ii](1) = temp_Ry;
				  
				  double temp_Px   =  -itsBunch->P[ii](0)*temp_sinPhi  + itsBunch->P[ii](1)*temp_cosPhi;
				  double temp_Py   =   itsBunch->P[ii](0)*temp_cosPhi  + itsBunch->P[ii](1)*temp_sinPhi;
				  
				  itsBunch->P[ii](0) = temp_Px;
				  itsBunch->P[ii](1) = temp_Py;

				  }

				*/
				for (int ii=0; ii<( itsBunch->getLocalNum()); ii++){

					double temp_Rx   =  -itsBunch->R[ii](1)*temp_sinTheta  + itsBunch->R[ii](0)*temp_cosTheta;
					double temp_Ry   =   itsBunch->R[ii](1)*temp_cosTheta  + itsBunch->R[ii](0)*temp_sinTheta;

					itsBunch->R[ii](0) = temp_Rx;
					itsBunch->R[ii](1) = temp_Ry;

					double temp_Px   =  -itsBunch->P[ii](1)*temp_sinTheta  + itsBunch->P[ii](0)*temp_cosTheta;
					double temp_Py   =   itsBunch->P[ii](1)*temp_cosTheta  + itsBunch->P[ii](0)*temp_sinTheta;

					itsBunch->P[ii](0) = temp_Px;
					itsBunch->P[ii](1) = temp_Py;

				}

				itsBunch->R += meanR;
				itsBunch->P += meanP;

				*gmsg <<" Phase space dumping in local frame after step "<<step_m+1<<endl;
				*gmsg <<" Dump step = "<<lastDumpedStep_m<<", time = " << t << " [ns]" << endl;
			}
			IpplTimings::stopTimer(DumpTimer_m);
			
		}
		if (!(step_m+1 % 1000))
			*gmsg << "Step " << step_m+1 << endl;

	}// end for: the integration is DONE after maxSteps_m steps!

	// some post-integration works
	*gmsg << "***---------------------------- PARTICLES TRACK DONE------ ----------------------------*** "<<endl;
	
	// calculate tunes after tracking.

	for(int ii=0; ii<( itsBunch->getLocalNum()); ii++) {
		if (itsBunch->ID[ii] == 0) {
			double FinalMomentum2= pow(itsBunch->P[ii](0),2.0) + pow(itsBunch->P[ii](1),2.0) + pow(itsBunch->P[ii](2),2.0);
			double FinalEnergy= (sqrt(1.0 +FinalMomentum2)-1.0)*itsReference.getM()*1.0e-6;
			*gmsgAll << "Final energy of reference particle = "<< FinalEnergy<<" [MeV]"<<endl;
			*gmsgAll << "Total phase space dump number(includes the initial distribution) = "<< lastDumpedStep_m+1 <<endl;
			*gmsgAll << "One can restart simulation from the last dump step ( -restart "<< lastDumpedStep_m <<" )"<<endl;
		}
	}

	Ippl::Comm->barrier();

	if (initialTotalNum_m == 2 ){
		*gmsg << "********** The result for tuning calulation (NO space charge) ********** "<<endl
			  << "Total turn number is "<< TturnNumber.back() <<endl;
		bool tuneflag = false;
		double nur,nuz;
		tuneflag = getTunes(Ttime,Tdeltr,Tdeltz,TturnNumber.back(),nur,nuz);
	}
	else if(flagNoDeletion && initialTotalNum_m > 2){
		if(myNode_m == 0 && flagDoTune_m) {
			*gmsg << "********** The result for tuning calulation (WITH space charge) ********** "<<endl
				  << "Total turn number is "<< TturnNumber.back() <<endl;
			bool tuneflag = false;
			double nur,nuz;
			int  No_element =(int)(TdeltrSC[0].size());
			for (int ii=0; ii< NoSamplePart; ii++){
			 	
				double dR_total = accumulate(TdeltrSC[ii].begin(), TdeltrSC[ii].end(), 0); 
				double dR_average = (dR_total / No_element ) ;
				*gmsg<<" average dR = "<<dR_average<<endl;
				double dZ_total = accumulate(TdeltzSC[ii].begin(), TdeltzSC[ii].end(), 0); 
				double dZ_average = (dZ_total / No_element );
				*gmsg<<" average dZ = "<<dZ_average<<endl;
				vector<double>::iterator tempIt; 
				
				for ( tempIt = TdeltrSC[ii].begin(); tempIt != TdeltrSC[ii].end(); tempIt++ ){
					*tempIt= *tempIt - dR_average;
				}
				
				for ( tempIt = TdeltzSC[ii].begin(); tempIt != TdeltzSC[ii].end(); tempIt++ ){
					*tempIt= *tempIt - dZ_average;
				}
				
				tuneflag = getTunes(Ttime,TdeltrSC[ii], TdeltzSC[ii], TturnNumber.back(),nur,nuz);
			
			}
		}
	}
	else{
		// not for multibunch
		if (!(itsBunch->weHaveBins()))
			*gmsg << "Total finished turn number (not correct for restart mode) = "<< turnnumber <<endl;
	}
	
	Ippl::Comm->barrier();

	if(myNode_m == 0) outf.close();

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
bool ParallelCyclotronTracker::derivate(double *y, double t, double *yp, int Pindex )
{
    Vector_t externalE, externalB, tempR;

    double partR;

    // a flag of visit flag, if cyclotron is not visited, set false.
    // bool visitflag = false;


    externalB = Vector_t(0.0, 0.0, 0.0);
    externalE = Vector_t(0.0, 0.0, 0.0);

    for(int i=0;i<3;i++) 
		tempR(i)=y[i];

    partR = sqrt( tempR(0)*tempR(0) + tempR(1)*tempR(1) );

    beamline_list::iterator sindex = FieldDimensions.begin();

    //debug
    if ( ((*sindex)->first) !="CYCLOTRON"){
		*gmsg<<"Error in ParallelCyclotronTracker, the CYCLOTRON object is not the first element in the beamline_list! "
			 << endl;
		exit(1);
    }
	
    if (((((*sindex)->second).first)[0] <= partR ) &&  ( (((*sindex)->second).first)[1] >= partR))
		( ((*sindex)->second).second )->apply(tempR, Vector_t(0.0), t, externalE, externalB);
	else 
		return false;
    
    for(int i=0;i<3;i++) externalB(i) = externalB(i) * 0.10;  //[kGauss] -> [T]
    for(int i=0;i<3;i++) externalE(i) = externalE(i) * 1.0e6;  //[kV/mm ] -> [V/m]
	
    // for working modes without space charge effects, override this step to save time
    if(itsBunch->hasFieldSolver()) {
		if ( itsBunch->ID[Pindex] != 0 ){
			// add external Field and self space charge field
			externalE += itsBunch->Ef[Pindex];
			externalB += itsBunch->Bf[Pindex];
		}
	    
	}

    double tempgamma = sqrt(1 + (y[3]*y[3]+y[4]*y[4]+y[5]*y[5]));

    /* 
	   d2x/dt2 = q/m * ( E + v x B )
       dx/dt =vx
       units: mm, and []
    */
	
    yp[0] = c_mmtns/tempgamma * y[3];  // [mn/ns]
    yp[1] = c_mmtns/tempgamma * y[4];  // [mm/ns]
    yp[2] = c_mmtns/tempgamma * y[5];  // [mm/ns]

    yp[3] = (externalE(0)*chtmc_m  + (externalB(2)*y[4] - externalB(1)*y[5])*chtm_m/tempgamma )*1.0e-9; // [1/ns]
    yp[4] = (externalE(1)*chtmc_m  - (externalB(2)*y[3] - externalB(0)*y[5])*chtm_m/tempgamma )*1.0e-9; // [1/ns];
    yp[5] = (externalE(2)*chtmc_m  + (externalB(1)*y[3] - externalB(0)*y[4])*chtm_m/tempgamma )*1.0e-9; // [1/ns];
	
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
bool ParallelCyclotronTracker::rk4(double x[],double t,double tau,int Pindex)
{
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

	visitflag=derivate( x, t, deriv1 ,Pindex);

    if (!visitflag) return false;

    // Evaluate f2 = f( x+tau*f1/2, t+tau/2 ).
	double half_tau = 0.5*tau;
    double t_half = t + half_tau;

    for(int i=0; i<PSdim; i++ )
		xtemp[i] = x[i] + half_tau*deriv1[i];

    visitflag=derivate( xtemp, t_half, deriv2 ,Pindex);

    if (!visitflag) return false;

    // Evaluate f3 = f( x+tau*f2/2, t+tau/2 ).
	for(int i=0; i<PSdim; i++ )
		xtemp[i] = x[i] + half_tau*deriv2[i];

    visitflag=derivate( xtemp, t_half, deriv3 ,Pindex);
	
    if (!visitflag) return false;

    // Evaluate f4 = f( x+tau*f3, t+tau ).
	double t_full = t + tau;
    for(int i=0; i<PSdim; i++ )
		xtemp[i] = x[i] + tau*deriv3[i];

    visitflag=derivate( xtemp, t_full, deriv4 ,Pindex);

    if (!visitflag) return false;

    // Return x(t+tau) computed from fourth-order R-K.
	for(int i=0; i<PSdim; i++ )
		x[i] += tau/6.*(deriv1[i] + deriv4[i] + 2.*(deriv2[i]+deriv3[i]));

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

bool ParallelCyclotronTracker::checkGapCross(Vector_t Rold, Vector_t Rnew, Component *elptr,double &Dold)
{
    bool flag = false;

    double sinx = (dynamic_cast<RFCavity*>(elptr))->getSinAzimuth();
    double cosx = (dynamic_cast<RFCavity*>(elptr))->getCosAzimuth();
    double PerpenDistance = (dynamic_cast<RFCavity*>(elptr))->getPerpenDistance();

    double distNew = (Rnew[0]*sinx - Rnew[1]*cosx) - PerpenDistance;
    double distOld = (Rold[0]*sinx - Rold[1]*cosx) - PerpenDistance;

    if (distOld > 0.0 && distNew <=0.0) flag = true;
    // This parameter is used correct cavity phase
    Dold = distOld;
    return flag;
}


struct adder : public unary_function<double, void>
{
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
bool ParallelCyclotronTracker::getTunes(vector<double> &t,vector<double> &r,vector<double> &z,
										int lastTurn,double &nur,double &nuz)
{
    Inform msg1("getTunes ");

    TUNE_class *tune;

    int Ndat = t.size();

    /*
      remove mean
    */
    double rsum =  for_each(r.begin(), r.end(), adder()).sum;

    for (int i=0; i< Ndat; i++)
		r[i] -= rsum;

    double zsum =  for_each(z.begin(), z.end(), adder()).sum;

    for (int i=0; i< Ndat; i++)
		z[i] -= zsum;
    double ti = *(t.begin());
    double tf = t[t.size()-1];
    double T = (tf-ti);

    t.clear();
    double dt = T/Ndat;
    double time=0.0;
    for (int i=0; i< Ndat; i++) {
		t.push_back(i);
		time+=dt;
    }
	
    T = t[Ndat-1];
	
    msg1 << " *************** nuR ***************" << endl;
    msg1 << endl<< "===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;
	
    int nhis_lomb = 10;
    int  stat = 0;
    // book tune class
    tune = new TUNE_class();
    stat = tune->LombAnalysis(t,r,nhis_lomb,T/lastTurn);
    if(stat != 0)
		msg1 << "TUNE: Lomb analysis failed" << endl;
    msg1 << " ***********************************" << endl << endl;
	
    delete tune;
    tune = NULL;
	// FixMe: need to come from the inputfile
    nhis_lomb = 100;
	
    if (zsum != 0.0) {
		msg1 << " *************** nuZ ***************" << endl;
		msg1 << endl<< "===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

		// book tune class
		tune = new TUNE_class();
		stat = tune->LombAnalysis(t,z,nhis_lomb,T/lastTurn);
		if(stat != 0) msg1 << "TUNE: Lomb analysis failed" << endl;
		msg1 << " ***********************************" << endl << endl;

		delete tune;
		tune = NULL;
    }
    return true;
}

/** 
 * 
 * 
 * @param x 
 * @param y 
 * 
 * @return angle range [0~2PI) degree
 */
double ParallelCyclotronTracker::calculateAngle(double x, double y)
{
    double thetaXY;

    if (x<0)                   thetaXY=pi+atan(y/x);
    else if ((x>0) && (y>=0))  thetaXY=atan(y/x);
    else if ((x>0) && (y<0))   thetaXY=2.0*pi+atan(y/x);
    else if ((x==0) && (y> 0)) thetaXY=pi/2.0;
    else if ((x==0) && (y< 0)) thetaXY=3.0/2.0*pi;

    return thetaXY;

}

// angle range [-PI~PI) degree
double ParallelCyclotronTracker::calculateAngle2(double x, double y)
{

    double thetaXY;
    if      (x>0)              thetaXY=atan(y/x);
    else if ((x<0)  && (y>0)) thetaXY=pi+atan(y/x);
    else if ((x<0)  && (y<=0)) thetaXY=-pi+atan(y/x);
    else if ((x==0) && (y> 0)) thetaXY=pi/2.0;
    else if ((x==0) && (y< 0)) thetaXY=-pi/2.0;

    return thetaXY;

}

/** 
 * need to go to DataSink
 * 
 * @param BinID 
 * @param step 
 */
void ParallelCyclotronTracker::readOneBunch(const int BinID, const int step)
{
    *gmsg<<"---------------- Start reading hdf5 file----------------"<<endl;
    H5PartFile *H5file;
    string fn;

    if (OPAL.hasRestartFile()){
		fn=OPAL.getRestartFileName();
		if (fn.find(string("_Part"),0) == string::npos){
			int pos = fn.find(string(".h5"),0);
			fn.replace( pos, fn.length(), "_Part00002.h5");
		}
		else{
			int pos=fn.find(string("_Part"),0);
			string numstr=fn.substr(pos+5,5);
			int numint = atoi(numstr.c_str())+1;
			char numCStr[5];
			sprintf(numCStr, "%05d", numint);
			string numstrnew(numCStr);
			fn.replace( pos+5, 5, numstrnew);
		}
		
    }else{
		fn= OPAL.getInputFn();
		int pos=fn.find(string("."),0);
		fn.erase(pos,fn.size()-pos);
		fn += string(".h5");
    }
    *gmsg<<"Read phase space data(Step#"<<step<<") for new bunch from hdf5 format file "<<fn<<endl;
	
#ifdef PARALLEL_IO
    H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif
	
    if(!H5file) {
		ERRORMSG("File open failed:  exiting!" << endl);
		exit(0);
    }

    H5PartSetStep(H5file, step);
    const int globalN=(int)H5PartGetNumParticles(H5file);

    h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);

    int numberOfParticlesPerNode = (int) floor((double) globalN / Ippl::getNodes());
    long long starti = myNode_m * numberOfParticlesPerNode;
    long long endi = 0;
    // ensure that we dont miss any particle in the end
    if(myNode_m == Ippl::getNodes() - 1)
		endi = -1;
    else
		endi = starti + numberOfParticlesPerNode;

    H5PartSetView(H5file,starti,endi);
    const int InjectN = (int)H5PartGetNumParticles(H5file);

    void *varray = malloc(InjectN*sizeof(double));
    double *farray = (double*)varray;

    const int LocalNum = itsBunch->getLocalNum();
    const int NewLocalNum = LocalNum + InjectN;

    itsBunch->create(InjectN);

    H5PartReadDataFloat64(H5file,"x",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++) {
		itsBunch->R[ii](0)=farray[ii-LocalNum]*1000.0; //m-->mm
		// unlike the process in of restart, here we set bin index forcely,
		// because this new bunch should  be a new bin with lowest energy.
		itsBunch->Bin[ii] = BinID;
    }
    H5PartReadDataFloat64(H5file,"y",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++)
		itsBunch->R[ii](1)=farray[ii-LocalNum]*1000.0;

    H5PartReadDataFloat64(H5file,"z",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++)
		itsBunch->R[ii](2)=farray[ii-LocalNum]*1000.0;

    H5PartReadDataFloat64(H5file,"px",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++)
		itsBunch->P[ii](0)=farray[ii-LocalNum];

    H5PartReadDataFloat64(H5file,"py",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++)
		itsBunch->P[ii](1)=farray[ii-LocalNum];

    H5PartReadDataFloat64(H5file,"pz",farray);
    for (int ii=LocalNum; ii < NewLocalNum; ii++)
		itsBunch->P[ii](2)=farray[ii-LocalNum];

    // update the bin status
    //    if(itsBunch->weHaveBins())
    //       itsBunch->pbin_m->updateStatus(BinID+1, InjectN);

    // free memory
    if(farray)
		free(farray);

    Ippl::Comm->barrier();
    H5PartCloseFile(H5file);

    // update statistics parameters of PartBunch
    // allocate ID for new particles
    itsBunch->boundpNoRep();
    itsBunch->boundp();
    *gmsg<<"----------------Finish reading hdf5 file----------------"<<endl;
}

/** 
 * 
 * 
 * @param para 
 */
void ParallelCyclotronTracker::setTrackCoeff(const double para)
{
    ratioCh_M_m = para * m_p ;
    *gmsg<<"* Charge to Mass ratio = "<<ratioCh_M_m<<" [proton unit]."<<endl;
    chtmc_m *= ratioCh_M_m; // s*C/kg*m
    chtm_m  *= ratioCh_M_m; // C/kg
}
