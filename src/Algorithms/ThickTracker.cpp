// ------------------------------------------------------------------------
// $RCSfile: ThickTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickTracker
//   The visitor class for building a map of given order for a beamline
//   using a finite-length lenses for all elements.
//   Multipole-like elements are done by expanding the Lie series.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------


/*
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
*/


#include <typeinfo>
#include <fstream>
#include "Algorithms/ThickTracker.h"
#include "Algorithms/OrbitThreader.h" 
#include "Algorithms/CavityAutophaser.h"

#include <cfloat>

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"

#include "Fields/BMultipoleField.h"
#include "Classic/FixedAlgebra/FTps.h"
#include "Classic/FixedAlgebra/FTpsMath.h"
#include "Classic/FixedAlgebra/FVps.h"

#include "Classic/Fields/BSingleMultipoleField.h"


#include "Classic/Algorithms/PartData.h"  //for the beam reference


#include "Physics/Physics.h"
//#include "Utilities/NumToStr.h"

#include "Elements/OpalBeamline.h"
#include <Classic/BeamlineCore/MultipoleRep.h>
#include "Distribution/MapGenerator.h"

#define DIM 3

class Beamline;
class PartData;
using Physics::c;

#define PSdim 6
typedef FVector<double, PSdim> Vector;
typedef FMatrix<double, PSdim, PSdim> Matrix;
typedef FTps<double, PSdim> Series;
typedef FVps<double, PSdim> Map, VSeries;
typedef FMatrix<FTps<double, PSdim>, PSdim, PSdim> MxSeries;


/*
namespace {
    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf1, double ds,
                           int nx = 20);
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds,
                        int nx = 20, int cx = 4);
    Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds,
                          int nx = 50);
    Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds,
                          int nx = 50, int cx = 4);
};
*/

//
// Class ThickTracker
// ------------------------------------------------------------------------
//
ThickTracker::ThickTracker(const Beamline &beamline,
                           const PartData &reference,
                           bool revBeam, bool revTrack):
  Tracker(beamline, reference, revBeam, revTrack),
  itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getCoordTransformationTo()),
  pathLength_m(0.0),
  zStop_m(),
  dtCurrentTrack_m(0.0),
  dtAllTracks_m(),
  localTrackSteps_m()
{
}


ThickTracker::ThickTracker(const Beamline &beamline,
                           PartBunchBase<double, 3> *bunch,
			   DataSink &ds,
			   const PartData &reference,	  
                           bool revBeam, bool revTrack,
			   const std::vector<unsigned long long> &maxSteps,
			   double zstart,
			   const std::vector<double> &zstop,
			   const std::vector<double> &dt):
  Tracker(beamline, bunch, reference, revBeam, revTrack),
  itsDataSink_m(&ds), 
  itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getCoordTransformationTo()),
  pathLength_m(0.0),
  zstart_m(zstart),
  zStop_m(),
  dtCurrentTrack_m(0.0),
  dtAllTracks_m(),
  localTrackSteps_m()
{
  CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
				 beamline.getCoordTransformationTo());
  referenceToLabCSTrafo_m = labToRef.inverted();

  for (std::vector<unsigned long long>::const_iterator it = maxSteps.begin(); it != maxSteps.end(); ++ it) {
    localTrackSteps_m.push(*it);
  }
  for (std::vector<double>::const_iterator it = dt.begin(); it != dt.end(); ++ it) {
    dtAllTracks_m.push(*it);
  }
  for (std::vector<double>::const_iterator it = zstop.begin(); it != zstop.end(); ++ it) {
    zStop_m.push(*it);
  }
}


ThickTracker::~ThickTracker()
{}



void ThickTracker::visitBeamline(const Beamline &bl) {

    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    if (fbl->getRelativeFlag()) {
        *gmsg << " do stuff" << endl;
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getCoordTransformationTo());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
}


void ThickTracker::updateRFElement(std::string elName, double maxPhase) {

}


void ThickTracker::prepareSections() {
    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();
}



void ThickTracker::saveCavityPhases() {
    itsDataSink_m->storeCavityInformation();
}

void ThickTracker::restoreCavityPhases() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() ||
        OpalData::getInstance()->inRestartRun()) {
        iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++ it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}


void ThickTracker::autophaseCavities(const BorisPusher &pusher) {

    double t = itsBunch_m->getT();
    Vector_t nextR = RefPartR_m / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, RefPartP_m, itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(nextR));

    for (auto element: elementSet) {
        if (element->getType() == ElementBase::TRAVELINGWAVE) {
            const TravelingWave *TWelement = static_cast<const TravelingWave *>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }

        } else if (element->getType() == ElementBase::RFCAVITY) {
            const RFCavity *RFelement = static_cast<const RFCavity *>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }
        }
    }

}


void ThickTracker::updateReferenceParticle(const BorisPusher &pusher) {
    const double dt = std::min(itsBunch_m->getT(), itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;
    Vector_t Ef(0.0), Bf(0.0);

    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;

    IndexMap::value_t elements = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(RefPartR_m));
    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        CoordinateSystemTrafo refToLocalCSTrafo = itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * referenceToLabCSTrafo_m;

        Vector_t localR = refToLocalCSTrafo.transformTo(RefPartR_m);
        Vector_t localP = refToLocalCSTrafo.rotateTo(RefPartP_m);
        Vector_t localE(0.0), localB(0.0);

        if ((*it)->applyToReferenceParticle(localR,
                                            localP,
                                            itsBunch_m->getT() - 0.5 * dt,
                                            localE,
                                            localB)) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            globalEOL_m = true;
        }

        Ef += refToLocalCSTrafo.rotateFrom(localE);
        Bf += refToLocalCSTrafo.rotateFrom(localB);
    }
    pusher.kick(RefPartR_m, RefPartP_m, Ef, Bf, dt);
    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;
}




void ThickTracker::selectDT() {
    if (itsBunch_m->getIfBeamEmitting()) {
        double dt = itsBunch_m->getEmissionDeltaT();
        itsBunch_m->setdT(dt);
    } else {
        double dt = dtCurrentTrack_m;
        itsBunch_m->setdT(dt);
    }
}

void ThickTracker::changeDT() {
    selectDT();
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}
void ThickTracker::findStartPosition(const BorisPusher &pusher) {

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = dtAllTracks_m.front();
    changeDT();

    if (Util::getEnergy(RefPartP_m, itsBunch_m->getM()) < 1e-3) {
        double gamma = 0.1 / itsBunch_m->getM() + 1.0;
        RefPartP_m = sqrt(std::pow(gamma, 2) - 1) * Vector_t(0, 0, 1);
    }

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        Vector_t oldR = RefPartR_m;
        updateReferenceParticle(pusher);
        pathLength_m += euclidean_norm(RefPartR_m - oldR);

        if (pathLength_m > zStop_m.front()) {
            if (localTrackSteps_m.size() == 0) return;

            dtAllTracks_m.pop();
            localTrackSteps_m.pop();
            zStop_m.pop();

            changeDT();
        }

        double speed = euclidean_norm(RefPartP_m) * Physics::c / sqrt(dot(RefPartP_m, RefPartP_m) + 1);
        if (std::abs(pathLength_m - zstart_m) <=  0.5 * itsBunch_m->getdT() * speed) {
            double tau = (pathLength_m - zstart_m) / speed;

            t += tau;
            itsBunch_m->setT(t);

            RefPartR_m /= (Physics::c * tau);
            pusher.push(RefPartR_m, RefPartP_m, tau);
            RefPartR_m *= (Physics::c * tau);

            pathLength_m = zstart_m;

            CoordinateSystemTrafo update(RefPartR_m,
                                         getQuaternion(RefPartP_m, Vector_t(0, 0, 1)));
            referenceToLabCSTrafo_m = referenceToLabCSTrafo_m * update.inverted();

            RefPartR_m = update.transformTo(RefPartR_m);
            RefPartP_m = update.rotateTo(RefPartP_m);

            return;
        }
    }
}

/**
 * @brief Algorithm for Thick Map-Tracking
 */


void ThickTracker::execute() {

    std::ofstream outfile;
    std::ofstream tmap;
    outfile.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/generatedMaps.txt");
    tmap.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/TransferMap.txt");
    outfile << std::setprecision(20);

	Inform msg("ThickTracker", *gmsg);

	msg << "in execute " << __LINE__ << " " << __FILE__ << endl;

	/*
    First some setup and general preparation. Mostly copied from ParalellTTracker.
    Some of them we maybe do not need at all.
	 */

    OpalData::getInstance()->setInPrepState(true);

    BorisPusher pusher(itsReference);
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    dtCurrentTrack_m = itsBunch_m->getdT();

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
    Options::openMode = Options::APPEND;
    }

    prepareSections();


//  itsOpalBeamline_m.compute3DLattice();
//  itsOpalBeamline_m.save3DLattice();
//  itsOpalBeamline_m.save3DInput();
//
//  std::queue<double> timeStepSizes(dtAllTracks_m);
//  std::queue<unsigned long long> numSteps(localTrackSteps_m);
//  double minTimeStep = timeStepSizes.front();
//  unsigned long long totalNumSteps = 0;
//  while (timeStepSizes.size() > 0) {
//    if (minTimeStep > timeStepSizes.front()) {
//      totalNumSteps = std::ceil(totalNumSteps * minTimeStep / timeStepSizes.front());
//      minTimeStep = timeStepSizes.front();
//    }
//    totalNumSteps += std::ceil(numSteps.front() * timeStepSizes.front() / minTimeStep);
//
//    numSteps.pop();
//    timeStepSizes.pop();
//  }
//
//  itsOpalBeamline_m.activateElements();
//
//  if (OpalData::getInstance()->hasPriorTrack() ||
//      OpalData::getInstance()->inRestartRun()) {
//
//    referenceToLabCSTrafo_m = itsBunch_m->toLabTrafo_m;
//    RefPartR_m = referenceToLabCSTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
//    RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);
//
//    pathLength_m = itsBunch_m->get_sPos();
//    zstart_m = pathLength_m;
//
//    restoreCavityPhases();
//  } else {
//    RefPartR_m = Vector_t(0.0);
//    RefPartP_m = euclidean_norm(itsBunch_m->get_pmean_Distribution()) * Vector_t(0, 0, 1);
//
//    if (itsBunch_m->getTotalNum() > 0) {
//      if (!itsOpalBeamline_m.containsSource()) {
//	RefPartP_m = OpalData::getInstance()->getP0() / itsBunch_m->getM() * Vector_t(0, 0, 1);
//      }
//
//      if (zstart_m > pathLength_m) {
//	findStartPosition(pusher);
//      }
//
//      itsBunch_m->set_sPos(pathLength_m);
//    }
//  }
//
//  Vector_t rmin, rmax;
//  itsBunch_m->get_bounds(rmin, rmax);
//
//  OrbitThreader oth(itsReference,
//		    referenceToLabCSTrafo_m.transformTo(RefPartR_m),
//		    referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
//		    pathLength_m,
//		    -rmin(2),
//		    itsBunch_m->getT(),
//		    minTimeStep,
//		    totalNumSteps,
//		    zStop_m.back() + 2 * rmax(2),
//		    itsOpalBeamline_m);
//
//  oth.execute();
//
  /*
    End of setup and general preparation.
  */


    msg << *itsBunch_m << endl;


    /*
    This is an example how one can loop
    over all elements.
    */

//=======================================================================


    //generate the Hamiltonian
    typedef FTps<double, 2 * DIM> Series;
    typedef FVps<double, 2 * DIM> Map;

    Series x = Series::makeVariable(0);		//SIXVect::X);
    Series px = Series::makeVariable(1);		//SIXVect::PX);
    Series y = Series::makeVariable(2);		//SIXVect::Y);
    Series py = Series::makeVariable(3);		//SIXVect::PY);
    //Series z = Series::makeVariable(4);		//SIXVect::TT);
    Series delta = Series::makeVariable(5);	//SIXVect::PT);



//Diese Parameter muessen per & uebernommen werden
//======================================================

    int order = 2;  //actually already defined
//   define parameter for Hamiltonian
//
//  double b = 1.; 	// field gradient [T/m]
//  double E = 200.; //PartData::getE(); 	//beam energy [MeV]
//  double ds = 0.4;     //stepsize [m]
//
//  double EProt = refPart.getM();  //938.2720813; //MeV/c

    double  gamma0 =  (itsBunch_m->getInitialGamma())   ;    //(EProt + E) / EProt;
    double  beta0 =  (itsBunch_m->getInitialBeta());   //std::sqrt(gamma0 * gamma0 - 1.0) / gamma0;

    double  P0 =  (itsBunch_m-> getP()); //beta0 * gamma0 * EProt * 1e6 / Physics::c;
    double  q =  (itsBunch_m->getQ());	// particle change [e]

//=====================================================
    double E=(itsReference.getE() - itsReference.getM());
    E = std::round(E);
    double AMU = 1.66053873e-27;
    double EZERO = 1.602176462e-19 ;
    double CLIGHT= 2.99792458e8 ;
    double AMUMEV = AMU*(CLIGHT*CLIGHT)/EZERO ;
    double M0= 1.00727646688;


    double ETA = E/(M0*AMUMEV) ;
    double PHI = std::sqrt(ETA*(2+ETA)) ;

    P0 = (AMUMEV*M0)*PHI;
    gamma0=ETA+1;
    beta0= std::sqrt(1-1/(gamma0*gamma0));


//=====================================================


    beta0 = (PHI / (1 + ETA));
    gamma0 = PHI / beta0;

    msg << std::setprecision(20);
    msg << "BreamEnergy: "<< E << endl;
    msg << "P0 COSY:  " << P0 << endl;
    msg << "E0 COSY:  " << AMUMEV*M0 << endl;
//=====================================================
    Series::setGlobalTruncOrder(order);

//  double r0 = 1.; 					// get aperture, resp. the radius of magnet [m]


    ///Creates the Hamiltonian for beam line element
    /**\param element iterative pointer to the elements along the beam line
    * \return Hamiltonian times elemenrLength
    */
    auto Hamiltonian = [&](std::shared_ptr<Component> element) {

        Series H; //Hamiltonian
        Series Hs;//Hamiltonian times ElementLength

        switch(element->getType()) {

            /**Driftspace
             * \f[H_{Drift}= \frac{\delta}{\beta_0} -
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } \f]
             */
            case ElementBase::ElementType::DRIFT: {
                Drift* pDrift= dynamic_cast<Drift*> (element.get());
                //outfile <<  "0T,  ";
                //outfile << pDrift ->getElementLength() << "m" <<std::endl;
                H=( delta / beta0 )
                - sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
                        - ( px*px )
                        - ( py*py )
                        - 1./( beta0 * beta0 * gamma0 * gamma0 ),order
                );
                Hs = H * pDrift->getElementLength();
                break;
            }

            /**Rectangular Bend
             * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
             * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right) \f]
             */
            case ElementBase::ElementType::RBEND: {
                RBend* pRBend= dynamic_cast<RBend*> (element.get());

//		        double rho = P0 / (b * q); 		                        //bending radius [m]
                double h = 1. / pRBend ->getBendRadius();               //inverse bending radius [1/m]
                double K0= pRBend ->getB()*(Physics::c/itsReference.getP());

//		        double phi = 45 * 0.0174532925;                         // tilt angle for dipole fringe field [rad]


                //outfile << K0 << "T,  ";
                //outfile << pRBend ->getArcLength() << "m" <<std::endl;
                H=( delta / beta0 )
                - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order
                        ))
                - (h * x)
                * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order-1
                        ))
                + K0 * x * (1. + 0.5 * h* x);

                Hs = H * pRBend->getArcLength();

                break;
            }

            /**Sector Bend
             * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
             * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right) \f]
             */
            case ElementBase::ElementType::SBEND: {

                SBend* pSBend= dynamic_cast<SBend*> (element.get());
//			    double rho = P0 / (b * q); 		                        //bending radius [m]
                double h = 1. / pSBend ->getBendRadius();               //inverse bending radius [1/m]

                double K0= pSBend ->getB()*(Physics::c/itsReference.getP());
                msg<< "Magnetic Field Sbend:" << pSBend ->getB() << "T" << endl;
//			    double phi = 45 * 0.0174532925;                         // tilt angle for dipole fringe field [rad]
                msg<< "Effecgtive Length:" << pSBend ->getEffectiveLength() << "m" << endl;
                //outfile << pSBend ->getB() << "T,  ";
                //outfile << pSBend ->getEffectiveLength() << "m" <<std::endl;
                H=( delta / beta0 )
                - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order
                        ))
                - (h * x)
                * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order-1
                        ))
                + K0 * x * (1. + 0.5 * h* x);

                Hs= H * pSBend->getEffectiveLength();

                break;
            }

            /**Quadrupole "in Multipolegroup"
             * \f[H_{Quadrupole}= \frac{\delta}{\beta_0} -
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2  -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }	+
             * \frac{1}{2} k_1 \left( x^2 - y^2 \right) \f]
             */
            case ElementBase::ElementType::MULTIPOLE: {
                Multipole* pMultipole= dynamic_cast<Multipole*> (element.get());

                double K1= pMultipole->getField().getNormalComponent(2)*(Physics::c/P0);
                K1= std::round(K1*1e6)/1e6;
//                outfile << K1 << "T/m,  ";
//                outfile << pMultipole ->getElementLength() << "m,  " ;
//                outfile << P0 << "MeV"<<std::endl;
//		        K1 = (q * b) / (P0 * r0)

                H= ( delta / beta0 )
                - sqrt ((1./ beta0 + delta ) *(1./ beta0 + delta)
                        - ( px*px )
                        - ( py*py )
                        - 1./( beta0*beta0 * gamma0*gamma0 ),order
                )
                + 0.5 * (q * K1) *(Physics::c/P0) * (x*x - y*y);

                Hs=H* pMultipole ->getElementLength();

                break;
            }

            default:
            throw LogicalError("ThickTracker::exercute,",
                    "Please use already defined beam line element:  At this time just driftspace, multipoles and dipoles");
            break;

        }
        return Hs;
    };


//=======================================================================



    FieldList allElements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);
    struct sort_by_pos {
    bool operator()(const ClassicField &a, const ClassicField &b)
        {return a.getElement()->getElementPosition() < b.getElement()->getElementPosition();}
    };

    allElements.sort(sort_by_pos());



    FieldList::iterator it = allElements.begin();
    const FieldList::iterator end = allElements.end();
    if (it == end) msg << "No element in lattice" << endl;

    Map combinedMaps;





//    msg << "beam parameter " << itsReference.getE() - itsReference.getM() << endl;
//    msg << "beam line length: " << allElements.size() << endl;
//    //const int el = allElements.size();
    msg << std::setprecision(20);
    msg << "P0 OPAL:  " << itsBunch_m-> getP() << endl;
    msg << "E0 OPAL:  " << itsBunch_m-> getM() << endl;



//    msg << "q:  " << q << endl;


    std::map<std::string, double> elementDic;
    std::map<std::string, double>::iterator dicit;

    std::vector<Map> mapVec;
    double selement =0;
    for (; it != end; ++ it) {
        dicit = elementDic.begin();

        std::shared_ptr<Component> element = (*it).getElement();


        dicit= elementDic.find(element->getName());
        if (dicit !=elementDic.end()){
            throw LogicalError("ThickTracker::execute,",
                                "Same Element twice in beamline:"+ element->getName());

        }else{
            elementDic.insert(std::pair<std::string, double> (element->getName(), element->getElementPosition()));
        }

        if (selement > element->getElementPosition()){
            msg << "There is an overlap! @ Element: " << element ->getName() <<
                    "-> starts at: " << element->getElementPosition()<<
                    " the previous ends at: " << selement << endl;
//            throw LogicalError("ThickTracker::exercute,",
//                                            "Overlap at element:"+ element->getName());

        }else{
            selement=element->getElementPosition()+ element ->getElementLength();
        }

        if (element ->getType()==ElementBase::MONITOR){
            mapVec.insert(mapVec.end(), combinedMaps);
            combinedMaps.identity();

        }


    //    outfile << "Element name: " << element->getName() << std::endl;
    //    outfile << "Element type:  " << element->getType()<< std::endl;
        msg << "=============================="<< endl;
        msg << "Element name: " << element->getName() << endl;
        msg << "Element Length: " << element->getElementLength() << endl;
//        msg << "Element nummer:   " << element->getType() << endl;
        msg <<  "ElementPosition"<< element->getElementPosition()<< endl;
        msg << "EntryPoint" << element->getEntrance()<< endl;


//        outfile << "Element Hamiltonian: \n" << Hamiltonian(element) << std::endl;
//        outfile << element->getName() << ",   ";
//        outfile << E/1.0e6 << "MeV,  ";

        Map ElementMap=ExpMap(-Hamiltonian(element) ,order);
//        outfile << Hamiltonian(element);
//        outfile << "\n=========="<< std::endl;
        //outfile << ElementMap << std::endl;
        combinedMaps = ElementMap * combinedMaps;
        //outfile << "CombinedMap:  \n" << combinedMaps << std::endl;
        //outfile << "------------------------------------------------------------" << std::endl;


  }
  tmap << "A customized FODO" << std::endl;
  tmap << combinedMaps << std::endl;
//  outfile << combinedMaps << std::endl;

  FVector<double, 6> particle;
  OpalParticle part;
  outfile <<"P0:  " << itsBunch_m->getP()/(M0*AMUMEV) << "    "
		  << itsBunch_m->getInitialGamma()*itsBunch_m->getInitialBeta()<< std::endl;


  for (unsigned int idx=0; idx< itsBunch_m->getTotalNum(); idx++){
	  part =itsBunch_m ->get_part(idx);

	  for(int dim=0; dim<6; dim++){
		  particle[dim]=itsBunch_m ->get_part(idx)[dim];
	  }
	  outfile << "=======================" << std::endl;
	  outfile << "Particle#: " << idx << "   "<<itsBunch_m->getBeta(idx) << std::endl;
	  outfile << particle;
//	  outfile << "particle[5]:   " << particle[5] << std::endl;
//	  outfile << "1st term=   " << (particle[5]*M0*AMUMEV/ itsBunch_m->getP()) << std::endl;
//	  outfile << "1st term=   " << (particle[5]*itsBunch_m->getM()/ itsBunch_m->getP()) << std::endl;
	  particle[5]=
				(particle[5]*itsBunch_m->getM()/itsBunch_m->getP())
				* std::sqrt( 1./(particle[5]* particle[5]) +1)
				-1./itsBunch_m->getInitialBeta();

//	  outfile << "---> Unit Transformation <--- \n" << particle;
	  outfile << "---> The Map is applied <---"<< std::endl;
	  particle= combinedMaps * particle;
//	  outfile << "particle[5]" << particle[5] << std::endl;
//
//	  outfile << "divident:  " << particle[5]*M0*AMUMEV << std::endl;
//	  outfile << "divisor:  " << itsBunch_m->getP() << std::endl;
//	  outfile << "M0:   "  << M0*AMUMEV << "OPAL ->  "<< itsBunch_m->getM() << std::endl;
//	  outfile << "P0:   " << itsBunch_m->getP() << std::endl;
//	  outfile << "gamma0:   "  << itsBunch_m->getInitialGamma() << std::endl;
//	  outfile << "delta:  " << particle[5] << std::endl;
	  particle[4]= selement;
	  particle[5]= (particle[5] + 1./itsBunch_m->getInitialBeta()) * itsBunch_m->getP()/itsBunch_m->getM()
			  /std::sqrt( 1./(itsBunch_m ->get_part(idx)[5]* itsBunch_m ->get_part(idx)[5]) +1) ;

	  outfile << particle << std::endl;
  }




    msg << part.y() << endl;







	msg<< "number of Particles" << itsBunch_m->getTotalNum() << endl;


	for(int i=0; i<6; i++){
			msg << "OPALPart" << part[i]<< endl;
	    	msg << "OPALBart" << particle[i]<< endl;
	    }

	std::map <std::string, double>::iterator mi;
	  for (mi=elementDic.begin(); mi!=elementDic.end(); ++mi)
	      msg << mi->first << " = " << mi->second << endl;

}


/*

namespace {
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx) {
        //std::cerr << "==> In implicitInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
        // Default: nx = 20, cx = 4

        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
        // version of implicitInt2 to obtain zf accurate through fourth-order in the
        // step-size ds.  When f derives from a Hamiltonian---i.e., f = J.grad(H)---
        // then this routine performs symplectic integration.  The optional arguments
        // nx and cx have the same meaning as in implicitInt2().

        // Convergence warning flag.
        static bool cnvWarn = false;

        // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
        static const double yt = pow(2., 1 / 3.);
        static const double ya = 1 / (2. - yt);
        static const double yb = -yt * ya;

        // Build matrix grad(f).
        MxSeries gradf;
        for(int i = 0; i < PSdim; ++i)
            for(int j = 0; j < PSdim; ++j)
                gradf[i][j] = f[i].derivative(j);

        // Initialize accumulated length, current step-size, and number of cuts.
        double as = std::abs(s), st = 0., dsc = std::abs(ds);
        if(s < 0.) dsc = -dsc;
        int ci = 0;

        // Integrate each step.
        Vector zf = zin;
        while(std::abs(st) < as) {
            Vector zt;
            bool ok = true;
            try {
                if(std::abs(st + dsc) > as) dsc = s - st;
                zt = ::implicitIntStep(zf, f, gradf, ya * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, yb * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, ya * dsc, nx);
            } catch(ConvergenceError &cnverr) {
                if(++ci > cx) {
                    std::string msg = "Convergence not achieved within " + NumToStr<int>(cx) + " cuts of step-size!";
                    throw ConvergenceError("ThickTracker::implicitInt4()", msg);
                }
                if(!cnvWarn) {
                    std::cerr << " <***WARNING***> [ThickTracker::implicitInt4()]:\n"
                              << "   Cutting step size, a probable violation of the symplectic condition."
                              << std::endl;
                    cnvWarn = true;
                }
                dsc *= 0.5;
                ok = false;
            }
            if(ok) {zf = zt; st += dsc;}
        }

        //std::cerr << "==> Leaving implicitInt4(...)" << std::endl;
        return zf;
    }

    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx) {
        //std::cerr << "==> In implicitIntStep(zin,f,gradf,ds,nx) ..." << std::endl;
        //std::cerr << "  ds = " << ds << std::endl;
        //std::cerr << " zin =\n" << zin << std::endl;
        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a single step of size ds, using Newton's method to solve the
        // implicit equation zf = zin + ds*f((zin+zf)/2).  For reasons of efficiency,
        // its arguments include the matrix gradf = grad(f).  The (optional) argument
        // nx limits the number of Newton iterations.  This routine returns a result
        // zf accurate through second-order in the step-size ds.  When f derives from
        // a Hamiltonian---i.e., f=J.grad(H)---then this routine performs symplectic
        // integration.

        // Set up flags, etc., for convergence (bounce) test.
        FVector<bool, PSdim> bounce(false);
        Vector dz, dz_old;
        int bcount = 0;
        static const double thresh = 1.e-8;

        // Use second-order Runge-Kutta integration to determine a good initial guess.
        double ds2 = 0.5 * ds;
        Vector z = f.constantTerm(zin);
        z = zin + ds2 * (z + f.constantTerm(zin + ds * z));

        // Newton iterations:
        //   z :-> [I-ds/2.grad(f)]^{-1}.[zin+ds.f((zin+z)/2)-ds/2.grad(f).z]
        // (A possible method for speeding up this computation would
        //  be to recompute grad(f) every n-th step, where n > 1!)
        Vector zf;
        int ni = 0;
        while(bcount < PSdim) {
            if(ni == nx) {
                std::string msg = "Convergence not achieved within " + NumToStr<int>(nx) + " iterations!";
                throw ConvergenceError("ThickTracker::implicitIntStep()", msg);
            }

            // Build gf = -ds/2.grad(f)[(zin+z)/2] and idgf_inv = [I-ds/2.grad(f)]^{-1}[(zin+z)/2].
            Vector zt = 0.5 * (zin + z);
            Matrix gf, idgf, idgf_inv;
            for(int i = 0; i < PSdim; ++i)
                for(int j = 0; j < PSdim; ++j)
                    gf[i][j] = -ds2 * gradf[i][j].evaluate(zt);
            idgf = gf;
            for(int i = 0; i < PSdim; ++i) idgf[i][i] += 1.;
            FLUMatrix<double, PSdim> lu(idgf);
            idgf_inv = lu.inverse();

            // Execute Newton step.
            zf = idgf_inv * (zin + ds * f.constantTerm(zt) + gf * z);

            //std::cerr << " -(ds/2)grad(f) =\n" << gf << std::endl;
            //std::cerr << " f =\n" << f.constantTerm(zt) << std::endl;
            //std::cerr << "zk =\n" << zf << std::endl;

            // Test for convergence ("bounce" test).
            dz_old = dz;
            dz = zf - z;
            if(ni) { // (we need at least two iterations before testing makes sense)
                for(int i = 0; i < PSdim; ++i) {
                    if(!bounce[i] && (dz[i] == 0. || (std::abs(dz[i]) < thresh && std::abs(dz[i]) >= std::abs(dz_old[i]))))
                        {bounce[i] = true; ++bcount;}
                }
            }
            z = zf;
            ++ni;
        }

        //std::cerr << "  zf =\n" << zf << std::endl;
        //std::cerr << "==> Leaving implicitIntStep(zin,f,gradf,ds,nx)" << std::endl;
        return zf;
    }

    Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx) {
        //std::cerr << "==> In fixedPointInt2(zin,f,ds,nx) ..." << std::endl;
        // Default: nx = 50

        //std::cerr << "  ds = " << ds << std::endl;
        //std::cerr << " zin =\n" << zin << std::endl;
        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a single step of size ds by iterating the equation
        //         z = zin + ds * f((zin+z)/2)
        // to find a fixed-point zf for z.  It is accurate through second order in the
        // step size ds.

        // Set up flags, etc., for convergence (bounce) test.
        FVector<bool, PSdim> bounce(false);
        Vector dz, dz_old;
        int bcount = 0;
        static const double thresh = 1.e-8;

        // Iterate z :-> zin + ds * f( (zin + z)/2 ).
        Vector zf;
        Vector z = zin;
        int ni = 0;
        while(bcount < PSdim) {
            if(ni == nx) {
                std::string msg = "Convergence not achieved within " + NumToStr<int>(nx) + " iterations!";
                throw ConvergenceError("ThickTracker::fixedPointInt2()", msg);
            }

            // Do iteration.
            zf = zin + ds * f.constantTerm((zin + z) / 2.0);

            // Test for convergence.
            dz_old = dz;
            dz = zf - z;
            if(ni) { // (we need at least two iterations before testing makes sense)
                for(int i = 0; i < PSdim; ++i) {
                    if(!bounce[i] && (dz[i] == 0. || (std::abs(dz[i]) < thresh && std::abs(dz[i]) >= std::abs(dz_old[i]))))
                        {bounce[i] = true; ++bcount;}
                }
            }
            z = zf;
            ++ni;
        }
        //std::cerr << "  zf =\n" << zf << std::endl;
        //std::cerr << "==> Leaving fixedPointInt2(...)" << std::endl;
        return zf;
    }

    Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx) {
        //std::cerr << "==> In fixedPointInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
        // Default: nx = 50, cx = 4

        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
        // version of fixedPointInt2 to obtain zf accurate through fourth-order in the
        // step-size ds.  The optional arguments nx and cx have the same meaning as in
        // implicitInt2().

        // Convergence warning flag.
        static bool cnvWarn = false;

        // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
        static const double yt = pow(2., 1 / 3.);
        static const double ya = 1 / (2. - yt);
        static const double yb = -yt * ya;

        // Initialize accumulated length, current step-size, and number of cuts.
        double as = std::abs(s), st = 0., dsc = std::abs(ds);
        if(s < 0.) dsc = -dsc;
        int ci = 0;

        // Integrate each step.
        Vector zf = zin;
        while(std::abs(st) < as) {
            Vector zt;
            bool ok = true;
            try {
                if(std::abs(st + dsc) > as) dsc = s - st;
                zt = ::fixedPointInt2(zf, f, ya * dsc, nx);
                zt = ::fixedPointInt2(zt, f, yb * dsc, nx);
                zt = ::fixedPointInt2(zt, f, ya * dsc, nx);
            } catch(ConvergenceError &cnverr) {
                if(++ci > cx) {
                    std::string msg = "Convergence not achieved within " + NumToStr<int>(cx) + " cuts of step-size!";
                    throw ConvergenceError("ThickTracker::fixedPointInt4()", msg);
                }
                if(!cnvWarn) {
                    std::cerr << " <***WARNING***> [ThickTracker::fixedPointInt4()]:\n"
                              << "   Cutting step size, a probable violation of the symplectic condition."
                              << std::endl;
                    cnvWarn = true;
                }
                dsc *= 0.5;
                ok = false;
            }
            if(ok) {zf = zt; st += dsc;}
        }

        //std::cerr << "==> Leaving fixedPointInt4(...)" << std::endl;
        return zf;
    }
};
*/
