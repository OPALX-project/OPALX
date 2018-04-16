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
#include "Classic/FixedAlgebra/FDoubleEigen.h"


#include "Classic/Fields/BSingleMultipoleField.h"


#include "Classic/Algorithms/PartData.h"  //for the beam reference

#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"

#include "Physics/Physics.h"

#include "Elements/OpalBeamline.h"
#include <Classic/BeamlineCore/MultipoleRep.h>
#include "Distribution/MapGenerator.h"


#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#define DIM 3

class Beamline;
class PartData;
using Physics::c;


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
  localTrackSteps_m(),
  truncOrder_m(1) // linear
{
		 x=(series_t::makeVariable(0));
		 px = series_t::makeVariable(1);		//SIXVect::PX);
		 y = series_t::makeVariable(2);			//SIXVect::Y);
		 py = series_t::makeVariable(3);		//SIXVect::PY);
		//z = Series::makeVariable(4);		//SIXVect::TT);
		 delta = series_t::makeVariable(5);		//SIXVect::PT);
}


ThickTracker::ThickTracker(const Beamline &beamline,
               PartBunchBase<double, 3> *bunch,
			   DataSink &ds,
			   const PartData &reference,	  
               bool revBeam, bool revTrack,
			   const std::vector<unsigned long long> &maxSteps,
			   double zstart,
			   const std::vector<double> &zstop,
			   const std::vector<double> &dt,
               const int& truncOrder):
	  Tracker(beamline, bunch, reference, revBeam, revTrack),
  	  itsDataSink_m(&ds),
  	  itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getCoordTransformationTo()),
  	  pathLength_m(0.0),
  	  zstart_m(zstart),
  	  zStop_m(),
  	  dtCurrentTrack_m(0.0),
  	  dtAllTracks_m(),
  	  localTrackSteps_m(),
  	  truncOrder_m(truncOrder)
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
   x=(series_t::makeVariable(0));
   px = series_t::makeVariable(1);		//SIXVect::PX);
   y = series_t::makeVariable(2);			//SIXVect::Y);
   py = series_t::makeVariable(3);		//SIXVect::PY);
  //Series z = Series::makeVariable(4);		//SIXVect::TT);
   delta = series_t::makeVariable(5);		//SIXVect::PT);
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

/**Drift Space Hamiltonian
* \f[H_{Drift}= \frac{\delta}{\beta_0} -
* \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } \f]
*/
void ThickTracker::setHamiltonianDrift(series_t& H, double& beta0, double& gamma0){

			H=( delta / beta0 )
			- sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
					- ( px*px )
					- ( py*py )
					- 1./( beta0 * beta0 * gamma0 * gamma0 ),truncOrder_m+1
			);
}
/**Rectangular Bend Hamiltonian
 * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
 * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
 * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right) \f]
 */
void ThickTracker::setHamiltonianRBend(series_t& H, double& beta0, double& gamma0, double& q,  double& h, double& K0 ){



	H=( delta / beta0 )
	- (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
					- ( px*px )
					- ( py*py )
					- 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m+1
			))
	- (h * x)
	* (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
					- ( px*px )
					- ( py*py )
					- 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m
			))
	+ K0 * x * (1. + 0.5 * h* x);


}


/**Sector Bend Hamiltonian
 * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
 * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
 * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right) \f]
 */
void ThickTracker::setHamiltonianSBend(series_t& H, double& beta0, double& gamma0, double& q,  double& h, double& K0 ){


	H=( delta / beta0 )
                    - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                    - ( px*px )
                                    - ( py*py )
                                    - 1./( beta0*beta0 * gamma0*gamma0 ),(truncOrder_m+1)
                            ))
                    - (h * x)
                    * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                    - ( px*px )
                                    - ( py*py )
                                    - 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m
                            ))
                    + K0 * x * (1. + 0.5 * h* x);



}


/**Quadrupole "in Multipolegroup" Hamiltonian
 * \f[H_{Quadrupole}= \frac{\delta}{\beta_0} -
 * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2  -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }	+
 * \frac{1}{2} k_1 \left( x^2 - y^2 \right) \f]
 */
void ThickTracker::setHamiltonianQuadrupole(series_t& H, double& beta0, double& gamma0, double& q, double& K1 ){

	H= ( delta / beta0 )
	- sqrt ((1./ beta0 + delta ) *(1./ beta0 + delta)
			- ( px*px )
			- ( py*py )
			- 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m+1
	)
	+ 0.5 * K1 * (x*x - y*y);
}


///Fills undefined beam path with a Drift Space
/** \f[H_{Drift}= \frac{\delta}{\beta_0} -
* \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } \f]
*/
void ThickTracker::fillDrift(std::list<structMapTracking>& mapBeamLine,double& elementPos, double& undefSpace){

    Inform msg("ThickTracker", *gmsg);
    msg <<"filled in a drift in between:  " << elementPos-undefSpace << " and "
                            << elementPos << " length: "<< undefSpace << endl;

	series_t H;
	structMapTracking fillDrift;

	double  gamma0 =  (itsBunch_m->getInitialGamma())   ;   //(EProt + E) / EProt;
	double  beta0 =  (itsBunch_m->getInitialBeta());   		//std::sqrt(gamma0 * gamma0 - 1.0) / gamma0;

//=====================================================
//Redefine constants for comparison with COSY Infinity
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

	gamma0=ETA+1;
	beta0= std::sqrt(1-1/(gamma0*gamma0));

	beta0 = (PHI / (1 + ETA));
	gamma0 = PHI / beta0;

//=====================================================


	fillDrift.elementName="fillDrift";
	fillDrift.elementPos= elementPos;
	fillDrift.nSlices=1;
	fillDrift.stepSize=undefSpace;
	setHamiltonianDrift(H, beta0, gamma0);

	fillDrift.elementMap=ExpMap(-H * fillDrift.stepSize, truncOrder_m);


	mapBeamLine.push_back(fillDrift);
}



///Creates the Hamiltonian for beam line element
/** \param element iterative pointer to the elements along the beam line
*
*/

void ThickTracker::defMapTrackingElement(std::shared_ptr<Component> element, structMapTracking& elSrct, std::list<structMapTracking>& mBL){
    std::ofstream mbl;
    mbl.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/mapBeamLine.txt",std::ios::app);

	series_t H;
	double  gamma0 =  (itsBunch_m->getInitialGamma())   ;   //(EProt + E) / EProt;
	double  beta0 =  (itsBunch_m->getInitialBeta());   		//std::sqrt(gamma0 * gamma0 - 1.0) / gamma0;

	double  P0 =  (itsBunch_m-> getP()); 					//beta0 * gamma0 * EProt * 1e6 / Physics::c;
	double  q =  (itsBunch_m->getQ());						// particle change [e]

//=====================================================
//Redefine constants for comparison with COSY Infinity
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

	beta0 = (PHI / (1 + ETA));
	gamma0 = PHI / beta0;

//=====================================================


	//select the right Hamiltonian for the beam line
	switch(element->getType()) {
		case ElementBase::ElementType::DRIFT: {
			Drift* pDrift= dynamic_cast<Drift*> (element.get());

			elSrct.nSlices= pDrift->getNSlices();
			elSrct.stepSize= pDrift->getElementLength()/elSrct.nSlices;

			setHamiltonianDrift(H, beta0, gamma0);
			elSrct.elementMap=ExpMap(-H*elSrct.stepSize ,truncOrder_m);
			mBL.push_back(elSrct);
		break;
		}


		case ElementBase::ElementType::RBEND: {
			RBend* pRBend= dynamic_cast<RBend*> (element.get());
			elSrct.nSlices= pRBend->getNSlices();
			elSrct.stepSize= pRBend->getElementLength()/elSrct.nSlices;

			double h = 1. / pRBend ->getBendRadius();               //inverse bending radius [1/m]
			double K0= pRBend ->getB()*(Physics::c/itsReference.getP());
			K0=std::round(K0*1e6)/1e6 *q*(Physics::c/P0);

			setHamiltonianRBend(H, beta0, gamma0, q, h, K0);

			H= H /pRBend->getElementLength() *pRBend->getArcLength();
			elSrct.elementMap=ExpMap(-H*elSrct.stepSize ,truncOrder_m);

			mBL.push_back(elSrct);
        break;
		}


		case ElementBase::ElementType::SBEND: {
			SBend* pSBend= dynamic_cast<SBend*> (element.get());

			elSrct.nSlices= pSBend->getNSlices();
			if (pSBend->getEntranceAngle()>1e-6){

			}
			if (pSBend->getExitAngle()>1e-6){

			}

			//elSrct.stepSize= pSBend->getElementLength()/elSrct.nSlices;
			elSrct.stepSize= pSBend->getEffectiveLength()/elSrct.nSlices;

			double K0= pSBend ->getB()*(Physics::c/itsReference.getP());

			K0= std::round(K0*1e6)/1e6 *q*(Physics::c/P0);


			double h = 1./ pSBend ->getBendRadius();               //inverse bending radius [1/m]
			K0<0 ? h = -1*h : h = 1*h;

			mbl<< "design Energy="<< pSBend->getDesignEnergy()<< std::endl;
			mbl<< "design Radius="<< pSBend->getBendRadius() << std::endl;
			mbl << "K0:   "<< K0 << std::endl;
			mbl << "P0:   " << P0/Physics::c << std::endl;
			mbl << "h:   "<< h << std::endl;
			mbl << "entrAngle" << pSBend->getEntranceAngle() << std::endl;


			setHamiltonianSBend(H, beta0, gamma0, q, h, K0);

			mbl << "Hamiltonian"<< -H*elSrct.stepSize <<std::endl;


			Inform msg("ThickTracker", *gmsg);
			msg << "ChordLength: (->getElementLength()) "<< pSBend->getElementLength()<< "m"<< endl;
			msg << "ArcLength: (->getEffectiveLength()) "<< pSBend->getEffectiveLength()<< "m"<< endl;

			//H= H /pSBend->getElementLength() *pSBend->getEffectiveLength();

			elSrct.elementMap=ExpMap(-H*elSrct.stepSize, truncOrder_m);

			//Consider Entrance and Exit pole face rotation
			if (pSBend->getEntranceAngle()>1e-6){
			    structMapTracking entrAngle;
			    entrAngle.elementName=pSBend->getName()+"entrAngle";
			    entrAngle.elementPos=(pSBend->getElementPosition());
			    entrAngle.nSlices=1;

			    H=(x*x - y*y) *0.5 * K0 * std::tan(pSBend->getEntranceAngle());
			    entrAngle.elementMap=ExpMap(-H,truncOrder_m);
			    mBL.push_back(entrAngle);

			}

			mBL.push_back(elSrct);

			if (pSBend->getExitAngle()>1e-6){
			    structMapTracking extAngle;
			    extAngle.elementName=pSBend->getName()+"extAngle";
			    extAngle.elementPos=(pSBend->getElementPosition()+ pSBend->getElementLength());
			    extAngle.nSlices=1;


                H=(x*x - y*y) *0.5 * K0 * std::tan(pSBend->getExitAngle())*q;
                extAngle.elementMap=ExpMap(-H,truncOrder_m);
                mBL.push_back(extAngle);
			}
			break;
		}


		case ElementBase::ElementType::MULTIPOLE: {
			Multipole* pMultipole= dynamic_cast<Multipole*> (element.get());

			elSrct.nSlices= pMultipole->getNSlices();
			elSrct.stepSize= pMultipole->getElementLength()/elSrct.nSlices;

			double K1= pMultipole->getField().getNormalComponent(2)*(Physics::c/P0);


			K1= std::round(K1*1e6)/1e6 *q*(Physics::c/P0);
			Inform msg("ThickTracker", *gmsg);
            msg << "K1 Quad:  "<< K1<< endl;


			setHamiltonianQuadrupole(H, beta0, gamma0, q, K1);
			elSrct.elementMap=ExpMap(-H*elSrct.stepSize ,truncOrder_m);

			mBL.push_back(elSrct);

			break;
		}
		case ElementBase::ElementType::MONITOR: {
		    Inform msg("ThickTracker", *gmsg);
		    msg << "Skip Monitor: " << element->getName() << endl;
		    elSrct.stepSize=0;
		    elSrct.nSlices=0;
		    elSrct.elementMap.identity();

		    mBL.push_back(elSrct);





		    break;
		}

		default:{
			throw LogicalError("ThickTracker::exercute,",
					"Please use already defined beam line element:  At this time just driftspace, multipoles and dipoles");
			break;
		}

	}
}


/**
 * @brief Algorithm for Thick Map-Tracking
 */


void ThickTracker::execute() {

    Inform msg("ThickTracker", *gmsg);

    msg << "in execute " << __LINE__ << " " << __FILE__ << endl;


    OpalData::getInstance()->setInPrepState(true);

    BorisPusher pusher(itsReference);
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    dtCurrentTrack_m = itsBunch_m->getdT();

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
    Options::openMode = Options::APPEND;
    }

    prepareSections();

    msg << *itsBunch_m << endl;
    msg << std::setprecision(10);


    
    msg << "Tuncation order: " << this->truncOrder_m << endl;

    series_t::setGlobalTruncOrder(truncOrder_m+1);

    map_t combinedMaps; //Final Transfer map_t
    series_t H;         //createHamiltonian
    OpalParticle part;

    structMapTracking mapTrackingElement;

    std::list<structMapTracking> mapBeamLine;
    std::list<structMapTracking>::iterator mapBeamLineit;

    std::vector<map_t> mapVec; // not necessary

    double positionMapTracking = zstart_m;

    itsBunch_m->set_sPos(positionMapTracking);


    //files for analysis

#ifdef PHIL_WRITE

    msg<<"P: " <<itsBunch_m->getP()<< endl;
    msg<<"beta: " << itsBunch_m->getInitialBeta()<< endl;
    msg<<"gamma: " << itsBunch_m->getInitialGamma()<< endl;
    msg<<"E0:  " << itsBunch_m->getM() <<endl;
    std::ofstream outfile;
    outfile.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/generatedMaps.txt");
    outfile << std::setprecision(8);

    std::ofstream mbl;
    mbl.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/mapBeamLine.txt",std::ios::app);
    mbl << std::setprecision(16);

    std::ofstream tmap;
    tmap.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/TransferMap.txt");
    tmap << std::setprecision(16);


    std::ofstream twiss;
    twiss.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/twiss.txt");
    tmap << std::setprecision(16);
#endif
    
    FieldList allElements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);

	//sorts beamline according elementposition
	struct sort_by_pos {
		bool operator()(const ClassicField &a, const ClassicField &b)
			{return a.getElement()->getElementPosition() < b.getElement()->getElementPosition();}
	};
	allElements.sort(sort_by_pos());

    FieldList::iterator it = allElements.begin();
	const FieldList::iterator end = allElements.end();
	if (it == end) msg << "No element in lattice" << endl;


    //loop over beam line
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();

        mapTrackingElement.elementPos= std::round(element->getElementPosition()*1e6)/1e6;
        mapTrackingElement.elementName= element->getName();

        //check for double implementations
        for(mapBeamLineit=mapBeamLine.begin(); mapBeamLineit != mapBeamLine.end(); ++mapBeamLineit){
        	if(mapBeamLineit->elementName == mapTrackingElement.elementName){
        		throw LogicalError("ThickTracker::execute,",
        		                                "Same Element twice in beamline:"+ element->getName());
        	}
        }

        // Fill Drift , if necessary
        if (positionMapTracking < mapTrackingElement.elementPos -1e-6){
            double undefSpace=mapTrackingElement.elementPos - positionMapTracking;
            fillDrift(mapBeamLine, positionMapTracking, undefSpace);
        }

        if (element->getType() == ElementBase::ElementType::MONITOR){
            element->setElementLength(0.0);
        }

        double elementLength=element->getElementLength();

        if (element->getType() == ElementBase::ElementType::SBEND){
            SBend* pSBend= dynamic_cast<SBend*> (element.get());
            elementLength= pSBend->getEffectiveLength();
        }


        //check for overlap //TODO fix the bug with overlapping monitors!!!
        if ( positionMapTracking > mapTrackingElement.elementPos+ 1e-6){
            msg << "There is an overlap! @ Element: " << mapTrackingElement.elementName <<
                   "-> starts at: " << mapTrackingElement.elementPos<<
                   "m, overlap length: " << positionMapTracking-mapTrackingElement.elementPos<<"m" << endl;

//            throw LogicalError("ThickTracker::exercute,",
//                                            "Overlap at element:"+ element->getName());

            if(positionMapTracking < mapTrackingElement.elementPos+ elementLength){
                positionMapTracking=std::round((mapTrackingElement.elementPos + elementLength)*1e6)/1e6;
            }

        }else{
            positionMapTracking=std::round((mapTrackingElement.elementPos + elementLength)*1e6)/1e6; //<--DEBUG
        }




//        //combine maps in between Monitors
//
//        if (element ->getType()==ElementBase::MONITOR){
//            mapVec.insert(mapVec.end(), combinedMaps);
//            combinedMaps.identity();
//
//        }


        //Create mapTrackElement and insert it in mapBemLine
        defMapTrackingElement(element, mapTrackingElement, mapBeamLine);


#ifdef PHIL_WRITE
        mbl<<   "Name:  " << mapTrackingElement.elementName<<std::endl
                <<"Position:  " <<mapTrackingElement.elementPos <<std::endl
                <<"StepSize:  " <<mapTrackingElement.stepSize <<std::endl
                <<"NSLices:  " <<mapTrackingElement.nSlices<<std::endl
                <<"map:\n" <<mapTrackingElement.elementMap<<std::endl
                <<"---------------------------------" <<std::endl;
    //outfile << element->getName() << std::endl;
    //outfile << createHamiltonian(element, mapTrackingElement.stepSize, mapTrackingElement.nSlices);

#endif
    }
    //=================================
    // TODO: remove Messages later
    //=================================
    for (mapBeamLineit=mapBeamLine.begin(); mapBeamLineit != mapBeamLine.end(); ++mapBeamLineit) {
        msg << "=============================="<< endl
            << "Name: " <<  mapBeamLineit->elementName << endl
            << "InPosition:  "<< mapBeamLineit->elementPos <<endl
            << "FinPosition: "<< mapBeamLineit->elementPos + mapBeamLineit->nSlices*mapBeamLineit->stepSize << endl;

        mbl<<   "Name:  " << mapBeamLineit->elementName<<std::endl
           << "Position:  " <<mapBeamLineit->elementPos <<std::endl
           << "StepSize:  " <<mapBeamLineit->stepSize <<std::endl
           << "NSLices:  " <<mapBeamLineit->nSlices<<std::endl
           << "map:\n" <<mapBeamLineit->elementMap<<std::endl
           << "---------------------------------" <<std::endl;



    }
    //=================================

    std::size_t totalSlices=0;
    fMatrix_t tFMatrix;
    double position=0;
    //combined map_t
    for (mapBeamLineit=mapBeamLine.begin(); mapBeamLineit != mapBeamLine.end(); ++mapBeamLineit) {

        for (std::size_t slice=0; slice < mapBeamLineit->nSlices; slice++){
            combinedMaps= mapBeamLineit->elementMap * combinedMaps;
            combinedMaps=combinedMaps.truncate(truncOrder_m);
            tFMatrix= combinedMaps.linearTerms();
            linTAnalyze(tFMatrix);
            position+= mapBeamLineit->stepSize;
            //twiss << position << std::endl;
            totalSlices++;
        }
    }

#ifdef PHIL_WRITE
    //twiss << "=================================================" << std::endl;
    outfile <<"Total Particle Number:  " <<  itsBunch_m->getTotalNum() << "  Total number of Slices:   "<< totalSlices << std::endl;

    //eliminate higher order terms (occurring through multiplication)




	tmap << "Transfermap" << std::endl;
    tmap << "Gantry2" << std::endl;
    tmap << combinedMaps << std::endl;
#endif
    //track the Particles
	setTime();

	double t = itsBunch_m->getT();
	itsBunch_m->setT(t);
	OpalData::getInstance()->setInPrepState(false);
	selectDT();

//    trackParticles_m(
//#ifdef PHIL_WRITE
//        outfile,
//#endif
//        mapBeamLine);
    
//-------------------------



    //fMatrix_t sFMatrix=  itsBunch_m->getSigmaMatrix();
//    fMatrix_t tFMatrix= combinedMaps.linearTerms();
//    linTAnalyze(tFMatrix);

    /*
    
	FMatrix<double, 2 * DIM, 2 * DIM> sigmaSFMatrix= sFMatrix*skewMatrix;

	cfMatrix_t eigenValM, eigenVecM, invEigenValM;*/


    linSigAnalyze();


#ifdef PHIL_WRITE

    /*tmap<< "\n\n--------------------------"<< std::endl;
    tmap<< "the S*SigmaMatrix -> Bunch"<< std::endl;
    tmap<< "--------------------------"<< std::endl;
    tmap<< sigmaSFMatrix<<std::endl;
    tmap<< "--------------------------"<< std::endl;

	tmap<<"EigenValues D"<<std::endl;
	tmap<< eigenValM<<std::endl;
	tmap<<"EigenVectors E"<<std::endl;
	tmap<< eigenVecM<<std::endl;
	tmap<<"InvEigenVectors E-1"<<std::endl;
	tmap<< invEigenValM<<std::endl;
	tmap<<"================\n\n"<<std::endl;

	tmap<<"S*Sigma"<<std::endl;
	tmap<< sigmaSFMatrix<<std::endl;

	tmap<<"Compositon: E*D*E-1"<<std::endl;
	tmap<< eigenVecM * eigenValM * invEigenValM<< std::endl;*/



#endif

}


void ThickTracker::trackParticles_m(
#ifdef PHIL_WRITE
    std::ofstream& outfile,
#endif

    const std::list<structMapTracking>& mapBeamLine) {
    int sliceidx=0;
    
    FVector<double, 6> particle, partout;

    dumpStats(sliceidx, true, true);

    double betagamma=itsBunch_m->getInitialBeta() * itsBunch_m->getInitialGamma();

    //(1) Loop Beamline
    for(auto mapBeamLineit=mapBeamLine.begin(); mapBeamLineit != mapBeamLine.end(); ++mapBeamLineit) {
        //(2) Loop Slices
        for (std::size_t slice=0; slice < mapBeamLineit->nSlices; slice++){
            //(3) Loop Particles
            for (unsigned int partidx=0; partidx< itsBunch_m->getLocalNum(); ++partidx){


                for (int d = 0; d < 3; ++d) {
                    particle[2 * d] = itsBunch_m->R[partidx](d);
                    particle[2 *d + 1] = itsBunch_m->P[partidx](d);
                }
                
#ifdef PHIL_WRITE
                if (sliceidx==0) outfile << sliceidx <<"  "<< partidx << " ["<< particle;
#endif
                //Units
                particle[1]/=betagamma;
                particle[3]/=betagamma;

                particle[5] = (particle[5]*itsBunch_m->getM()/itsBunch_m->getP()) //TODO change P to P0
                                * std::sqrt( 1./(particle[5]* particle[5]) +1)
                                -1./itsBunch_m->getInitialBeta();

                //Apply map_t
                particle= mapBeamLineit->elementMap * particle;


                //Units back
                particle[1]*=betagamma;
                particle[3]*=betagamma;

                particle[5] = (particle[5] + 1./itsBunch_m->getInitialBeta()) * itsBunch_m->getP()/itsBunch_m->getM() //TODO change P to P0
                                /std::sqrt( 1./(itsBunch_m ->get_part(partidx)[5]* itsBunch_m ->get_part(partidx)[5]) +1) ; 

#ifdef PHIL_WRITE
                //Write in File
                partout=particle;
                partout[4]+=mapBeamLineit->elementPos + mapBeamLineit->stepSize * (slice+1);
                //outfile<< "partout" << mapBeamLineit->elementPos << "+" << mapBeamLineit->stepSize * slice+1<< std::endl;
                outfile << sliceidx+1 <<"  "<< partidx << " ["<< partout;
#endif


                itsBunch_m->set_part(particle, partidx);



            }
            bool const psDump = true;   //((itsBunch_m->getGlobalTrackStep() % Options::psDumpFreq) + 1 == Options::psDumpFreq);
            bool const statDump = true; //((itsBunch_m->getGlobalTrackStep() % Options::statDumpFreq) + 1 == Options::statDumpFreq);

            changeDT();
            itsBunch_m->set_sPos(mapBeamLineit->elementPos + mapBeamLineit->stepSize * (slice+1));
            dumpStats(sliceidx, psDump, statDump);
            sliceidx ++;
        }
    }    
}

//TODO: Write a nice comment
void ThickTracker::eigenDecomp(fMatrix_t& M, cfMatrix_t& eigenVal, cfMatrix_t& eigenVec, cfMatrix_t& invEigenVec){

	double data[4 * DIM * DIM];
	int idx, s;
	for (int i = 0; i < 2 * DIM; i++) {
		for (int j = 0; j < 2 * DIM; j++) {
			idx = i * 2 * DIM + j;
			data[idx] = M[i][j];
		}
	}


	gsl_matrix_view m = gsl_matrix_view_array(data, 2 * DIM, 2 * DIM);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(2 * DIM);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(2 * DIM, 2 * DIM);
	gsl_matrix_complex *eveci = gsl_matrix_complex_alloc(2 * DIM, 2 * DIM);
	gsl_permutation * p = gsl_permutation_alloc(2 * DIM);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(2 * DIM);

	//get Eigenvalues and Eigenvectors
	gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free(w);
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	for (int i = 0; i < 2 * DIM; ++i) {
		eigenVal[i][i] = complex<double>(
				GSL_REAL(gsl_vector_complex_get(eval, i)),
				GSL_IMAG(gsl_vector_complex_get(eval, i)));
		for (int j = 0; j < 2 * DIM; ++j) {
			eigenVec[i][j] = complex<double>(
					GSL_REAL(gsl_matrix_complex_get(evec, i, j)),
					GSL_IMAG(gsl_matrix_complex_get(evec, i, j)));
		}
	}

	//invert Eigenvectormatrix
	gsl_linalg_complex_LU_decomp(evec, p, &s);
	gsl_linalg_complex_LU_invert(evec, p, eveci);

	//Create invEigenVecMatrix
	for (int i = 0; i < 2 * DIM; ++i) {
		for (int j = 0; j < 2 * DIM; ++j) {
			invEigenVec[i][j] = complex<double>(
					GSL_REAL(gsl_matrix_complex_get(eveci, i, j)),
					GSL_IMAG(gsl_matrix_complex_get(eveci, i, j)));
		}
	}


	//free space
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_matrix_complex_free(eveci);


}

//Analyzes a TransferMap for the tunes, symplecticity and stability1
void ThickTracker::linTAnalyze(fMatrix_t& tMatrix){

	fMatrix_t blocktMatrix;

	cfMatrix_t eigenValM, eigenVecM, invEigenVecM;
	eigenDecomp(tMatrix,eigenValM, eigenVecM, invEigenVecM);


	cfMatrix_t cblocktMatrix = getBlockDiagonal(tMatrix, eigenVecM,invEigenVecM);


	std::ofstream twiss;
	twiss.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/twiss.txt",std::ios::app);
	twiss << std::setprecision(16);


	for (int i=0; i<2*DIM; i++){
			for (int j =0; j<2*DIM; j++){
				blocktMatrix[i][j]=cblocktMatrix[i][j].real();
			}

	}

	cfMatrix_t teigenValM, teigenVecM, tinvEigenVecM;
	eigenDecomp(blocktMatrix, teigenValM, teigenVecM, tinvEigenVecM);

	//twiss << "=====================" << std::endl;
	//twiss << tMatrix << std::endl;
	//twiss << cblocktMatrix;

	//twiss << "the EigenValues of Transfer" << std::endl;
	for (int i=0; i<2*DIM; i++){
	    //twiss << eigenValM[i][i] << "  " << std::endl;
	}

	//twiss << "the abs() of EigenValues of Transfer => Stability (equal 1)" << std::endl;
    for (int i=0; i<2*DIM; i++){
        //twiss << std::abs(eigenValM[i][i]) << "  " << std::endl;
    }

	FVector<std::complex<double>, DIM> betaTunes, betaTunes2, betaTunes3;
	//twiss << "Beta tunes in 2 ways" << std::endl;
	//twiss << tMatrix << std::endl;
	for (int i = 0; i < DIM; i++){
	    betaTunes[i]=std::log(eigenValM[i*2][i*2])/ (2*Physics::pi * complex<double>(0, 1));
	    betaTunes2[i]= std::acos(cblocktMatrix[i*2][i*2])/(complex<double>(2*Physics::pi, 0));
	    betaTunes3[i]= std::acos(eigenValM[i*2][i*2].real())/(2*Physics::pi);
	    twiss<<"1: "<<betaTunes[i] <<std::endl;
	    twiss<<"2: "<<betaTunes2[i]<< std::endl;
	    twiss<<"3: "<<betaTunes3[i]<< std::endl;
	}


	//TODO: do something with the tunes etc
}

//TODO: Write a nice comment
void ThickTracker::linSigAnalyze(){
	fMatrix_t sigmaBlockM;
	fMatrix_t sigMatrix = itsBunch_m->getSigmaMatrix();

	fMatrix_t skewMatrix;
	    for (int i=0; i <6 ; i=i+2){
	    	skewMatrix[0+i][1+i]=1.;
	    	skewMatrix[1+i][0+i]=-1.;
	    }

	sigMatrix=sigMatrix*skewMatrix;

	cfMatrix_t eigenValM, eigenVecM, invEigenVecM;

	eigenDecomp(sigMatrix,eigenValM, eigenVecM, invEigenVecM);

	cfMatrix_t cblocksigM = getBlockDiagonal(sigMatrix, eigenVecM,invEigenVecM);


	std::ofstream tmap;
	//tmap.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/TransferMap.txt",std::ios::app);
	//tmap << std::setprecision(16);



	for (int i=0; i<2*DIM; i++){
			for (int j =0; j<2*DIM; j++){
				sigmaBlockM[i][j]=cblocksigM[i][j].real();
				//if (j==2*DIM-1) tmap<<std::endl;
			}

	}
	//tmap<< "sigmaBlock"<< std::endl;
	//tmap<< sigmaBlockM;
//	tmap<< "\n\n\n" << std::endl;
//	for (int i=0; i<2*DIM; i++){
//				for (int j =0; j<2*DIM; j++){
//					tmap << cSigmaBlockM[i][j].imag() << "  ";
//					if (j==5) tmap<<std::endl;
//				}
//		}

	cfMatrix_t teigenValM, teigenVecM, tinvEigenVecM;
	eigenDecomp(sigmaBlockM, teigenValM, teigenVecM, tinvEigenVecM);

	//tmap << "the EigenValues of Sigmamap => Beam emitances" << std::endl;
	for (int i=0; i<2*DIM; i++){
	//tmap << eigenValM[i][i] << "  " << std::endl;
	}

//	tmap << "the EigenValues of RE(BlockTM)" << std::endl;
//	for (int i=0; i<2*DIM; i++){
//		tmap << teigenValM[i][i] << "  " << std::endl;
//		}

	//TODO: Where to go with EigenValues?
}

//TODO: Write a nice comment
ThickTracker::cfMatrix_t ThickTracker::getBlockDiagonal(fMatrix_t& M,
		cfMatrix_t& eigenVecM, cfMatrix_t& invEigenVecM){

	cfMatrix_t cM, qMatrix, invqMatrix, nMatrix, invnMatrix, rMatrix;

	std::ofstream tmap;
	//tmap.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/TransferMap.txt",std::ios::app);
	//tmap << std::setprecision(16);


	for (int i=0; i<2*DIM; i++){
		for (int j =0; j<2*DIM; j++){
			cM[i][j]=complex<double>(M[i][j],0);
		}
	}


	for (int i=0; i <6 ; i=i+2){
		qMatrix[0+i][0+i]=complex<double>(1.,0);
		qMatrix[0+i][1+i]=complex<double>(0,1.);
		qMatrix[1+i][0+i]=complex<double>(1.,0);
		qMatrix[1+i][1+i]=complex<double>(0,-1);

		invqMatrix[0+i][0+i]=complex<double>(1.,0);
		invqMatrix[0+i][1+i]=complex<double>(1.,0);
		invqMatrix[1+i][0+i]=complex<double>(0.,-1.);
		invqMatrix[1+i][1+i]=complex<double>(0,1.);
		}
	qMatrix/=std::sqrt(2.);
	invqMatrix/=std::sqrt(2);

	nMatrix=eigenVecM*qMatrix;
	invnMatrix= invqMatrix* invEigenVecM;

	rMatrix= invnMatrix * cM * nMatrix;

//	tmap<< "Qmatrix"<< std::endl;
//	tmap<< qMatrix<< std::endl;
//	tmap<< "invQ"<< std::endl;
//	tmap<< invqMatrix<< std::endl;
//
//	tmap<< "Q*Q-1"<< std::endl;
//	tmap<< qMatrix* invqMatrix<< std::endl;
//
//	tmap<< "NMatrix"<< std::endl;
//	tmap<< nMatrix<< std::endl;
//	tmap<< "invNMatrix"<< std::endl;
//	tmap<< invnMatrix<< std::endl;
//	tmap<< "N* invNMatrix"<< std::endl;
//	tmap<< nMatrix*invnMatrix << std::endl;


	return rMatrix;

}

void ThickTracker::dumpStats(long long step, bool psDump, bool statDump) {

    OPALTimer::Timer myt2;
    Inform msg("ThickTracker", *gmsg);

    std::size_t numParticlesInSimulation_m = itsBunch_m->getTotalNum();

    if (itsBunch_m->getGlobalTrackStep() % 10 + 1 == 10) {
        msg << level1;
    } else if (true){ //itsBunch_m->getGlobalTrackStep() % 1 + 1 == 1) {
        msg << level2;
    } else {
        msg << level3;
    }

    if (numParticlesInSimulation_m == 0) {
        msg << myt2.time() << " "
            << "Step " << std::setw(6) <<  itsBunch_m->getGlobalTrackStep() << "; "
            << "   -- no emission yet --     "
            << "t= "   << Util::getTimeString(itsBunch_m->getT())
            << endl;
        return;
    }

    itsBunch_m->calcEMean();
    //size_t totalParticles_f = numParticlesInSimulation_m;

    if (std::isnan(pathLength_m) || std::isinf(pathLength_m)) {
        throw OpalException("ParallelTTracker::dumpStats()",
                            "there seems to be something wrong with the position of the bunch!");
    } else {

        msg << myt2.time() << " "
            << "Step " << std::setw(6)  <<  itsBunch_m->getGlobalTrackStep() << " "
            << "at "   << Util::getLengthString(pathLength_m) << ", "
            << "t= "   << Util::getTimeString(itsBunch_m->getT()) << ", "
            << "E="    << Util::getEnergyString(itsBunch_m->get_meanKineticEnergy())
            << endl;

        writePhaseSpace(step, psDump, statDump);
    }
}

void ThickTracker::writePhaseSpace(const long long step, bool psDump, bool statDump) {
    extern Inform *gmsg;
    Inform msg("OPAL ", *gmsg);
    Vector_t externalE, externalB;
    Vector_t FDext[2];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);

    if (psDump || statDump) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(referenceToLabCSTrafo_m.transformTo(RefPartR_m),
                                     referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
                                     itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(),
                                     externalE,
                                     externalB);
        FDext[0] = referenceToLabCSTrafo_m.rotateFrom(externalB);
        FDext[1] = referenceToLabCSTrafo_m.rotateFrom(externalE * 1e-6);
    }

    if (statDump) {
        std::vector<std::pair<std::string, unsigned int> > collimatorLosses;
        FieldList collimators = itsOpalBeamline_m.getElementByType(ElementBase::COLLIMATOR);
        if (collimators.size() != 0) {
            for (FieldList::iterator it = collimators.begin(); it != collimators.end(); ++ it) {
                Collimator* coll = static_cast<Collimator*>(it->getElement().get());
                std::string name = coll->getName();
                unsigned int losses = coll->getLosses();
                collimatorLosses.push_back(std::make_pair(name, losses));
            }
            std::sort(collimatorLosses.begin(), collimatorLosses.end(),
                      [](const std::pair<std::string, unsigned int>& a, const std::pair<std::string, unsigned int>& b) ->bool {
                          return a.first < b.first;
                      });
            std::vector<unsigned int> bareLosses(collimatorLosses.size(),0);
            for (size_t i = 0; i < collimatorLosses.size(); ++ i){
                bareLosses[i] = collimatorLosses[i].second;
            }

            reduce(&bareLosses[0], &bareLosses[0] + bareLosses.size(), &bareLosses[0], OpAddAssign());

            for (size_t i = 0; i < collimatorLosses.size(); ++ i){
                collimatorLosses[i].second = bareLosses[i];
            }
        }
        
        // Write statistical data.
        itsDataSink_m->writeStatData(itsBunch_m, FDext, collimatorLosses);

        msg << level3 << "* Wrote beam statistics." << endl;
    }

    if (psDump && (itsBunch_m->getTotalNum() > 0)) {
        // Write fields to .h5 file.
        const size_t localNum = itsBunch_m->getLocalNum();
        double distToLastStop = zStop_m.back() - pathLength_m;
        Vector_t driftPerTimeStep = itsBunch_m->getdT() * Physics::c * RefPartP_m / Util::getGamma(RefPartP_m);
        bool driftToCorrectPosition = std::abs(distToLastStop) < 0.5 * euclidean_norm(driftPerTimeStep);
        Ppos_t stashedR;

        if (driftToCorrectPosition) {
            const double tau = distToLastStop / euclidean_norm(driftPerTimeStep) * itsBunch_m->getdT();
            if (localNum > 0) {
                stashedR.create(localNum);
                stashedR = itsBunch_m->R;

                for (size_t i = 0; i < localNum; ++ i) {
                    itsBunch_m->R[i] += tau * (Physics::c * itsBunch_m->P[i] / Util::getGamma(itsBunch_m->P[i]) -
                                               driftPerTimeStep / itsBunch_m->getdT());
                }
            }

            itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m + tau * driftPerTimeStep / itsBunch_m->getdT());
            CoordinateSystemTrafo update(tau * driftPerTimeStep / itsBunch_m->getdT(),
                                         Quaternion());
            itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m * update.inverted();

            itsBunch_m->set_sPos(zStop_m.back());

            itsBunch_m->calcBeamParameters();
        }
        if (!statDump && !driftToCorrectPosition) itsBunch_m->calcBeamParameters();

        msg << *itsBunch_m << endl;
        itsDataSink_m->writePhaseSpace(itsBunch_m, FDext);

        if (driftToCorrectPosition) {
            if (localNum > 0) {
                itsBunch_m->R = stashedR;
                stashedR.destroy(localNum, 0);
            }

            itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
            itsBunch_m->set_sPos(pathLength_m);

            itsBunch_m->calcBeamParameters();
        }

        msg << level2 << "* Wrote beam phase space." << endl;
    }
}

void ThickTracker::setTime() {
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}


