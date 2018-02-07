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
#include <cmath>

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


    msg << *itsBunch_m << endl;
    msg << std::setprecision(20);

//=======================================================================


    //generate the Hamiltonian
    typedef FTps<double, 2 * DIM> Series;
    typedef FVps<double, 2 * DIM> Map;

    Series x = Series::makeVariable(0);			//SIXVect::X);
    Series px = Series::makeVariable(1);		//SIXVect::PX);
    Series y = Series::makeVariable(2);			//SIXVect::Y);
    Series py = Series::makeVariable(3);		//SIXVect::PY);
    //Series z = Series::makeVariable(4);		//SIXVect::TT);
    Series delta = Series::makeVariable(5);		//SIXVect::PT);



    int order = 8;  //actually already defined
    bool ModeSlices= true;

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


    Series::setGlobalTruncOrder(order+1);

    ///Creates the Hamiltonian for beam line element
    /**\param element iterative pointer to the elements along the beam line
    * \return Hamiltonian times elemenrLength
    */
    auto Hamiltonian = [&](std::shared_ptr<Component> element) {

        Series H; //Hamiltonian

        switch(element->getType()) {

            /**Driftspace
             * \f[H_{Drift}= \frac{\delta}{\beta_0} -
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } \f]
             */
            case ElementBase::ElementType::DRIFT: {
                //Drift* pDrift= dynamic_cast<Drift*> (element.get());
                H=( delta / beta0 )
                - sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
                        - ( px*px )
                        - ( py*py )
                        - 1./( beta0 * beta0 * gamma0 * gamma0 ),order+1
                );
                
                Drift* drift_p = dynamic_cast<Drift*> (element.get());
                
                msg << "drift NSLICES: " << drift_p->getNSlices() << endl;

                break;
            }

            /**Rectangular Bend
             * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
             * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right) \f]
             */
            case ElementBase::ElementType::RBEND: {
                RBend* pRBend= dynamic_cast<RBend*> (element.get());

                double h = 1. / pRBend ->getBendRadius();               //inverse bending radius [1/m]
                double K0= pRBend ->getB()*(Physics::c/itsReference.getP());

                H=( delta / beta0 )
                - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order+1
                        ))
                - (h * x)
                * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order
                        ))
                + K0 * x * (1. + 0.5 * h* x);

                H= H /pRBend->getElementLength() * pRBend->getArcLength();

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
                H=( delta / beta0 )
                - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order+1
                        ))
                - (h * x)
                * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                - ( px*px )
                                - ( py*py )
                                - 1./( beta0*beta0 * gamma0*gamma0 ),order
                        ))
                + K0 * x * (1. + 0.5 * h* x);


                H= H /pSBend->getElementLength() *pSBend->getEffectiveLength();

                break;
            }

            /**Quadrupole "in Multipolegroup"
             * \f[H_{Quadrupole}= \frac{\delta}{\beta_0} -
             * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2  -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }	+
             * \frac{1}{2} k_1 \left( x^2 - y^2 \right) \f]
             */
            case ElementBase::ElementType::MULTIPOLE: {
                Multipole* pMultipole= dynamic_cast<Multipole*> (element.get());

                msg << "NSLICES = " << pMultipole->getNSlices() << endl;
                double K1= pMultipole->getField().getNormalComponent(2)*(Physics::c/P0);
                K1= std::round(K1*1e6)/1e6;

                H= ( delta / beta0 )
                - sqrt ((1./ beta0 + delta ) *(1./ beta0 + delta)
                        - ( px*px )
                        - ( py*py )
                        - 1./( beta0*beta0 * gamma0*gamma0 ),order+1
                )
                + 0.5 * (q * K1) *(Physics::c/P0) * (x*x - y*y);

                break;
            }

            default:
            throw LogicalError("ThickTracker::exercute,",
                    "Please use already defined beam line element:  At this time just driftspace, multipoles and dipoles");
            break;

        }
        return H;
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

    Map combinedMaps, ElementMap;

    FVector<double, 6> particle;
    OpalParticle part;



    msg << "P0 OPAL:  " << itsBunch_m-> getP() << endl;
    msg << "E0 OPAL:  " << itsBunch_m-> getM() << endl;


    std::map<std::string, double> elementDic;
    std::list<std::pair<int, Map>> mapBeamLine;

    std::map<std::string, double>::iterator dicit;

    std::vector<Map> mapVec;

    double selement =0;
    double elementpos =0;
    double stepsize= itsBunch_m->getdT()* Physics::c * itsBunch_m->getInitialBeta();
    msg << "getdT():  " << stepsize << endl;
    unsigned int nSlices, totalSlices=0;
    double residualSliceL;
    //bool initializedBunch=false;


    //files for analysis
    std::ofstream outfile;
    outfile.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/generatedMaps.txt");
    outfile << std::setprecision(8);

    std::ofstream tmap;
    tmap.open ("/home/phil/Documents/ETH/MScProj/OPAL/src/tests/Maps/TransferMap.txt");
    tmap << std::setprecision(20);

    //loop over beam line
    for (; it != end; ++ it) {

        dicit = elementDic.begin();

        std::shared_ptr<Component> element = (*it).getElement();

        elementpos= std::round(element->getElementPosition() * 1e6)/1e6;

        dicit= elementDic.find(element->getName());

        //order beam line according the element position
        if (dicit !=elementDic.end()){
            throw LogicalError("ThickTracker::execute,",
                                "Same Element twice in beamline:"+ element->getName());

        }else{
            elementDic.insert(std::pair<std::string, double> (element->getName(), elementpos));
        }

        if (selement > elementpos){
            msg << "There is an overlap! @ Element: " << element ->getName() <<
                    "-> starts at: " << elementpos<<
                    " the previous ends at: " << selement << endl;
//            throw LogicalError("ThickTracker::exercute,",
//                                            "Overlap at element:"+ element->getName());

        }else{
            selement=std::round(elementpos+ element ->getElementLength()*1e6)/1e6;
        }


        //combine maps in between Monitors

        if (element ->getType()==ElementBase::MONITOR){
            mapVec.insert(mapVec.end(), combinedMaps);
            combinedMaps.identity();

        }


//=================================
// TODO: remove Messages later
//=================================

        msg << "=============================="<< endl;
        msg << "Element name: " << element->getName() << endl;
        msg << "Element Length: " << element->getElementLength() << endl;

        msg <<  "ElementPosition"<< elementpos<< endl;
        msg << "EntryPoint" << element->getEntrance()<< endl;

//===================================




        //Use Slices for a better visualization
        if (ModeSlices){
        	nSlices = element->getElementLength()/ stepsize;
        	residualSliceL =  std::fmod(element->getElementLength(), stepsize);

        	mapBeamLine.push_back(std::pair<int, Map> (nSlices, ExpMap(-Hamiltonian(element) * stepsize  ,order)));
        	mapBeamLine.push_back(std::pair<int, Map> (1, ExpMap(-Hamiltonian(element) * residualSliceL ,order)));

//=================================
// TODO: remove Messages later
//=================================

        	msg << "stepsize:   " << stepsize << endl;
        	msg << "nSlices:   " << nSlices << endl;
        	msg << "residualSliceL:   " << residualSliceL << endl;

//================================
//

    	}else{

    		ElementMap=ExpMap(-Hamiltonian(element) * element->getElementLength()  ,order);

    		mapBeamLine.push_back(std::pair<int, Map> (1, ElementMap));

        }
    }


	for(std::pair<int,Map>  &elementidx : mapBeamLine) {
		for (int slice=0; slice < elementidx.first; slice++){
			combinedMaps= elementidx.second * combinedMaps;
			totalSlices++;
		}
	}

	outfile.seekp(0);
	outfile <<"Total Particle Number:  " <<  itsBunch_m->getTotalNum() << "  Total number of Slices:   "<< totalSlices << std::endl;



    msg << "totalSlices:   " << totalSlices << endl;
    msg << "totalDistance  " << totalSlices * stepsize << endl;

    //eliminate higher order terms (occurring through multiplication)
	combinedMaps=combinedMaps.truncate(order);

	tmap << "A customized FODO" << std::endl;
	tmap << combinedMaps << std::endl;


	//track the Particles

	int sliceidx;

	//(1) Loop Particles
	for (unsigned int partidx=0; partidx< itsBunch_m->getTotalNum(); ++partidx){
		msg << "inParticleLoop:   " << partidx << endl;
		sliceidx=0;

		part =itsBunch_m ->get_part(partidx);
		for(int dim=0; dim<6; dim++){
			particle[dim]=itsBunch_m ->get_part(partidx)[dim];
			if (dim<5) particle[dim]= std::round(particle[dim] * 1e6)/1e6;
		}
		outfile << sliceidx <<"  "<< partidx << " ["<< particle;

		//(2) Loop Elements
		for(std::pair<int,Map>  &elementidx : mapBeamLine) {

			//(3) Loop Slices
			for (int slice=0; slice < elementidx.first; slice++){
				combinedMaps= elementidx.second * combinedMaps;

				//Units

				particle[5] = (particle[5]*itsBunch_m->getM()/itsBunch_m->getP())
					* std::sqrt( 1./(particle[5]* particle[5]) +1)
					-1./itsBunch_m->getInitialBeta();

				//Apply Map
				particle= elementidx.second * particle;
				sliceidx ++;

				//Units back
				particle[4]+= stepsize;

				particle[5] = (particle[5] + 1./itsBunch_m->getInitialBeta()) * itsBunch_m->getP()/itsBunch_m->getM()
						/std::sqrt( 1./(itsBunch_m ->get_part(partidx)[5]* itsBunch_m ->get_part(partidx)[5]) +1) ;

				//Write in File
				outfile << sliceidx <<"  "<< partidx << " ["<< particle;
			}

		}

	}





	std::map <std::string, double>::iterator mi;
	  for (mi=elementDic.begin(); mi!=elementDic.end(); ++mi)
	      msg << mi->first << " = " << mi->second << endl;

}
