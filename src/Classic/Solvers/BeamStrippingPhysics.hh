#ifndef BEAMSTRIPPINGPHYSICS_HH
#define BEAMSTRIPPINGPHYSICS_HH
//Class:BeamStrippingPhysics
//  Defines the beam stripping physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang, Stachel, Adelmann$
//-------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamStripping.h"
#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include <vector>

#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

#include "Utility/IpplTimings.h"

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;

class LogicalError;
class LossDataSink;
class Inform;

class BeamStrippingPhysics: public ParticleMatterInteractionHandler {
public:
    BeamStrippingPhysics(const std::string &name, ElementBase *element, std::string &mat);
    ~BeamStrippingPhysics();

    void apply(PartBunchBase<double, 3> *bunch,
               const std::pair<Vector_t, double> &boundingSphere,
               size_t numParticlesInSimulation = 0);

    virtual const std::string getType() const;

    void print(Inform& msg);
    bool stillActive();
    bool stillAlive(PartBunchBase<double, 3> *bunch);

    inline double getTime() {return T_m;}
    std::string getName() {return FN_m;}
    size_t getParticlesInMat() {return locPartsInMat_m;}
    unsigned getRediffused() {return rediffusedStat_m;}

    inline void doPhysics(PartBunchBase<double, 3> *bunch);


private:

    void Material();
    void GasDensity(double &pressure, double &temperature, int &iComp);
    void CrossSection(double &Eng);
    void FractionLost(double &Eng);
    bool GasStripping(double &r);
    double RandomGenerator();

    double  T_m;                     // own time, maybe larger than in the bunch object
    double dT_m;                     // dt from bunch
    
    double mass;
    double gasDensity[3];
    double CS[3];
    double fg;

//    gsl_rng *rGen_m;

    std::string material_m;
    std::string FN_m;
    ElementBase::ElementType bstpshape_m;
    std::string strippingStr_m;

	int NbComponents;
	static const double fMolarFraction[3];
	static const double fCrossSectionSingle[3][48];
	static const double fEnergyCSSingle[3][48];
	static const double fCrossSectionDouble[3][40];
	static const double fEnergyCSDouble[3][40];

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    std::unique_ptr<LossDataSink> lossDs_m;
};

#endif //BEAMSTRIPPINGPHYSICS_HH
