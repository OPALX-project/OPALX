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
#include "AbsBeamline/Cyclotron.h"
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
class Cyclotron;

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

    void setCyclotron(Cyclotron* cycl) { cycl_m = cycl; };

private:

    void Material();

    void MolecularDensity(const double &pressure, const double &temperature, int &iComp);
    void CrossSection(double &Eng);

    double CSAnalyticFunction(double Eng, double Eth,
    		double a1, double a2, double a3, double a4, double a5, double a6);

    double RandomGenerator();

    bool GasStripping(double &deltas, double &r);

    bool LorentzStripping(double &gamma, double &E, double &r);

    double  T_m;                     // own time, maybe larger than in the bunch object
    double dT_m;                     // dt from bunch
    
    double mass;
    double charge;
    double molecularDensity[3];
    double CS[3];

    std::string material_m;
    std::string FN_m;
    ElementBase::ElementType bstpshape_m;
    std::string strippingStr_m;

	int NbComponents;
	static const double fMolarFraction[3];

	static const double CSCoefSingle_Hminus[3][7];
	static const double CSCoefDouble_Hminus[3][7];
	static const double CSCoefSingle_Hplus[3][9];
	static const double CSCoefDouble_Hplus[3][9];

	/*
	static const double fCrossSectionSingle[3][48];
	static const double fEnergyCSSingle[3][48];
	static const double fCrossSectionDouble[3][40];
	static const double fEnergyCSDouble[3][40];
	*/

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    Vector_t extE_m, extB_m;
    Cyclotron *cycl_m;
};

#endif //BEAMSTRIPPINGPHYSICS_HH
