#ifndef BEAMSTRIPPINGPHYSICS_HH
#define BEAMSTRIPPINGPHYSICS_HH
//Class:BeamStrippingPhysics
//  Defines the beam stripping physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
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
class BeamStripping;

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

    bool GasStripping(double &deltas);

    bool LorentzStripping(double &gamma, double &E);

    void SecondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i);

    Cyclotron *cycl_m;
    BeamStripping *bstp_m;

    gsl_rng * r_m;

    std::string material_m;
    std::string FN_m;
    ElementBase::ElementType bstpshape_m;

    double  T_m;
    double dT_m;

    double mass_m;
    double charge_m;

	double m_h;

	int NbComponents;
	//    double totalmolecularDensity_m;
	double molecularDensity[3];

    std::unique_ptr<LossDataSink> lossDs_m;

    double NCS_single_all;
    double NCS_double_all;
    double NCS_total_all;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    static const double fMolarFraction[3];
	static const double CSCoefSingle_Hminus[3][7];
	static const double CSCoefDouble_Hminus[3][7];
	static const double CSCoefSingle_Hplus[3][9];
	static const double CSCoefDouble_Hplus[3][9];
	static const double CSCoefSingleLoss_H[3][7];
	static const double CSCoefSingleCapt_H[3][9];
};

#endif //BEAMSTRIPPINGPHYSICS_HH
