#ifndef BEAMSTRIPPINGPHYSICS_HH
#define BEAMSTRIPPINGPHYSICS_HH

// ------------------------------------------------------------------------
//
// Class: BeamStrippingPhysics
//   Defines the beam stripping physics models
//
// ------------------------------------------------------------------------
// Class category: ParticleMatterInteraction
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

    BeamStrippingPhysics(const std::string &name, ElementBase *element);
    ~BeamStrippingPhysics();
    
    void setCyclotron(Cyclotron* cycl) { cycl_m = cycl; };
    
    void apply(PartBunchBase<double, 3> *bunch,
               const std::pair<Vector_t, double> &boundingSphere,
               size_t numParticlesInSimulation = 0);

    virtual const std::string getType() const;
    void print(Inform& msg);
    bool stillActive();
    bool stillAlive(PartBunchBase<double, 3> *bunch);

    inline double getTime() {return T_m;}
    std::string getName() {return element_ref_m->getName();}
    size_t getParticlesInMat() {return locPartsInMat_m;}
    unsigned getRediffused() {return rediffusedStat_m;}
    unsigned int getNumEntered() {return bunchToMatStat_m;}
    inline void doPhysics(PartBunchBase<double, 3> *bunch);
    
private:

    void CrossSection(const Vector_t &R, double &Eng);

    double CSAnalyticFunctionNakai(double Eng, double Eth, int &i);
    
    double CSAnalyticFunctionTabata(double Eng, double Eth,
                            double a1, double a2, double a3, double a4, double a5, double a6);
                              
    double CSChebyshevFitting(double Eng, double Emin, double Emax);

    bool GasStripping(double &deltas);

    bool LorentzStripping(double &gamma, double &E);

    void SecondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i, bool pdead_LS);
    void TransformToProton(PartBunchBase<double, 3> *bunch, size_t &i);
    void TransformToHydrogen(PartBunchBase<double, 3> *bunch, size_t &i);
    void TransformToHminus(PartBunchBase<double, 3> *bunch, size_t &i);
    void TransformToH3plus(PartBunchBase<double, 3> *bunch, size_t &i);

    bool computeEnergyLoss(double &Eng,
                           const double deltat,
                           bool includeFluctuations = true) const { return false;}

    Cyclotron *cycl_m;
    BeamStripping *bstp_m;
    ElementBase::ElementType bstpshape_m;

    gsl_rng * r_m;
    
    double T_m;
    double dT_m;

    double mass_m;
    double charge_m;

    double m_h;
    
    double totalmolecularDensity_m;
    double molecularDensity[3];
    std::string gas_m;
    double pressure_m;
    ///@{ size limits took from cyclotron
    double rhoinit_m; 
    double rhofinal_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    double NCS_a;
    double NCS_b;
    double NCS_c;
    double NCS_total;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    static const double CSCoefSingle_Hminus[3][9];
    static const double CSCoefDouble_Hminus[3][9];
    static const double CSCoefSingle_Hplus[3][9];
    static const double CSCoefDouble_Hplus[3][9];
    static const double CSCoefSingleLoss_H[3][9];
    static const double CSCoefSingleCapt_H[3][9];
    
    static const double CSCoefHminusProduction_H_Tabata[13];
    static const double CSCoefProtonProduction_H_Tabata[9];
    static const double CSCoefProtonProduction_H2plus_Tabata[11];
    static const double CSCoefH3plusProduction_H2plus_Tabata[7];
    
    static const double CSCoefSingle_Hminus_Chebyshev[11];
    static const double CSCoefDouble_Hminus_Chebyshev[11];
    static const double CSCoefSingle_Hplus_Chebyshev[11];
    static const double CSCoefDouble_Hplus_Chebyshev[11];
    static const double CSCoefHydrogenProduction_H2plus_Chebyshev[11];
    
    static double a_m[9];
    static double b_m[3][9];
};

#endif //BEAMSTRIPPINGPHYSICS_HH
