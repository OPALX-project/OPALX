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
    void calcNumMolecules(double &pressure, double &temperature);
    void CrossSection(double &Eng, vector<double> &sigma, vector<double> &energycs);
    void FractionLost(double &Eng);
    bool GasStripping(double &r);
    double RandomGenerator();

    double  T_m;                     // own time, maybe larger than in the bunch object
    double dT_m;                     // dt from bunch
    
    double mass;
    double NumMolecules_m;
    double CS;
    double fg;

//    gsl_rng *rGen_m;

    std::string material_m;
    std::string FN_m;
    ElementBase::ElementType bstpshape_m;
    std::string strippingStr_m;

    double Z_m;
    double A_m;
    double A2_c;
    double A3_c;
    double A4_c;
    double A5_c;
    double rho_m;
    double X0_m;
    double I_m;
    double n_m;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    std::unique_ptr<LossDataSink> lossDs_m;
};

#endif //BEAMSTRIPPINGPHYSICS_HH
