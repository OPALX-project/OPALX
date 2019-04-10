#ifndef COLLIMATORPHYSICS_HH
#define COLLIMATORPHYSICS_HH
//Class:CollimatorPhysics
//  Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang, Stachel, Adelmann$
//-------------------------------------------------------------------------
#include <vector>
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Degrader.h"
#include <gsl/gsl_rng.h>

#include "Utility/IpplTimings.h"

#ifdef OPAL_DKS
#include "DKSOPAL.h"
#endif

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;
class Inform;

#ifdef OPAL_DKS
typedef struct __align__(16) {
    int label;
    unsigned localID;
    Vector_t Rincol;
    Vector_t Pincol;
    long IDincol;
    int Binincol;
    double DTincol;
    double Qincol;
    long LastSecincol;
    Vector_t Bfincol;
    Vector_t Efincol;
} PART;

typedef struct {
    int label;
    unsigned localID;
    Vector_t Rincol;
    Vector_t Pincol;
} PART_DKS;

#else
typedef struct {
    int label;
    unsigned localID;
    Vector_t Rincol;
    Vector_t Pincol;
    long IDincol;
    int Binincol;
    double DTincol;
    double Qincol;
    long LastSecincol;
    Vector_t Bfincol;
    Vector_t Efincol;
} PART;
#endif


class CollimatorPhysics: public ParticleMatterInteractionHandler {
public:
    CollimatorPhysics(const std::string &name, ElementBase *element, std::string &mat);
    ~CollimatorPhysics();

    void apply(PartBunchBase<double, 3> *bunch,
               const std::pair<Vector_t, double> &boundingSphere,
               size_t numParticlesInSimulation = 0);

    virtual const std::string getType() const;

    void print(Inform& os);
    bool stillActive();
    bool stillAlive(PartBunchBase<double, 3> *bunch);

    inline double getTime() {return T_m;}
    std::string getName() { return FN_m;}
    size_t getParticlesInMat() { return locPartsInMat_m;}
    unsigned getRediffused() { return rediffusedStat_m;}

    inline void doPhysics(PartBunchBase<double, 3> *bunch);


private:

    void Material();
    void CoulombScat(Vector_t &R, Vector_t &P, const double &deltat);
    void EnergyLoss(double &Eng, bool &pdead, const double &deltat);
    bool EnergyLoss(double &Eng, const double &deltat);

    void Rot(double &px, double &pz, double &x, double &z, double xplane, double Norm_P,
             double thetacou, double deltas, int coord);

    void copyFromBunch(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere);
    void addBackToBunch(PartBunchBase<double, 3> *bunch, unsigned i);

#ifdef OPAL_DKS
  void copyFromBunchDKS(PartBunchBase<double, 3> *bunch,
                        const std::pair<Vector_t, double> &boundingSphere);
    void addBackToBunchDKS(PartBunchBase<double, 3> *bunch, unsigned i);

    void setupCollimatorDKS(PartBunchBase<double, 3> *bunch, size_t numParticlesInSimulation);
    void clearCollimatorDKS();

    void applyDKS();
    void applyHost(PartBunchBase<double, 3> *bunch, Degrader *deg, Collimator *coll);
    void deleteParticleFromLocalVectorDKS();

#endif


    void deleteParticleFromLocalVector();

    void calcStat(double Eng);

    double  T_m;                     // own time, maybe larger than in the bunch object

    double dT_m;                     // dt from bunch

    gsl_rng *rGen_m;

    std::string material_m;
    std::string FN_m;
    ElementBase::ElementType collshape_m;
    std::string collshapeStr_m;

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

    // some statistics

    double Eavg_m;
    double Emax_m;
    double Emin_m;



    std::vector<PART> locParts_m;

    std::unique_ptr<LossDataSink> lossDs_m;

#ifdef OPAL_DKS
    DKSOPAL dksbase;
    int curandInitSet;

    int ierr;
    int maxparticles;
    int numparticles;
    int numlocalparts;
    void *par_ptr;
    void *mem_ptr;

    std::vector<PART_DKS> dksParts_m;

    static const int numpar;
#endif

    IpplTimings::TimerRef DegraderApplyTimer_m;
    IpplTimings::TimerRef DegraderLoopTimer_m;
    // :FIXME: remove unused declaration
    // IpplTimings::TimerRef DegraderInitTimer_m;
    IpplTimings::TimerRef DegraderDestroyTimer_m;
};

inline
void CollimatorPhysics::calcStat(double Eng) {
    Eavg_m += Eng;
    if (Emin_m > Eng)
        Emin_m = Eng;
    if (Emax_m < Eng)
        Emax_m = Eng;
}

inline
void CollimatorPhysics::EnergyLoss(double &Eng, bool &pdead, const double &deltat) {
    pdead = EnergyLoss(Eng, deltat);
}

#endif //COLLIMATORPHYSICS_HH

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
