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
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Degrader.h"
#include <gsl/gsl_rng.h>

#include "Utility/IpplTimings.h"

#ifdef OPAL_DKS
#include "DKSOPAL.h"
#endif

class ElementBase;
class PartBunch;
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


class CollimatorPhysics: public SurfacePhysicsHandler {
public:
    CollimatorPhysics(const std::string &name, ElementBase *element, std::string &mat,
		      bool enableRutherfordScattering, double lowEnergyThr);
    ~CollimatorPhysics();

    void apply(PartBunch &bunch, size_t numParticlesInSimulation = 0);

    virtual const std::string getType() const;

    void print(Inform& os);
    bool stillActive();
    bool stillAlive(PartBunch &bunch);

    inline double getTime() {return T_m;}
    std::string getName() { return FN_m;}
    size_t getParticlesInMat() { return globalPartsInMat_m;}
    unsigned getRedifused() { return redifusedStat_m;}

    void doPhysics(PartBunch &bunch, Degrader *deg, Collimator *col);


private:

    void Material();
    void applyCoulombScat(Vector_t &R, Vector_t &P, double &deltat);
    bool computeEnergyLoss(double &Eng, double &deltat);

    void applyDKS(PartBunch &bunch, size_t numParticlesInSimulation);

    void computeScatteringEffect(Vector_t &P, Vector_t &R,
                                 double xplane, double Norm_P,
                                 double thetacou, double deltas, unsigned char coord);

    void copyFromBunch(PartBunch &bunch);
    void addBackToBunch(PartBunch &bunch, unsigned i);

#ifdef OPAL_DKS
    void copyFromBunchDKS(PartBunch &bunch);
    void addBackToBunchDKS(PartBunch &bunch, unsigned i);

    void setupCollimatorDKS(PartBunch &bunch, Degrader *deg, size_t numParticlesInSimulation);
    void clearCollimatorDKS();

    void applyDKS();
    void applyHost(PartBunch &bunch, Degrader *deg, Collimator *coll);
    void deleteParticleFromLocalVectorDKS();

#endif


    void deleteParticleFromLocalVector();

    bool checkHit(Vector_t R, Vector_t P, double dt, Degrader *deg, Collimator *coll);

    inline void calcStat(double Eng) {
        Eavg_m += Eng;
        if (Emin_m > Eng)
            Emin_m = Eng;
        if (Emax_m < Eng)
            Emax_m = Eng;
    }

    double  T_m;                     // own time, maybe larger than in the bunch object

    double dT_m;                     // dt from bunch

    gsl_rng *rGen_m;

    enum SHAPE {
    CYCLCOLLIMATOR,
        COLLIMATOR,
        DEGRADER,
        NOSHAPE
    };

    std::string material_m;
    std::string FN_m;
    std::string collShapeStr_m;
    SHAPE collShape_m;
    double Z_m; // target atomic number
    double A_m; // target atomic weight
    double A2_c; // ??
    double A3_c; // ??
    double A4_c; // ??
    double A5_c; // ??
    double rho_m; // target mass density
    double X0_m; // radiation length
    double I_m; // ?target mean excitation?
    // double n_m;

    bool enableRutherfordScattering_m;
    double lowEnergyThr_m;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned redifusedStat_m;
    size_t localPartsInMat_m;
    size_t globalPartsInMat_m;

    // some statistics

    double Eavg_m;
    double Emax_m;
    double Emin_m;

    std::vector<PART> locParts_m;
    size_t nextLocalID_m;

    std::unique_ptr<LossDataSink> lossDs_m;

#ifdef OPAL_DKS
    DKSOPAL dksbase_m;
    int curandInitSet_m;

    int ierr_m;
    int maxparticles_m;
    int numparticles_m;
    void *par_mp;
    void *mem_mp;

    std::vector<PART_DKS> dksParts_m;

    static const int numpar_ms = 13;
#endif

    IpplTimings::TimerRef DegraderApplyTimer_m;
    IpplTimings::TimerRef DegraderLoopTimer_m;
    IpplTimings::TimerRef DegraderInitTimer_m;


};

inline
const std::string CollimatorPhysics::getType() const {
    return "CollimatorPhysics";
}

#endif //COLLIMATORPHYSICS_HH