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

class ElementBase;
class PartBunch;
class LossDataSink;
class Inform;

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


class CollimatorPhysics: public SurfacePhysicsHandler {
public:
    CollimatorPhysics(const std::string &name, ElementBase *element, std::string &mat);
    ~CollimatorPhysics();

    void apply(PartBunch &bunch);

    virtual const std::string getType() const;

    void print(Inform& os);
    bool stillActive() { return bunchToMatStat_m != 0;}

   inline double getTime() {return T_m;}

private:

    void Material();
    void CoulombScat(Vector_t &R, Vector_t &P, double &deltat);
    void EnergyLoss(double &Eng, bool &pdead, double &deltat);

    void Rot(double &px, double &pz, double &x, double &z, double xplane, double Norm_P, double thetacou, double deltas, int coord);

    void copyFromBunch(PartBunch &bunch);
    void addBackToBunch(PartBunch &bunch, unsigned i);
    void deleteParticleFromLocalVector();

    bool checkHit(Vector_t R, Vector_t P, double dt, Degrader *deg, Collimator *coll); 

    inline void calcStat(double Eng) {
      Eavg_m += Eng;
      if (Emin_m > Eng)
	Emin_m = Eng;
      if (Emax_m < Eng)
	Emax_m = Eng;
    }

    bool allParticlesIn_m;
  
    double  T_m;                     // own time, maybe larger than in the bunch object
                                    
    double dT_m;                     // dt from bunch

    gsl_rng *rGen_m;

    std::string material_m;
    std::string FN_m;
    std::string collshape_m;
    int Z_m;
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
    unsigned redifusedStat_m;

    // some statistics

    double Eavg_m;
    double Emax_m;
    double Emin_m;

    std::vector<PART> locParts_m;

    std::unique_ptr<LossDataSink> lossDs_m;

};

#endif //COLLIMATORPHYSICS_HH
