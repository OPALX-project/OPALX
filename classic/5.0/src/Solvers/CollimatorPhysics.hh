#ifndef COLLIMATORPHYSICS_HH
#define COLLIMATORPHYSICS_HH
//Class:CollimatorPhysics
//  Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang $
//-------------------------------------------------------------------------
#include <vector>
#include "Ippl.h"
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"

class RANLIB_class;
class ElementBase;
class PartBunch;
class LossDataSink;

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
    CollimatorPhysics(const std::string &name, ElementBase *element, const double &major, const double &minor, std::string &mat);

    void apply(PartBunch &bunch);

    virtual const std::string getType() const;

    ~CollimatorPhysics();
    void Material();
    void EnergyLoss(double &Eng, bool &pdead);
    void CoulombScat(Vector_t &R, Vector_t &P, double &deltat);
    void CoulombScat();
    void EnergyLoss(double &Eng, bool &pdead, double &deltat);
    // void  Rot(Vector_t &P, Vector_t Prot, double normP);
    // double Rot(double &p1, double &p2, double &scatang);
    void Rot(double &px, double &pz, double &x, double &z, double xplane, double Norm_P, double thetacou, double deltas, int coord);


private:

    RANLIB_class *rGen_m;
    double a_m;
    double b_m;
    double xp_m;
    double yp_m;
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double zstart_m;
    double zend_m;
    double width_m;

    double Begin_m;
    double End_m;
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

    unsigned matToBunchStat_m;
    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned redifusedStat_m;

    // some statistics

    double Eavg_m;
    double Emax_m;
    double Emin_m;


public:
    double dT_m;
    int N_m;

private:
    bool incoll_m;
    Point  geom_m[5];

    double time_m;

    std::vector<PART> locParts_m;
  
    void setCColimatorGeom();
    int  checkPoint( const double & x, const double & y );
    LossDataSink *lossDs_m;

    bool checkInColl(Vector_t R);

    void copyFromBunch(PartBunch &bunch);

    void addBackToBunch(PartBunch &bunch, unsigned i);

    void deleteParticleFromLocalVector();

public:
    void print(Inform& os);

    bool stillActive() { return bunchToMatStat_m != 0;}
};

#endif //COLLIMATORPHYSICS_HH
