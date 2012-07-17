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

#include "Solvers/SurfacePhysicsHandler.hh"
//#include "Algorithms/PBunchDefs.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"

class RANLIB_class;
class ElementBase;
class PartBunch;
class LossDataSink;

class CollimatorPhysics: public SurfacePhysicsHandler {
public:
    CollimatorPhysics(const std::string &name, ElementBase *element, const double &major, const double &minor, std::string &mat);

    void apply(PartBunch &bunch);

    virtual const std::string getType() const;

    ~CollimatorPhysics();
    void Material();
    void EnergyLoss(double &Eng, bool &pdead);
    void CoulombScat(Vector_t &R, Vector_t &P, double &deltat, double scalefactor);
    void CoulombScat();
    void EnergyLoss(double &Eng, bool &pdead, double &deltat);
    void  Rot(Vector_t &P, Vector_t Prot, double normP);
    double Rot(double &p1, double &p2, double &scatang);

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

    bool incoll_m;
    Point  geom_m[5];

    int index_m;
    std::vector<unsigned> label_m;
    std::vector<Vector_t> Rincol_m;
    std::vector<Vector_t> Pincol_m;
    std::vector<long> IDincol_m;
    std::vector<int> Binincol_m;
    std::vector<double> DTincol_m;
    std::vector<double> Qincol_m;
    std::vector<long> LastSecincol_m;
    std::vector<Vector_t> Bfincol_m;
    std::vector<Vector_t> Efincol_m;
    std::vector<double> time_m;
    std::vector<int> steps_m;
  
    void setCColimatorGeom();
    int  checkPoint( const double & x, const double & y );
    LossDataSink *lossDs_m;
};

#endif //COLLIMATORPHYSICS_HH
