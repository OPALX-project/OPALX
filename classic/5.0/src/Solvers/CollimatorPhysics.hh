#ifndef COLLIMATORPHYSICS_HH
#define COLLIMATORPHYSICS_HH
//Class:CollimatorPhysics
//  Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: bi $
//-------------------------------------------------------------------------
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Algorithms/PartBunch.h"
#include "Distribution/ranlib.h"
#include <vector>

class ElementBase;

class CollimatorPhysics: public SurfacePhysicsHandler {
public:
    CollimatorPhysics(const string &name, ElementBase *element, const double &major, const double &minor, string &mat);

    void apply(PartBunch &bunch);

    virtual const string getType() const;

    ~CollimatorPhysics();
    void Material();
    void EnergyLoss(double &Eng, bool &pdead);
    void CoulombScat(Vector_t &R, Vector_t &P, double &deltat, double scalefactor);
    void CoulombScat();
    void EnergyLoss(double &Eng, bool &pdead, double &deltat);
    void  Rot(Vector_t &P, Vector_t Prot, double normP);
    double Rot(double &p1, double &p2, double &scatang);
    double calculateAngle(double x, double y);
private:
    RANLIB_class *rGen_m;
    double a_m;
    double b_m;
    double xp_m;
    double yp_m;
    double angstart_m;
    double angend_m;
    double rstart_m;
    double rend_m;
    double width_m;

    double Begin_m;
    double End_m;
    string material_m;
    string FN_m;
    string collshape_m;
    int Z_m;
    double A_m;
    double rho_m;
    double X0_m;
    double I_m;
    double n_m;

    bool incoll_m;

    int index_m;
    vector<unsigned> label_m;
    vector<Vector_t> Rincol_m;
    vector<Vector_t> Pincol_m;
    vector<long> IDincol_m;
    vector<int> Binincol_m;
    vector<double> DTincol_m;
    vector<double> Qincol_m;
    vector<long> LastSecincol_m;
    vector<Vector_t> Bfincol_m;
    vector<Vector_t> Efincol_m;
    vector<double> time_m;
    vector<int> steps_m;

    LossDataSink *lossDs_m;
};

#endif //COLLIMATORPHYSICS_HH
