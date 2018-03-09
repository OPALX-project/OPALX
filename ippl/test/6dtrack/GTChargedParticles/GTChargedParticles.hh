#ifndef CHARGEDPARTICLES_HH
#define CHARGEDPARTICLES_HH
#include "Ippl.h"
#include "GTConst.hh"
#include "GTConfigure.hh"
#include "GTElemData.hh"

#include "GTPhysics.h"

#include "Algorithms/PartData.h"

typedef ParticleSpatialLayout<double, 6> playout_t;

template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {


  typedef UniformCartesian<6,double>               Mesh_t;
  typedef Cell                                     Center_t;
  typedef CenteredFieldLayout<6, Mesh_t, Center_t> FieldLayout_t;
  typedef Field<double, 6, Mesh_t, Center_t>       Field_t;


  typedef IntCIC IntrplCIC_t;
  typedef IntNGP IntrplNGP_t;

 public:
  typedef Vektor<double,6> Vector_t;
  typedef Vektor<double,3> Vector3;



 private:

  Field_t  rho_m;
  BConds<double,6,Mesh_t,Center_t> bc_m;

  Vector_t rmin_m;
  Vector_t rmax_m;
  Vector_t hr_m;
  Vector_t nr_m;

  double centroid_m[2*3];
  double moments_m[2*3][2*3];

  Vector3 rmean_m, pmean_m;
  Vector3 rrms_m, prms_m, eps_m, rprms_m;
  Vector3 csbeta_m,csgamma_m,csalpha_m;

  // macroparticle charge in Cb
  double q_m;
  double spos_m;
  double actTime_m;

  ofstream statof_m;
  ofstream partof_m;
  unsigned int partofCall_m;

public:

  // the attributes for this set of particles (atoms).

  ParticleAttrib< double > Q;

  PartData &pdata_m;
  SimCfgData &simCfg_m;

    ChargedParticles(PartData &pd, SimCfgData &simCfg) :
    pdata_m(pd),
    simCfg_m(simCfg),
    q_m(simCfg.Qtot/simCfg.nEInit),  // macroparticle charge in Cb
    partofCall_m(0)
  {
    // initialize the base class, by creating a new layout
    Inform msg("ChargedParticles ");
    Mesh_t *mesh;
    FieldLayout_t *FL;
    NDIndex<6> domain;
    e_dim_tag decomp[6];

    Vektor<int,6> nr(simCfg.nx,simCfg.ny,simCfg.nz,
		     simCfg.nx,simCfg.ny,simCfg.nz);

    for(int i=0; i<6; i++)
	domain[i] = domain[i] = Index(nr[i] + 1);

    for (int d=0; d < 6; ++d) {
	if (d<4)
	    decomp[d] = PARALLEL;
	else
	    decomp[d] = SERIAL;
    }
    // create mesh and layout objects for this problem domain

    mesh = new Mesh_t(domain);
    FL   = new FieldLayout_t(*mesh, decomp);

    initialize(new PL(*FL,*mesh));

    // register our attributes with the base class
    addAttribute(Q);

    for (int i=0; i < 2*6; ++i) {
	bc_m[i] = new ZeroFace<double,6,Mesh_t,Center_t>(i);
	getBConds()[i] = ParticleNoBCond;
    }

    rho_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<6>(1), bc_m);

    msg << " Field and domain construction done = " << domain << endl;

  }

  /**
   * Interpolates the charge in Q onto a grid rho_m
   * @param none
   * @return none
   */
  void scatterQ() {
    rho_m = 0.0;
    Q.scatter(rho_m, this->R, IntrplNGP_t());
  }

  /**
   * Shows how to work with fields and validated the
   * charge deposition method. Is used if FIELDCHEK
   * is defined
   * @param none
   * @return none
   */
    void checkScatteredQ() {

      Inform m("checkScatteredQ ");

      double s1Q = sum(Q);
      double s2Q = sum(rho_m);
      m << "s1Q= " << s1Q << " s2Q= " << s2Q << " delta= " << abs(s1Q-s2Q) << endl;
  }


  double getLocalVolumes(Vector_t r) {
    double ret=0;
    NDIndex<6> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<6> elem = getMesh().getCellContaining(r);
    if (lIdx.contains(elem)) {
      ret = getMesh().getCellVolume(elem);
    }
    return ret;
  }


  double getLocalValues(Vector_t r) {
    double ret=0;
    NDIndex<6> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<6> elem = getMesh().getCellContaining(r);
    if (lIdx.contains(elem)) {
      ret = rho_m.localElement(elem);
    }
    return ret;
  }




  /**
   * Calculats the local density at a given particle position r
   * @param r position of the particle
   * @return local density in [Q m**-3]
   */
  double getLocalDensity(Vector_t r) {
    double den=0;
    NDIndex<6> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<6> elem = getMesh().getCellContaining(r);
    if (lIdx.contains(elem)) {
      double lDen = rho_m.localElement(elem);
      double lVol = getMesh().getCellVolume(elem);
      den = lDen/lVol;
    }
    return den;
  }

  /**
   * Calculats the analytic local density at a given particle position r
   * @param r position of the particle
   * @param p momentum of the particle
   * @return local density in [Q m**-3]
   */
  double getLocalDensityAnalytic(Vector_t r, Vector_t p) {
    Vector_t e = get_emit();
    Vector_t a = get_csalpha();
    Vector_t b = get_csbeta();
    Vector_t g = get_csgamma();

    Vector_t rho;

    for(unsigned i = 0; i < 3; i++)
      rho(i) = 1 / (2 * Physics::pi * e(i)) * exp(- abs(1/(2 * e(i)) * (g(i) * r(i) * r(i) + 2 * a(i) * r(i) * p(i) + b(i) * p(i) * p(i))));

    return rho(0) * rho(1) * rho(2);
  }
 /**
   * Prototypes implemented in GTChargedParticles.cpp
   *
   *
   *
   */

  void openFiles(string baseFn, string title,bool restartMode);
  void closeFiles();
  void calcBeamParameters();
  void calcMoments();

  void writeStatSDDSHeader(string title, unsigned int N, double qTot);
  void writeStatistics();

  void writePhaseSpaceSDDS(string baseName);
  void writePartSDDSHeader(string title, unsigned int N, double qTot);

  /**
   * Resets if the check if a particle has collided
   */
  void resetCollisions() { collided = 0; }
  /**
   * Checks if a particle has collided.
   * @param i ID of the particle to check for collision
   * @return true if it did, false otherwise
   */
  inline bool hasCollided(unsigned i) {return collided[i]!=0;}
  /**
   * collides 2 particles
   * @param i first particle to collide
   * @param j second particle to collide
   */
  inline void collide(unsigned i, unsigned j) {
    collided[i] = 1;
    collided[j] = 1;
  }

  inline const Mesh_t& getMesh() const { return getLayout().getLayout().getMesh(); }
  inline       Mesh_t& getMesh() { return getLayout().getLayout().getMesh(); }
  inline const FieldLayout_t& getFieldLayout() const {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }
  inline       FieldLayout_t& getFieldLayout() {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }

  inline bool isRoot() { return Ippl::myNode()==0; }

  void boundp() {

    bounds(R,rmin_m,rmax_m);

    NDIndex<6> domain = getFieldLayout().getDomain();
    for(int i=0; i<6; i++)
      nr_m[i] = domain[i].length();

    Vektor<double,6> len = (rmax_m - rmin_m);

    for(int i=0; i<6; i++)
      hr_m[i]    = (len[i] / (nr_m[i]-1));

    // rescale mesh
    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin( rmin_m );

    //    BinaryRepartition(*this);
    update();
    rho_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<6>(1), bc_m);
  }

 /**
   *
   *
   *
   * @return the macroparticle charge in Cb
   */
  inline double getMoment(int i, int j) {return moments_m[i][j];}


 /**
   *
   *
   *
   * @return the macroparticle charge in Cb
   */
  //void writePartSDDSHeader(string title, unsigned int N, double qTot);


 /**
   *
   *
   *
   * @return the macroparticle charge in Cb
   */
  inline double getQ()    { return q_m;}


  /**
   *
   *
   *
   * @return the actual simulation time [s]
   */
  inline double getTime()  { return spos_m; }


  /**
   *
   *
   *
   * @return increment the simulation time [s]
   */
  inline void incTime(double tInc) { actTime_m += tInc;}


  /**
   *
   *
   *
   * @return the relativistic factor []
   */
  inline double getGamma() {return pdata_m.getGamma();}


  /**
   *
   *
   *
   * @return beta v/c []
   */
  inline double getBeta()  {return pdata_m.getBeta();}


  /**
   *
   *
   *
   * @return get the actual position along the design orbit [m]
   */
  inline double set_spos(double s) { spos_m = s; }

  /**
   *
   *
   *
   * @return get the actual position along the design orbit [m]
   */
  inline double get_spos() { return spos_m; }



  /**
   *
   *
   *
   * @return get the total charge of the bunch [Cb]
   */
  inline double getQTot() { return getTotalNum()*q_m; }



  /**
   *    * --------- max
   *    *     0/0
   *    * --------- *
   *
   * @return get the maximal (max) extent of the bunch [m]
   */
  inline Vector_t get_maxExtend() { return rmax_m; }


  /**
   *    * --------- *
   *    *     0/0
   *  min --------- *
   *
   * @return get the minimal (min) extent of the bunch [m]
   */
  inline Vector_t get_origin() { return rmin_m; }


  inline Vector_t get_hr() { return hr_m; }



  inline Vector_t get_emit() const { return eps_m; }
  inline Vector_t get_rprms() const { return rprms_m; }
  inline Vector_t get_csbeta() const { return csbeta_m; }
  inline Vector_t get_csgamma() const { return csgamma_m; }
  inline Vector_t get_csalpha() const { return csalpha_m; }

  inline Vector_t get_rmean() const { return rmean_m; }
  inline Vector_t get_pmean() const { return pmean_m; }
  inline Vector_t get_rrms() const { return rrms_m; }
  inline Vector_t get_prms() const { return prms_m; }

  inline Vektor<double,6> getR(unsigned i) { return this->R[i]; }
  inline void setR(unsigned i, Vektor<double,6> x) { this->R[i] = x; }

  /**
   * @return the turn number to restart
   */
  int readRestartInfo(string Fn);
  void writeRestartInfo(string basename, unsigned turn);

  void getInbalance(string fn, bool first) {

    unsigned long int locN = getLocalNum();
    unsigned long int idealN = getTotalNum()/Ippl::getNodes();
    ofstream of;

    /*
       Get the particle inbalance from all the nodes
    */
    double  *locBal  = (double*) malloc(Ippl::getNodes()*sizeof(double));
    double  *globBal = (double*) malloc(Ippl::getNodes()*sizeof(double));

    for(int i=0; i<Ippl::getNodes(); i++)
      locBal[i]=globBal[i]=0.0;

    locBal[Ippl::myNode()] = locN;

    reduce(locBal, locBal + Ippl::getNodes(), globBal, OpAddAssign());

    if (first && isRoot()) {
      of.open(fn.c_str(),ios::out);
      of << "# inbalance N_i, i=1.." << Ippl::getNodes() << " and Ntotal " << endl;
      for(int i=0; i<Ippl::getNodes(); i++)
	of << globBal[i] << "\t";
      of << getTotalNum() << endl;
      of.close();
    }
    else if (isRoot()) {
      of.open(fn.c_str(),ios::app);
      for(int i=0; i<Ippl::getNodes(); i++)
	of << globBal[i] << "\t";
      of << getTotalNum() << endl;
      of.close();
    }
    if (locBal)
      free(locBal);
    if (globBal)
      free(globBal);
  }
};
#include "GTChargedParticles.cpp"
#endif