#ifndef CHARGEDPARTICLES_HH
#define CHARGEDPARTICLES_HH
#include "Ippl.h"
#include "GTConst.hh"
#include "GTConfigure.hh"
#include "GTElemData.hh"

#include "GTPhysics.h"

#include "Algorithms/PartData.h"

typedef ParticleInteractLayout<double, 3> playout_t;

template<class PL>
class ChargedParticles : public ParticleBase<PL> {

  
  typedef UniformCartesian<3,double>               Mesh_t;
  typedef Cell                                     Center_t;
  typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
  typedef Field<double, 3, Mesh_t, Center_t>       Field_t;


  typedef IntCIC IntrplCIC_t;
  typedef IntNGP IntrplNGP_t;

 public:
  typedef Vektor<double,3> Vector_t;



 private:
    
  Field_t  rho_m;
  BConds<double,3,Mesh_t,Center_t> bc_m;

  Vektor<double,3> rmin_m;
  Vektor<double,3> rmax_m;
  Vektor<double,3> pmin_m;
  Vektor<double,3> pmax_m;
  Vektor<double,3> hr_m;
  Vektor<int,3> nr_m;

  double centroid_m[2*3];
  double moments_m[2*3][2*3];

  Vector_t rmean_m, pmean_m;
  Vector_t rrms_m, prms_m, eps_m, rprms_m;
  Vector_t csbeta_m,csgamma_m,csalpha_m;
  
  // macroparticle charge in Cb
  double q_m;
  double spos_m;
  double actTime_m;

  ofstream statof_m;
  ofstream partof_m;
  unsigned int partofCall_m;

  // electron statistics
  unsigned neAtWall_m;
  double   eSecElec_m;
  double   maxEsecElec_m;
  unsigned nsTot_m;

  // collision statistics
  unsigned collision_num;
  int lost_num;
  int lost2_num;

  double interrad_m;

  // crosssection statistics
  double deltav;
  double deltav2;
  double deltau;
  double deltau_analytic;
  double deltau_semianalytic;
  unsigned events;
  
public:

  // the attributes for this set of particles (atoms).  

  ParticleInteractAttrib< Vektor<double,3> > P;
  ParticleInteractAttrib< double > Q;
  ParticleInteractAttrib< unsigned > collided;
  ParticleInteractAttrib< unsigned > lost;

  PartData &pdata_m;
  SimCfgData &simCfg_m;

  ChargedParticles(double interrad, PartData &pd, SimCfgData &simCfg, int parallelDim, int meshMultFactor):
    pdata_m(pd),
    simCfg_m(simCfg),
    q_m(simCfg.Qtot/simCfg.nEInit),  // macroparticle charge in Cb
    partofCall_m(0),
    interrad_m(interrad)
  {
    // initialize the base class, by creating a new layout
    Inform msg("ChargedParticles ");

    IpplRandom.SetSeed(static_cast<unsigned long>(23131719));
	
    Mesh_t *mesh;
    FieldLayout_t *FL;
    NDIndex<3> domain;
    e_dim_tag decomp[3];

    Vektor<int,3> nr(simCfg.nx,simCfg.ny,simCfg.nz);
    
    for(int i=0; i<3; i++)
      domain[i] = domain[i] = Index(nr[i] + 1);
  
    for (int d=0; d < 3; ++d) 
      decomp[d] = (d == parallelDim) ? PARALLEL : SERIAL; 
  
    // create mesh and layout objects for this problem domain
    
    mesh = new Mesh_t(domain);
    FL   = new FieldLayout_t(*mesh, decomp);

    initialize(new PL(*FL,*mesh));
    
    // register our attributes with the base class
    addAttribute(P);
    addAttribute(Q);
    addAttribute(collided);
    
    // set the interaction radius for the atoms.
    getLayout().setInteractionRadius(interrad_m);
    
    lost_num = 0;
    lost2_num = 0;
    
    for (int i=0; i < 2*3; ++i) {
      bc_m[i] = new ZeroFace<double,3,Mesh_t,Center_t>(i);
      getBConds()[i] = ParticleNoBCond;
    }
    
    rho_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
    
    msg << " Field and domain construction done = " << domain << endl;

    deltau = 0;
    deltav = 0;
    deltau_analytic = 0;
    deltav2 = 0;
    deltau_semianalytic = 0;
    events = 0;
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
    
    NDIndex<Dim> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem;
    double lQ = 0.0;
    double lV = 0.0;
    double gQ = 0.0;
    double gV = 0.0;

    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
      elem[0]=Index(i,i); 
      for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
	elem[1]=Index(j,j);
	for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
	  elem[2]=Index(k,k);	
	  lQ += rho_m.localElement(elem);	
	  lV += getMesh().getCellVolume(elem);
	}
      }
    }
    reduce(lQ,gQ,OpAddAssign());
    reduce(lV,gV,OpAddAssign());
    double s1Q = sum(Q);
    double s2Q = sum(rho_m);

    m << "gQ= " << gQ << " gV= " << gV << " hr= " << hr_m << endl;
    m << "s1Q= " << s1Q << " s2Q= " << s2Q << " delta= " << abs(s1Q-s2Q) << endl;
  }
    

  double getLocalVolumes(Vector_t r) {
    double ret=0;
    NDIndex<Dim> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem = getMesh().getCellContaining(r);
    if (lIdx.contains(elem)) {
      ret = getMesh().getCellVolume(elem);
    }    
    return ret;
  }


  double getLocalValues(Vector_t r) {
    double ret=0;
    NDIndex<Dim> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem = getMesh().getCellContaining(r);
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
    NDIndex<Dim> lIdx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem = getMesh().getCellContaining(r);
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
      rho(i) = 1 / (2 * Physics::pi * e(i)) * exp(- abs(1/(2 * e(i)) * (g(i) * r(i) * r(i) + (i == 2 ? 2 * a(i) * r(i) * p(i) : 0) + b(i) * p(i) * p(i))));
    
    return rho(0) * rho(1) * rho(2);
  }

  /**
   * Calculates the bunch volume
   * @return the bunch volume
   */
  double getBunchVolume() {
    return Physics::pi * 4 / 3 * pow(5.0,3.0) * get_rrms()(0) * get_rrms()(1) * get_rrms()(2);
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
   * Sets the interaction radius. This is the sphere in which getPairList
   * will be looking for particles
   * @param r radius of tghe sphere
   * 
   */  
  void setInteractionRadius(double r) {getLayout().setInteractionRadius(r);}

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
    bounds(P,pmin_m,pmax_m);
    
    NDIndex<3> domain = getFieldLayout().getDomain();
    for(int i=0; i<3; i++)
      nr_m[i] = domain[i].length();
    
    Vektor<double,3> len = (rmax_m - rmin_m);
      
    for(int i=0; i<3; i++)
      hr_m[i]    = (len[i] / (nr_m[i]-1));
      
    // rescale mesh
    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin( rmin_m );

    //    BinaryRepartition(*this);
    update();

    rho_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);                     
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

  inline Vector_t get_R_beam(unsigned i) { return Vector_t(this->R[i](0),this->R[i](1),this->R[i](2) * getGamma());}
  inline Vector_t get_P_beam(unsigned i) { return Vector_t(this->P[i](0),this->P[i](1),this->P[i](2) * getGamma());}

  inline Vektor<double,3> getR(unsigned i) { return this->R[i]; }
  inline void setR(unsigned i, Vektor<double,3> x) { this->R[i] = x; }

  inline Vektor<double,3> getP(unsigned i) { return this->P[i]; }
  inline void setP(unsigned i, Vektor<double,3> p) { this->P[i] = p; }

  // Collision functions

  /**
   * @param the number of collisions during the computation
   */  
  inline void set_collisionnum(unsigned n) {collision_num = n;}
  /** 
   * add_lost increases the number of particles lost by touschek scattering in the
   * longitudinal directions
   * @param i number of lost particles
   */
  inline void add_lost(int i) {lost_num+=i;}
  /** 
   * add_lost2 increases the number of particles lost by touschek scattering in the
   * transversal directions
   * @param i number of lost particles
   */  
  inline void add_lost2(int i) {lost2_num =+ i;}
  /**
   * @return the number of particles lost due to longitudinal scattering
   */
  inline int get_lost_num() {return lost_num;}
  /**
   * @return the number of particles lost due to transversal scattering
   */
  inline int get_lost2_num() {return lost2_num;}
  /**
   * @return the number of collisions
   */
  inline unsigned get_collision_num() {return collision_num;}

  /**
   * adding up for crosssection stuff
   */
  inline void add_delta_v(double d) {deltav += d;}
  inline void add_delta_v2(double d) {deltav2 += d;}
  inline void add_delta_u_semianalytic(double d) {deltau_semianalytic += d;}
  inline void add_delta_u_analytic(double d) {deltau_analytic += d;}
  inline void add_delta_u(double d) {deltau += d;}
  inline void add_events(unsigned i) {events += i;}

  inline double get_delta_v() {return deltav;}
  inline double get_delta_v2() {return deltav2;}
  inline double get_delta_u() {return deltau;}
  inline double get_delta_u_analytic() {return deltau_analytic;}
  inline double get_delta_u_semianalytic() {return deltau_semianalytic;}
  inline unsigned get_events() {return events;}

  /**
   * @return the debye length
   */
  inline double get_debye_length() {return sqrt(Physics::epsilon_0 * Physics::EMASS * dot(get_pmean(),get_pmean()) * pow(Physics::c * getBeta(),2) / ((getQTot() / Physics::q_e)/sqrt(dot(get_rmean(), get_rmean())) * pow(Physics::q_e,2)));}

  /**
   * @return deviation from the courant snyder invariant
   */
  inline double get_csfault(unsigned i) {return get_csgamma()(i) * get_rrms()(i) * get_rrms()(i) + 2 * get_csalpha()(i) * get_rrms()(i) * get_prms()(i) + get_csbeta()(i) * get_prms()(i) * get_prms()(i) - get_emit()(i);}

  /**
   * lose particles
   * @param id local id of particle to lose
   */
  inline  void lose_particle(unsigned id) {lost[id] = 1;}
  void send_lose_particle();
  void write_lose_particle(string file);
  void init_lose_particle(string file);

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
