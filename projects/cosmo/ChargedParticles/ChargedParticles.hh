
#ifndef CHARGED_PARTICLES_HH
#define CHARGED_PARTICLES_HH

#include "Ippl.h"
#include "FieldSolver.hh"

#include "Parser/FileStream.h"
#include "SimData.hh"
#include "PwrSpec.hh"

#include <iostream>
#include <fstream>
#include <strstream>


template <class T, unsigned int Dim>
class ChargedParticles : public ParticleBase< ParticleSpatialLayout<T,Dim> > {

private:
  
public:
  
  typedef IntCIC IntrplCIC_t; 

  typedef typename ParticleSpatialLayout<T,Dim>::ParticlePos_t Ppos_t;
  typedef typename ParticleSpatialLayout<T,Dim>::ParticleIndex_t PID_t;
  
  typedef UniformCartesian<Dim,T> Mesh_t;
  
  typedef typename ParticleSpatialLayout<T,Dim>::SingleParticlePos_t Vector_t;
  
  typedef ParticleSpatialLayout<T, Dim, Mesh_t>      Layout_t;
  
  //typedef Cell                                       Center_t;
  typedef Vert                                       Center_t;
  
  typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
  typedef Field<T, Dim, Mesh_t, Center_t>            Field_t;
  typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
  
  ParticleAttrib<T>          M; // mass
  ParticleAttrib<Vector_t>   V;      // particle velocity
  ParticleAttrib<Vector_t>   F;      // force field at particle position

  ParticleAttrib<Vector_t>   rtmp;   

  T ti_m;                               // actual integration time
  T dt_m;                               // timestep

  Vector_t rmin_m, rmax_m;
  
  // Grid stuff
  Field_t  rho_m;                       // mass density on grid and scalar field
  VField_t gphi_m;                      // gradient of the scalat field on the grid
  
  Vector_t hr_m;
  Vektor<int,Dim> nr_m;
  
  BConds<T,Dim,Mesh_t,Center_t> bc_m;
  BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;
  
  FieldSolver<T,Dim> *fs_m;
  
  bool fieldNotInitialized_m;


  //constructor
  ChargedParticles(Layout_t *playout) : 
    ParticleBase< ParticleSpatialLayout<T,Dim> >(playout)
  {
    // register the particle attributes
    addAttribute(M);
    addAttribute(V);
    addAttribute(F); 
    addAttribute(rtmp); 
    
    setBCAllPeriodic();
    
    // random number check
    IpplRandom.SetSeed(234244117);
    double dummy = IpplRandom();
    
    fieldNotInitialized_m = true;

  }

  void genIC(SimData simData) {
    Inform msg ("genIC ");

    msg << "Start generating initial conditions" <<  endl;
 
    unsigned int count = 0;
    unsigned int NP = simData.np*simData.np*simData.np;
    
    if (isRoot()) {
      for (unsigned int i=0; i<NP;i++) { 
	create(1);
	double arg1 = 1.0;
	this->R[count](0) = arg1;
	this->R[count](1) = arg1;
	this->R[count](2) = arg1;
	this->V[count](0) = arg1;
	this->V[count](1) = arg1;
	this->V[count](2) = arg1;
	count++;
      }
    }
   
    this->M = 1.0;
    this->rho_m;
    this->M.scatter(this->rho_m, this->R, IntCIC());

    msg << "Particles created: " << getLocalNum() << endl;
    //    rescaleAll();
    // msg2all << "Particles after update: " << getLocalNum() << endl;
  } 



  inline const Mesh_t& getMesh() const { return getLayout().getLayout().getMesh(); }
 
  inline Mesh_t& getMesh() { return getLayout().getLayout().getMesh(); }
 
  inline const FieldLayout_t& getFieldLayout() const {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }
 
  inline FieldLayout_t& getFieldLayout() {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }

  inline void set_FieldSolver(FieldSolver<T,Dim> *fs) { fs_m = fs; }

  inline T getTime() {return ti_m; }
  inline void setTime(T ti) {ti_m = ti; }

  inline T getdT() {return dt_m; }
  inline void setdT(T dt) {dt_m = dt; }
  
  inline Vector_t getGridSize() { return nr_m; }
  inline Vector_t getMeshSpacing() { return hr_m; }
  inline Vector_t getRmin() { return rmin_m; }
  inline Vector_t getRmax() { return rmax_m; }
  
  inline bool isRoot() { return (Ippl::Comm->myNode() == 0); }
    
  void rescaleAll() {  
    
    bounds(this->R, rmin_m, rmax_m);

    Vector_t len;
    len = rmax_m - rmin_m;
    
    NDIndex<Dim> domain = getFieldLayout().getDomain(); 
    for(int i=0; i<Dim; i++) {
      nr_m[i] = domain[i].length();
      hr_m[i] = (len[i] / (nr_m[i]-1));
    }
    
    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(rmin_m);
    
    //    gphi_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
    
    gphi_m.initialize(getMesh(), getFieldLayout(), vbc_m);
    rho_m.initialize(getMesh(), getFieldLayout(),  bc_m);
    update();
  }
  
  
  
  inline void do_binaryRepart() {
    BinaryRepartition(*this);
    boundp();
  }

  inline void setBCAllPeriodic() {
    for (int i=0; i < 2*Dim; ++i) {
      bc_m[i] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
      getBConds()[i] = ParticlePeriodicBCond;
    }
  }

  inline void setBCAllOpen() {
    for (int i=0; i < 2*Dim; ++i) {
      bc_m[i] = new ZeroFace<T,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
      getBConds()[i] = ParticlePeriodicBCond;
    }
  }
  
   
  void readInputDistribution(string fn) {
    Inform msg2all("readInputDistribution ", INFORM_ALL_NODES);
   
    if (isRoot()) {
      double  arg1, arg2, arg3, arg4, arg5, arg6;
      ifstream *ifstr = new ifstream;
      ifstr->open(fn.c_str(), ios::in );
      unsigned int count = 0;
      
      *ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;
      
      while(*ifstr) {
	create(1);
	this->R[count](0) = arg1;
	this->R[count](1) = arg3;
	this->R[count](2) = arg5;
	this->V[count](0) = arg2;
	this->V[count](1) = arg4;
	this->V[count](2) = arg6;
	count++;
	*ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;
      }
      ifstr->close();
    }
    msg2all << "Particles created: " << getLocalNum() << endl;
    rescaleAll();
    msg2all << "Particles after update: " << getLocalNum() << endl;
    //    do_binaryRepart();
    // msg2all << "Particles after BinaryRepartition: " << getLocalNum() << endl;
  }
    

  void dumpParticles(string fname) {
  
    Inform msg("dumpParticles ");
  
    unsigned int pwi = 14;
    ofstream of;
    unsigned int ok = 0;
    int toSave=0;

    int tag1 = IPPL_APP_TAG1;
    int tag2 = IPPL_APP_TAG2;
  
    msg << "dumping " << getTotalNum() << endl;
  
    if (singleInitNode()) {
      of.open(fname.c_str(),ios::out);
      for (unsigned int ii = 0 ; ii < getLocalNum(); ++ii) {
	for (int jj = 0 ; jj < Dim ; ++jj) {
	  of << setw(pwi) << this->R[ii](jj) << ' ';
	  of << setw(pwi) << this->V[ii](jj) << ' ';
	}
	of << endl;
      }
      of.close(); 
      
      for (toSave=1; toSave <= Ippl::getNodes()-1; toSave++) {
	// notify toSave that he can save
	Message *msg = new Message();
	msg->put(toSave);
	bool res = Ippl::Comm->send(msg, toSave, tag1); 
	if (! res) 
	  ERRORMSG("Ippl::Comm->send(smsg, 0, tag1) failed " << endl;);
	// wail til toSave reports done with save
	Message* doneMsg;
	doneMsg = Ippl::Comm->receive_block(toSave,tag2);
	doneMsg->get(&ok);
	delete doneMsg;
	if (ok != toSave) {
	  ERRORMSG("Node0 got wrong msg from node " << toSave);
	  ERRORMSG("received ID = " << ok  << endl);
	}
      }
    }
    else {
      // Wail til we get ready to save
      Message* readyToSaveMesg;
      int z=0;
      
      readyToSaveMesg = Ippl::Comm->receive_block(z, tag1);
      readyToSaveMesg->get(&ok);
      delete readyToSaveMesg;
      if (ok == Ippl::myNode()) {
	of.open(fname.c_str(),ios::app);
	for (unsigned int ii = 0 ; ii < getLocalNum(); ++ii) {
	  for (int jj = 0 ; jj < Dim ; ++jj) {
	    of << setw(pwi) << this->R[ii](jj) << ' ';
	    of << setw(pwi) << this->V[ii](jj) << ' ';
	  }
	  of << endl;
	}
	of.close();
	// Notify node0 we are done by sending the node id
	Message *msg = new Message();
	msg->put(ok);
	bool res = Ippl::Comm->send(msg, 0, tag2);
	if (! res) 
	  ERRORMSG("Ippl::Comm->send(smsg, 0, tag2) failed " << endl;);
      }
      else {
	ERRORMSG("Node" << Ippl::myNode() << "received ID = " << ok  << endl);
      }
    }
    Ippl::Comm->barrier(); 
  }
  
  
  // CALCULATES AND OUTPUTS PHYSICAL DATA OF BEAM
  inline void writeDiagnostics() {

    
  } 
  
  inline void map2(double tau, double omegatot,double alpha, double pp) {
    /// adot in units of H_0^-1
    Inform msg("map2: ");
    if (fs_m) {
      fs_m->computeForceField(this);
    
      double adot = sqrt((omegatot+(1.0-omegatot)*pow(pp,(3.0/alpha))) / pow(pp,1.0/alpha) );
      double phiscale = 1.5*omegatot/pow(pp,1.0/alpha);
      double fscale   = phiscale/(alpha*adot*pow(pp,(1.0 - 1.0/alpha)));
      msg << "adot= " << adot << " phiscale= " << phiscale << " fscale= " << fscale << endl; 
      assign(V, V + fscale*tau*F);
    }
  }  
  
  inline void map1(double tau, double omegatot,double alpha, double pp) {
    Inform msg("map1: ");
    /// adot in units of H_0^-1
    double adot = sqrt((omegatot+(1.0-omegatot)*pow(pp,(3.0/alpha))) / pow(pp,1.0/alpha) );
    double pre  = 1.0/( alpha*adot*pow(pp,(1.0 + 1.0/alpha)));

    assign(this->R, this->R + tau * pre * this->V);
    update();
  }
  
  void testConvolution() { 
    fs_m->testConvolution(this);
  }

  void print(ostream &os) 
  {
    os << " ------ Conditions @ t= " << getTime() << " ------ " << endl;
    os << " N= " << getTotalNum() << " Grid= " << getGridSize() << endl;
    os << " h= " << getMeshSpacing() << " []" << endl ;
    os << " rmax=" << getRmax() << " []" << endl;
    os << " rmin=" << getRmin() << " []" << endl;
  }
};

inline ostream &operator<<(ostream &os, ChargedParticles<double,3> &data)
{
  data.print(os);
  return os;
}
#endif 

