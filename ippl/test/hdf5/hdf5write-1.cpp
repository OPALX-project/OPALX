// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.
 *
 ***************************************************************************/

/***************************************************************************

 This test program sets up a simple sine-wave electric field in 3D, 
   creates a population of particles with random q/m values (charge-to-mass
   ratio) and velocities, and then tracks their motions in the static
   electric field using nearest-grid-point interpolation. 

   After each timestep the full phase space is saved
Usage:

 $MYMPIRUN -np 4 hdf5write-1 128 128 128 10000 10 NGP OOP GUARDCELLS  --processes 4 --commlib mpi

***************************************************************************/

#include "Ippl.h"

<<<<<<< .mine

=======
>>>>>>> .r6004
#include <hdf5.h>
#include "H5Part.hh"
<<<<<<< .mine
// #include "H5Block.hh"
=======
#include "H5Block.hh"
>>>>>>> .r6004


#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>


// dimension of our positions
const unsigned Dim = 3;

// some typedefs 
typedef ParticleSpatialLayout<double,Dim>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout<double,Dim> playout_t;
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
typedef IntCIC IntrplCIC_t;
typedef IntNGP IntrplNGP_t;

#define GUARDCELL 1


enum BC_t {OOO,OOP,PPP};
enum InterPol_t {NGP,CIC};

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

template<class PL>
class ChargedParticles : public ParticleBase<PL> {
public:
  ParticleAttrib<double>     Q; // charge-to-mass ratio
  typename PL::ParticlePos_t P;  // particle velocity
  typename PL::ParticlePos_t E;  // electric field at particle position
  typename PL::ParticlePos_t B;  // magnetic field at particle position

  ChargedParticles(PL* pl, BC_t bc, InterPol_t interpol, e_dim_tag decomp[Dim], bool gCells) : 
    ParticleBase<PL>(pl),
    bco_m(bc),
    interpol_m(interpol),
    fieldNotInitialized_m(true),
    doRepart_m(true),
    withGuardCells_m(gCells)
  {
    // register the particle attributes
    this->addAttribute(Q);
    this->addAttribute(P);
    this->addAttribute(E);
    this->addAttribute(B);
    setupBCs();
    for(int i=0; i<Dim; i++)
	decomp_m[i]=decomp[i];

    openDataSink(string("dataH5.dat"));

  }
  
  ~ChargedParticles() {
    H5PartCloseFile(file_m);
  }
  
  /*
    In case we have OOP or PPP boundary conditions
    we must define the domain, i.e can not be deduced from the
    particles as in the OOO case.
  */

  ChargedParticles(PL* pl, BC_t bc, InterPol_t interpol, Vector_t hr, Vector_t rmin, Vector_t rmax, e_dim_tag decomp[Dim], bool gCells) : 
    ParticleBase<PL>(pl),
    bco_m(bc),
    interpol_m(interpol),
    hr_m(hr),
    rmin_m(rmin),
    rmax_m(rmax),
    fieldNotInitialized_m(true),
    doRepart_m(true),
    withGuardCells_m(gCells)
  {
    // register the particle attributes
    this->addAttribute(Q);
    this->addAttribute(P);
    this->addAttribute(E);
    this->addAttribute(B);
    setupBCs();
    for(int i=0; i<Dim; i++)
	decomp_m[i]=decomp[i];
    openDataSink(string("dataH5.dat"));
    
  }
  
  void setupBCs() {
    if (bco_m == OOO)
      setBCAllOpen();
    else if (bco_m == PPP)
      setBCAllPeriodic();
    else 
      setBCOOP();
  }

  inline const Mesh_t& getMesh() const { return this->getLayout().getLayout().getMesh(); }
  
  inline Mesh_t& getMesh() { return this->getLayout().getLayout().getMesh(); }
    
  inline const FieldLayout_t& getFieldLayout() const {
    return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
  }

  inline FieldLayout_t& getFieldLayout() {
    return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
  }

  void gather() {
    if (interpol_m==CIC)
      gatherCIC();
    else
      gatherNGP();
  }

  double scatter() {
    Inform m("scatter ");
    double initialQ = sum(Q);
    if (interpol_m==CIC)
      scatterCIC();
    else
      scatterNGP();
    
    /*
      now sum over all gridpoints ... a bit nasty !

    */
    
    Field<double,Dim> tmpf;
    NDIndex<Dim> domain = getFieldLayout().getDomain();

    FieldLayout_t *FL  = new FieldLayout_t(getMesh(), decomp_m);
    tmpf.initialize(getMesh(), *FL);

    tmpf[domain] = EFDMag_m[domain];

    NDIndex<Dim> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> idxdom = getFieldLayout().getDomain();
    NDIndex<Dim> elem;

    double Q = 0.0;
    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
        elem[0] = Index(i,i);
        for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
	    elem[1] = Index(j,j);
	    for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
		elem[2] = Index(k,k);
		Q +=  tmpf.localElement(elem);
	    }
        }
	
    }
    reduce(Q,Q,OpAddAssign());
    m << "sum(Q)= " << initialQ << " sum(EFDMag)= " << sum(EFDMag_m) << endl;
    return initialQ-Q;
  }
  
  void myUpdate() {
    
    double hz   = hr_m[2];
    double zmin = rmin_m[2];
    double zmax = rmax_m[2];

    if (bco_m != PPP) {
      bounds(this->R, rmin_m, rmax_m);
                  
      NDIndex<Dim> domain = this->getFieldLayout().getDomain(); 
      
      for(int i=0; i<Dim; i++)
	nr_m[i] = domain[i].length();
	
      for(int i=0; i<Dim; i++)
	  hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1.0); 

      if (bco_m == OOP) {
	rmin_m[2] = zmin;
	rmax_m[2] = zmax;
	hr_m[2] = hz;
      }

      getMesh().set_meshSpacing(&(hr_m[0]));
      getMesh().set_origin(rmin_m);

      if(withGuardCells_m) {
	  EFD_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
	  EFDMag_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
      }
      else {
	  EFD_m.initialize(getMesh(), getFieldLayout(), vbc_m);
	  EFDMag_m.initialize(getMesh(), getFieldLayout(), bc_m);
      }
    } 
    else {
      if(fieldNotInitialized_m) {
	fieldNotInitialized_m=false;
	getMesh().set_meshSpacing(&(hr_m[0]));
	getMesh().set_origin(rmin_m);
	if(withGuardCells_m) {
	    EFD_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
	    EFDMag_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
	}
	else {
	    EFD_m.initialize(getMesh(), getFieldLayout(), vbc_m);
	    EFDMag_m.initialize(getMesh(), getFieldLayout(), bc_m);
	}
      }
    }
    this->update();
  }

  void initFields() {
    Inform m("initFields ");
    
    NDIndex<Dim> domain = getFieldLayout().getDomain();
    
    for(int i=0; i<Dim; i++)
      nr_m[i] = domain[i].length();
  
    int nx = nr_m[0];
    int ny = nr_m[1];
    int nz = nr_m[2];

    double phi0 = 0.1*nx;            

    m << "rmin= " << rmin_m << " rmax= " << rmax_m << " h= " << hr_m << " n= " << nr_m << endl;
    
    Index I(nx), J(ny), K(nz);
  
    assign(EFD_m[I][J][K](0), -2.0*pi*phi0/nx * cos(2.0*pi*(I+0.5)/nx) * cos(4.0*pi*(J+0.5)/ny) * cos(pi*(K+0.5)/nz));

    assign(EFD_m[I][J][K](1),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

    assign(EFD_m[I][J][K](2),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

    assign(EFDMag_m[I][J][K], 
	   EFD_m[I][J][K](0) * EFD_m[I][J][K](0) +
	   EFD_m[I][J][K](1) * EFD_m[I][J][K](1) +
	   EFD_m[I][J][K](2) * EFD_m[I][J][K](2));
  }

  Vector_t getRMin() { return rmin_m;}
  Vector_t getRMax() { return rmax_m;}
  Vector_t getHr() { return hr_m;}

  void setRMin(Vector_t x) { rmin_m = x; }
  void setHr(Vector_t x) { hr_m = x; }

  void savePhaseSpace() {

    Inform msg("savePhaseSpaceData ");
    Inform msg2all("savePhaseSpaceData ",INFORM_ALL_NODES);
    
    long long  N = this->getLocalNum();
    
    void *varray = malloc(N*sizeof(double));
    double *farray = (double*)varray;
<<<<<<< .mine
    //    h5part_int64_t *larray = (h5part_int64_t *)varray;
    long long int *larray = varray;
=======
    h5part_int64_t *larray = (h5part_int64_t *)varray;
>>>>>>> .r6004
  
    /* ------------------------------------------------------------------------ 
       Get the particle decomposition from all the nodes
    */
    long long *locN = (long long *) malloc(Ippl::getNodes()*sizeof(long long));
    long long  *globN = (long long*) malloc(Ippl::getNodes()*sizeof(long long));
    
    for(int i=0; i<Ippl::getNodes(); i++) {
      globN[i] = locN[i]=0;
    }
    locN[Ippl::myNode()] = N;
    reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());
    
    /* ------------------------------------------------------------------------ */
    
    double actPos           = 0.0;
    double structLenght     = 0.0;
    Vector_t org(0.0);
    Vector_t maxX(0.0);
    Vector_t minX(0.0);
    Vector_t maxP(0.0);
    Vector_t minP(0.0);
    unsigned long protons   = N;
    unsigned long electrons = 0;
    Vector_t centroid(0.0);
    unsigned nTot           = this->getTotalNum();

    H5PartSetStep(file_m,idx_m);  
    H5PartSetNumParticles(file_m,N); 
  
    /* write scalar data i.e the header */
    long long step = idx_m;
    H5PartWriteStepAttrib(file_m,"Step", H5T_NATIVE_INT64,&step,1);
  
    H5PartWriteStepAttrib(file_m,"Spos",     H5T_NATIVE_DOUBLE,&actPos,1);
    H5PartWriteStepAttrib(file_m,"structLen",H5T_NATIVE_DOUBLE,&structLenght,1);  
    H5PartWriteStepAttrib(file_m,"org",      H5T_NATIVE_DOUBLE,&org,3);
    H5PartWriteStepAttrib(file_m,"maxX",     H5T_NATIVE_DOUBLE,&maxX,3);
    H5PartWriteStepAttrib(file_m,"minX",     H5T_NATIVE_DOUBLE,&minX,3);
    H5PartWriteStepAttrib(file_m,"maxP",     H5T_NATIVE_DOUBLE,&maxP,3);
    H5PartWriteStepAttrib(file_m,"minP",     H5T_NATIVE_DOUBLE,&minP,3);
    H5PartWriteStepAttrib(file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);
  
    H5PartWriteStepAttrib(file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

    for (long long i=0; i<N;i++)
      farray[i] =  this->R[i](0);
    H5PartWriteDataFloat64(file_m,"x",farray); 
  
    for (long long i=0; i<N;i++)
      farray[i] =  this->R[i](1);
    H5PartWriteDataFloat64(file_m,"y",farray); 

    for (long long i=0; i<N;i++)
      farray[i] =  this->R[i](2);
    H5PartWriteDataFloat64(file_m,"z",farray); 
    
    for (long long i=0; i<N;i++)
      farray[i] =  P[i](0);
    H5PartWriteDataFloat64(file_m,"px",farray); 

    for (long long i=0; i<N;i++)
      farray[i] =  P[i](1);
    H5PartWriteDataFloat64(file_m,"py",farray); 

    for (long long i=0; i<N;i++)
      farray[i] =  P[i](2);
    H5PartWriteDataFloat64(file_m,"pz",farray); 
  
    for (long long i=0; i<N;i++) {
      if (Q[i] > 0.0) 
	larray[i] =  this->ID[i];
      else           
	larray[i] =  -1*(long long int)this->ID[i];
    }
    H5PartWriteDataInt64(file_m,"id",larray);  
<<<<<<< .mine

    // ada save block data
    /*  
  NDIndex<Dim> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem;
    h5part_int64_t herr = H5BlockDefine3DFieldLayout (
				       file_m,
				       idx[0].min(), idx[0].max(),
				       idx[1].min(), idx[1].max(),
				       idx[2].min(), idx[2].max());

    h5part_float64_t *data;
    data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );
    int ii = 0;
    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
      elem[0] = Index(i,i);
      for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
	elem[1] = Index(j,j);
	for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
	  elem[2] = Index(k,k);
	  data[ii] = EFDMag_m.localElement(elem);
	  ii++;
	}
      }
    } 
    herr = H5Block3dWriteScalarField ( file_m, "EFmag", data );
    
    if(data)
      free(data);

*/

    if(varray)  
=======

    // ada save block data

    h5part_int64_t l[6];
    
    NDIndex<Dim> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> elem;
    h5part_int64_t herr = H5BlockDefine3DFieldLayout (
				       file_m,
				       idx[0].min(), idx[0].max(),
				       idx[1].min(), idx[1].max(),
				       idx[2].min(), idx[2].max());
    if (herr < 0)
      msg << "H5BlockDefine3DFieldLayout err " << herr << endl;

    h5part_float64_t *data;
    data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );
    int ii = 0;
    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
      elem[0] = Index(i,i);
      for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
	elem[1] = Index(j,j);
	for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
	  elem[2] = Index(k,k);
	  data[ii] = EFDMag_m.localElement(elem);
	  ii++;
	}
      }
    } 

    herr = H5Block3dWriteScalarField ( file_m, "EFmag", data );
    if (herr < 0)
      msg << "H5Block3dWriteScalarField err " << herr << endl;

    msg << "Wrote Block data i.e. scalar field: sum(f)= " << sum(EFDMag_m) << endl;

    if (Ippl::myNode() == 0) {
      for (h5part_int64_t p=0; p<Ippl::getNodes();p++) {
	herr = H5Block3dGetPartitionOfProc(file_m, p, &l[0], &l[1], &l[2], &l[3], &l[4], &l[5]);
	stringstream lstr;
	lstr << "layout" << p;
	H5BlockWriteFieldAttrib (file_m,"EFmag", lstr.str().c_str(), H5PART_INT64,l,6);       
      }
    }


    if(data)
      free(data);
  if(varray)  
>>>>>>> .r6004
      free(varray);
    
    idx_m++;
  }

private:

  void openDataSink(string hdf5fn) {

#ifdef PARALLEL_IO
    file_m=H5PartOpenFileParallel(hdf5fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
    file_m=H5PartOpenFile(hdf5fn.c_str(),H5PART_WRITE);
#endif
    if(!file_m) {
      INFOMSG("File open failed:  exiting!" << endl);
      exit(0);
    }   

    idx_m = 0;
  }

  inline void setBCAllOpen() {
    for (int i=0; i < 2*Dim; i++) {
      this->getBConds()[i] = ParticleNoBCond;
      bc_m[i]  = new ZeroFace<double  ,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
    }
  }
  
  inline void setBCAllPeriodic() {
    for (int i=0; i < 2*Dim; i++) {
      this->getBConds()[i] = ParticlePeriodicBCond;
      bc_m[i]  = new PeriodicFace<double  ,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
    }  
  }

  inline void setBCOOP() {
    for (int i=0; i < 2*Dim - 2; i++) {
      bc_m[i]  = new ZeroFace<double  ,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
      this->getBConds()[i] = ParticleNoBCond;
    }
    for (int i= 2*Dim - 2; i < 2*Dim; i++) {
      bc_m[i]  = new PeriodicFace<double  ,Dim,Mesh_t,Center_t>(i);
      vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
      this->getBConds()[i] = ParticlePeriodicBCond;    
    }  
  }

  inline void gatherNGP() {
    // create interpolater object (nearest-grid-point method)
    IntNGP myinterp;
    E.gather(EFD_m, this->R, myinterp);
  }

  inline void gatherCIC() {
    // create interpolater object (cloud in cell method)
    IntCIC myinterp;
    E.gather(EFD_m, this->R, myinterp);
  }

  inline void scatterCIC() {
    // create interpolater object (cloud in cell method)
    IntCIC myinterp;
    Q.scatter(EFDMag_m, this->R, myinterp);
  }

  inline void scatterNGP() {
    // create interpolater object (cloud in cell method)
    IntNGP myinterp;
    Q.scatter(EFDMag_m, this->R, myinterp);
  }

  Field<Vektor<double,Dim>,Dim> EFD_m;
  Field<double,Dim> EFDMag_m;
  
  BConds<double,Dim,Mesh_t,Center_t> bc_m;
  BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;

  Vektor<int,Dim> nr_m;
    
  Vector_t hr_m;
  Vector_t rmin_m;
  Vector_t rmax_m;

  BC_t bco_m;
  InterPol_t interpol_m;
  bool fieldNotInitialized_m;
  bool doRepart_m;
  bool withGuardCells_m;  
  e_dim_tag decomp_m[Dim];
  
#ifdef GTHDF5
  H5PartFile *file_m;
#endif

  int idx_m;
};


int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);

  H5PartSetVerbosityLevel(0);

  Vektor<int,Dim> nr(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

  const unsigned int totalP = atoi(argv[4]);
  const int nt              = atoi(argv[5]);

  InterPol_t myInterpol;
  if (string(argv[6])==string("CIC"))
    myInterpol = CIC;
  else
    myInterpol = NGP;
  
  bool gCells;
  gCells =  (string(argv[8])==string("GUARDCELLS"));
  

  msg << "Particle test PIC3d " << endl;
  msg << "nt " << nt << " Np= " << totalP << " grid = " << nr <<endl;
  
  if (myInterpol==CIC)
    msg << "Cloud in cell (CIC) interpolation selected" << endl; 
  else
    msg << "Nearest grid point (NGP) interpolation selected" << endl; 
  
  if (gCells)
      msg << "Using guard cells" << endl;
  else
      msg << "Not using guard cells" << endl;

  BC_t myBC;  
  if (string(argv[7])==string("OOO")) {
    myBC = OOO; // open boundary 
    msg << "BC == OOO" << endl;
  }
  else if (string(argv[7])==string("OOP")) {
    myBC = OOP; // open boundary in x and y, periodic in z
    msg << "BC == OOP" << endl;
  }
  else {
    myBC = PPP; // all periodic 
    msg << "BC == PPP" << endl;
  }

  e_dim_tag decomp[Dim];
  int serialDim = 2;

  msg << "Serial dimension is " << serialDim  << endl;

  Mesh_t *mesh;
  FieldLayout_t *FL;
  ChargedParticles<playout_t>  *P;

  NDIndex<Dim> domain;
  if (gCells) {
      for(int i=0; i<Dim; i++)
	  domain[i] = domain[i] = Index(nr[i] + 1);
  }
  else {
      for(int i=0; i<Dim; i++)
	  domain[i] = domain[i] = Index(nr[i]);
  }

  for (int d=0; d < Dim; ++d)
      decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;
  
  // create mesh and layout objects for this problem domain
  mesh          = new Mesh_t(domain);
  FL            = new FieldLayout_t(*mesh, decomp);
  playout_t* PL = new playout_t(*FL, *mesh);

  msg << "Mesh, FL and playout created " << endl;

  if (myBC==OOO)
    P = new ChargedParticles<playout_t>(PL,myBC,myInterpol,decomp,gCells);
  else {
    /*
      In case of periodic BC's define
      the domain with hr and rmin
    */

    Vector_t hr(1.0);
    Vector_t rmin(0.0);
    Vector_t rmax(nr);

    P = new ChargedParticles<playout_t>(PL,myBC,myInterpol,hr,rmin,rmax,decomp,gCells);
  }


  // initialize the particle object: do all initialization on one node,
  // and distribute to others
  P->create(totalP / Ippl::getNodes());
    
  msg << "particles created " << endl;
    
  // quiet start for particle positions
  assign(P->R(0), IpplRandom * nr[0]);
  assign(P->R(1), IpplRandom * nr[1]);
  assign(P->R(2), IpplRandom * nr[2]);

  double q = 1.0/totalP;

  // random initialization for charge-to-mass ratio
  assign(P->Q,q);
    
  msg << "initial conditions assigned " << endl;

  // redistribute particles based on spatial layout
  P->myUpdate(); 

  msg2all << "Nlocal= " << P->getLocalNum() << endl;
  Ippl::Comm->barrier();

  msg << "initial update and initial mesh done .... Q= " << sum(P->Q) << endl;
  
  msg << P->getMesh() << endl;
  msg << P->getFieldLayout() << endl;

  msg << "scatter test done delta= " <<  P->scatter() << endl;
    
  P->initFields();
  msg << "P->initField() done " << endl;

  Timer mytimer;

  // begin main timestep loop
  msg << "Starting iterations ..." << endl;
  for (unsigned int it=0; it<nt; it++) {
    // advance the particle positions
    // basic leapfrogging timestep scheme.  velocities are offset
    // by half a timestep from the positions.
    assign(P->R, P->R + dt * P->P);
    
    // update particle distribution across processors
    P->myUpdate();
    
    // gather the local value of the E field
    P->gather();
      
    // advance the particle velocities
    assign(P->P, P->P + dt * P->Q * P->E);
    msg << "Finished iteration " << it << " - min/max r and h " << P->getRMin() << P->getRMax() << P->getHr() << endl;
    mytimer.clear();    
    mytimer.start();
    P->savePhaseSpace();
    mytimer.stop();

    long long  N = P->getLocalNum();
    double rate = ((6*N*sizeof(double) + N*sizeof(unsigned))/1000000.0) / mytimer.clock_time(); 
    reduce(rate,rate,OpAddAssign());      
    reduce(N,N,OpAddAssign());      
    msg << "Number of particles (x,y,z,px,py,pz,id) " << N << " in file set Write to disk bw= " << rate << " [MB/s] "<< endl; 
  }
  
  Ippl::Comm->barrier();
  msg << "Particle test PIC3d: End." << endl;

  delete P;

  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

