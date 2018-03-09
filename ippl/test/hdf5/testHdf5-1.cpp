/*
  testHdf5-1.cpp 

  Ok just to see if hdf5 does not interfere with all 
  the other stuff I compile and link with

  The HDF5 stuff is based on some cool "jump start help" by John
  Shalf. Only serial and with some ugly buffer copy. I dump all 
  random generated particles at this time just to keep thing busy.

  Fri Mar 14 (the windiest day in 2003 so far ... and I could not go out)
  
*/
#include <fstream>
#include <new>
#include <exception>
#include <vector>
#include <map>
#include <set> 
#include <iostream>
#include <sstream>

#include "Ippl.h"
#include "extpde.h"
#include "H5Part.hh"

using namespace std;

typedef Vektor<double,3> Vector_t;

class ChargedParticles: public IpplParticleBase<double,3> {
private:
  typedef IpplParticleBase<double,3>::ParticlePos_t ParticlePos_t;

  Grid *grid_m;

  H5PartFile *file;

  long long *hdf5BuffID;
  double *hdf5BuffX;
  double *hdf5BuffY;
  double *hdf5BuffZ;
  

public:
  ParticlePos_t P;
   
  ChargedParticles(Grid *g) :
    grid_m(g)
  { 
    // particle base knows all about the domain
    addAttribute(grid_m);
    
    // add more variables
    addAttribute(P);

    file=H5PartOpenFile("parttest.h5",H5PART_WRITE);
    if(!file) {
      perror("File open failed:  exiting!");
      exit(0);
    }
  }

  void allocateHdf5Buffers() {
    hdf5BuffID = (long long *)malloc(getLocalNum()*sizeof(long long));
    hdf5BuffX = (double*)malloc(getLocalNum()*sizeof(double));
    hdf5BuffY = (double*)malloc(getLocalNum()*sizeof(double));
    hdf5BuffZ = (double*)malloc(getLocalNum()*sizeof(double));
  }

  void dump2Hdf5(unsigned int step) {
    
    for(unsigned i=0;i<getLocalNum();i++) {
      hdf5BuffX[i]=R[i](0);
      hdf5BuffY[i]=R[i](1);
      hdf5BuffZ[i]=R[i](2);
      hdf5BuffID[i]=ID[i];
    }
    H5PartSetStep(file,step); /* must set the current timestep in file */
    H5PartSetNumParticles(file,getLocalNum()); /* then set number of particles to store */
    /* now write different tuples of data into this timestep of the file */
    H5PartWriteDataFloat64(file,"x",hdf5BuffX); 
    H5PartWriteDataFloat64(file,"y",hdf5BuffY);
    H5PartWriteDataFloat64(file,"z",hdf5BuffZ);
    H5PartWriteDataInt64(file,"id",hdf5BuffID);
  }


  ~ChargedParticles() 
  {
    H5PartCloseFile(file);
    delete(hdf5BuffID);
    delete(hdf5BuffZ);
    delete(hdf5BuffY);
    delete(hdf5BuffX);
  }
};


int main(int argc, char ** argv)
{

  int mynodes,myrank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mynodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  Ippl ippl(argc, argv);
  Inform msg2all("testPart-2 ",INFORM_ALL_NODES);  
  Inform msg("testPart-2 ");

  IpplRandom.SetSeed(static_cast<unsigned long>(23131719*(Ippl::myNode()+13)));

 /*
    read in from paramX.inp
  */
  double a; double b; double l;   // geometry
  double maxErr;
  int iter_relax;
  double o_square, s_square; 
  int p_test;
  int r_parallel;
  long nPart;
  double normalConst;
  double r,hInit;
  
  ifstream PARAMETER;
  string filename(argv[1]);
  PARAMETER.open(filename.c_str(),std::ios::in);
  PARAMETER >> hInit  >> maxErr >> iter_relax >> r_parallel 
            >> o_square >> s_square >> p_test >> a >> b 
	    >> l >> nPart >> r >> normalConst;
  PARAMETER.close();


  unsigned long nPartPerNode = 
    static_cast<unsigned long>(nPart/Ippl::getNodes());
  unsigned long pDelta = nPart - nPartPerNode*Ippl::getNodes();
  if (Ippl::myNode() == 0) { // ajust on one node
    if (nPart < 0)
      nPartPerNode += -1*pDelta;
    else
      nPartPerNode -= pDelta;
  }
  
  msg << "Running test number " << p_test << endl;

  Grid_gen_parameters ggen_parameter;
  ggen_parameter.Set_r_parallel(r_parallel);    
  ggen_parameter.Set_offset_square(o_square);
  ggen_parameter.Set_stretch_square(s_square);          
  
  /*
    Domain is a cylinder stretched 
    by Vs and translated by Vm
  */ 

  D3vector Vs(a,b,l);
  D3vector Vm(-a/2.0,-b/2.0,-l/2.0);

  Cylinder   cylinder;
  Ball       ball;
  Square     square;
  //  Domain geomDomain(&cylinder*Vs+Vm);
 
  // Domain geomDomain(&square*Vs+Vm); 
  Domain geomDomain(&square);
  Grid *grid;

  Variable *rhs;
  Variable *tmp;
  Variable *f;
  Variable *err;
  Variable *u;
  Res_stencil_boundary *ResS;
    
  if (hInit>0.0)
    grid = new Grid(hInit,&geomDomain,ggen_parameter,MPI_COMM_WORLD); // Construction of grid
  else {
    int nMax = (int) abs(hInit);
    grid = new Grid(nMax,&geomDomain,ggen_parameter,MPI_COMM_WORLD); // Construction of grid
  }
  
  rhs  = new Variable(grid);
  tmp  = new Variable(grid);
  f    = new Variable(grid);
  err  = new Variable(grid);
  u    = new Variable(grid);  
  ResS = new Res_stencil_boundary(grid);           // definition of a stencil 
    
  grid->Initialize();

  /*
    Now the particle stuff 
  */

  ChargedParticles *bunch = new ChargedParticles(grid);
  
  // create uniformly distributed particles in the computational domain

  for(unsigned p=0; p<nPartPerNode; p++) {
    Vektor<double,3> x = Vektor<double,3>(a*(-0.5 + IpplRandom()),
					  b*(-0.5 + IpplRandom()),
					  l*(-0.5 + IpplRandom()));
    if (bunch->create(x,1.0)) {
      bunch->R[p] = x;
      bunch->P[p] = Vektor<double,3>(Ippl::myNode()+1);
    }
    else 
      p--;
  }
  

  bunch->update();
  msg << "max(R)= " << max(bunch->R) << " min(R)= " << min(bunch->R) << endl;
  msg << "Ntotal= " << bunch->getTotalNum() << " sum(P)= " << sum(bunch->P) << endl;
  msg2all << "Nlocal= " << bunch->getLocalNum() << endl;
 
  Vector_t s;
  for(unsigned i=0; i<bunch->getLocalNum(); i++)
    s+=bunch->P[i];
  msg << "Ntotal= " << bunch->getTotalNum() << " sum2(P)= " << s << endl;
  
  
  bunch->allocateHdf5Buffers();
  bunch->dump2Hdf5(0);
    

 /*
     f = Helm_FE(bunch,bunch->q);

     for(int i=0; i<bunch->totoalNum(); i++) {
      
     Index3D cellCoord = getCell(bunch->R[i]);
       if (validCell(cellCoord) { 
         thet = getThet(cellCoord);
	 addContib(bunch->q[i],thet, f)

       }
     }

     if (grid.Give_my_index() <= cellCoord)

     gridbase.h, Give_cell_type(cellCoord)

     bunch->grid_m->Give_variable(I.neighbour(WSTd),maxLevel)[f.Number_Variable()]=+alpha*bunch->q[i];

     alpha, barycentric variable von WSTd in I. p45 FEM

  */

  msg << "All done bye bye ...... " << endl;
  Ippl::Comm->barrier();
  return 0;
}

