/***************************************************************************
 
Units:
 
     R, position in units of mega paresec
     V, velocity in units of redshift
     F, force fiels in units of ???
 
 
Usage: cosmo <inputfile> [<distrInputFile>]
 
 
     mpirun -np 2 cosmo <inputfile> [<distrInputFile>] --commlib mpi

And in real:

     mpirun -np 1 cosmo indat fort.66 --commlib mpi

And for lesser output try:

     mpirun -np 1 cosmo indat fort.66 --info 0--commlib mpi

***************************************************************************/

#include "ChargedParticles.hh"

#include <ctime>

#include "FFTPoissonSolverPeriodic.hh"
#include <iostream>
#include <sstream>

#include "Parser/FileStream.h"
#include "SimData.hh"
#include "PwrSpec.hh"

using namespace std;

#define T double

/************************************/
 const string version = "1.1";
/************************************/

const unsigned int Dim = 3;


typedef ChargedParticles<T,Dim>::Mesh_t Mesh_t;
typedef ChargedParticles<T,Dim>::Layout_t Layout_t;
typedef ChargedParticles<T,Dim>::Vector_t Vector_t;  
typedef ChargedParticles<T,Dim>::FieldLayout_t FieldLayout_t;


int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg ("cosmos ");
  
  const double Pi = 4.0*atan(1.0);

  /// here we read in all data from the input
  /// file and compute simple derifed quantities
  
  SimData simData(argc, argv);
  msg << simData << endl;
  
  // create layout objects
  Index I(simData.ng), J(simData.ng), K(simData.ng);
  NDIndex<Dim> domain(I,J,K);
  
  // For efficiency in the FFT's, we will use a parallel decomposition
  // which is serial in the first dimension.
  int serialDim = 0;
  e_dim_tag decomp[Dim];
  for (int d=0; d < Dim; ++d)
    decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;
  
  // create mesh and layout objects for this problem domain
  Mesh_t *mesh;
  FieldLayout_t *FL;
  mesh = new Mesh_t(domain);
  FL = new FieldLayout_t(*mesh, decomp);
  Layout_t * PL = new Layout_t(*FL, *mesh); 

  // create an empty ChargedParticles object and set the BC's
  ChargedParticles<T,Dim>* univ = new ChargedParticles<T,Dim>(PL);

  PwrSpec<T,Dim>* pSpec = new PwrSpec<T,Dim>(univ);
  
  if (simData.jinit == 1)
    univ->readInputDistribution(string(argv[2]));
  else
    univ->genIC(simData);


  // SOME GENERAL INITIAL CONDITION SETTINGS
  univ->setTime(0.0);
  
  univ->set_FieldSolver(new FFTPoissonSolverPeriodic<T,Dim>(univ));
  
  int ng        = simData.ng/2 + 1;
  double gpscal = (1.0*ng)/(1.0*simData.ng);  // will be needed to scale IC's to the grid

  double ain  = 1.0 /(1+simData.zin);
  double afin = 1.0 /(1+simData.zfin);
  double apr  = 1.0 /(1+simData.zpr);

  double omegab   = simData.deut / simData.hubble / simData.hubble;
  double omegatot = simData.omegadm + omegab;

  /// Convert to `symplectic velocities' at a=a_in 
  univ->V *= Vector_t(ain*ain,ain*ain,ain*ain);
  
  double pp   = pow(ain ,simData.alpha);
  double pfin = pow(afin,simData.alpha);
  double ppr  = pow(apr ,simData.alpha);
 
  double tau = (pfin-pp)/(simData.nsteps);
  
  double pmass = 2.77536627e11*pow(simData.rL,3.0)*(omegatot)/simData.hubble/univ->getTotalNum();

  msg << "Total Mass = " << pmass << endl;

  if(simData.pbchoice==2) {
    /// input date is in physical units rescale to box units
    
    univ->R *= Vector_t(1.0/simData.rL*(simData.ng));
    univ->V *= Vector_t(1.0/simData.rL*(simData.ng)/100);

    univ->rescaleAll();
    msg << "We are in box units now: min(R) " << min(univ->R) << " max(R)= " << max(univ->R) << endl;
    msg << "We are in box units now: min(V) " << min(univ->V) << " max(V)= " << max(univ->V) << endl;
  }
  
  univ->dumpParticles(string("InitialPhaseSpaceBoxUnits.dat"));  
  
  // / initial power spectra is calculated and saved
  
  // / Initialize mass per simulation particle
  
  univ->M = 1.0;
  univ->rho_m = 0.0;
  univ->M.scatter(univ->rho_m, univ->R, IntCIC());
 
  msg << "read rho field from fort13.dat " << endl; 
  ifstream ifstr;
  univ->rho_m = 0.0;
  if (Ippl::myNode() == 0) {
    ifstr.open(string("fort.13").c_str(),ios::in );
    ifstr.precision(8);
    int dummy;
    NDIndex<3> elem;
    for (int i=domain[0].min(); i<=domain[0].max(); ++i) {
      elem[0]=Index(i,i);
      for (int j=domain[1].min(); j<=domain[1].max(); ++j) {
	elem[1]=Index(j,j);
	for (int k=domain[2].min(); k<=domain[2].max(); ++k) {
	  elem[2]=Index(k,k);
	  ifstr >> dummy >> dummy >> dummy >> univ->rho_m.localElement(elem);
	}
      }
    }
    ifstr.close();
  }

  msg << "max(rho)= " << max(univ->rho_m) 
      << " min(rho)= " << min(univ->rho_m)
      << " sum(rho)= " << sum(univ->rho_m) << endl;

  
  msg << "write rho field to rho.dat " << endl; 
  ofstream o;
  if (Ippl::myNode() == 0) {
    o.precision(8);
    o.open(string("rho.dat").c_str()); 
    o.setf(ios_base::scientific);
    
    for (int i=domain[0].min(); i<=domain[0].max(); ++i)
      for (int j=domain[1].min(); j<=domain[1].max(); ++j)
	for (int k=domain[2].min(); k<=domain[2].max(); ++k)
	  o << i << "  " << j << "  " << k << "  " << univ->rho_m[i][j][k] << endl;          
    o.close();
  }

  univ->rho_m = univ->rho_m-1.0;

  /* 
  NEEEDS WORK 
  */

  pSpec->cacl1D(univ);
  msg << "1DSpectra norm: tpiL= " << 2*Pi/simData.rL << " rscale= " << simData.rL/simData.ng << endl;
  pSpec->save1DSpectra(string("InitialPwrSpec.dat")); //,2*Pi/simData.rL ,simData.rL/simData.ng);

   
  univ->testConvolution();

  msg << *univ;
  
  double at=   pow(pp,1.0/simData.alpha);   
  double zact= 1.0/pow(pp,1.0/simData.alpha) - 1.0; 

  msg << endl << endl << "Start Integration " << endl;
  msg << "a(t)= \t" << at << " z= \t" << zact << " step= \t" << 0 << endl;
  for (int i=0; i<simData.nsteps; i++){    

    //univ->rtmp = univ->R;
                                
    univ->map1(0.5*tau,omegatot,simData.alpha,pp);       // R = R + t*V
    pp += 0.5*tau;

    //univ->map2(tau,omegatot,simData.alpha,pp);           // V update with old R (rtmp)
    pp += 0.5*tau; 

    univ->map1(0.5*tau,omegatot,simData.alpha,pp);       // last stream with R

    at  = pow(pp,1.0/simData.alpha);
    zact= 1.0/pow(pp,1.0/simData.alpha) - 1.0; 

    if(simData.entest == 0) {
      
      // scatter R -> rho

      // FFT rho

      // bin to get 1D spectrum

    }
    
    if (!(i % simData.iprint)) {
      string fn;
      char numbuf[6];
      sprintf(numbuf, "%05d", i);  
      fn = string("phasespace-") + string(numbuf) + string(".dat");
      univ->dumpParticles(fn);  
    }
    INFOMSG("a(t)= \t" << at << " z= \t" << zact << " step= \t" << i << endl);
  }
  univ->dumpParticles(string("FinalPhaseSpaceBoxUnits.dat"));  
  Ippl::Comm->barrier();
  msg << "All nodes done dude ..." << endl;
  return 0;
}









