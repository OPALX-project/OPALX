// -*- C++ -*-
/***************************************************************************
 *
 * 
 * FFTPoissonSolverPeriodic.cc
 * 
 * 
 * 
 *
 * 
 * 
 *
 ***************************************************************************/
// include files
#include <iostream>
using namespace std;


////////////////////////////////////////////////////////////////////////////
// constructor
template <class T, unsigned int Dim>
FFTPoissonSolverPeriodic<T,Dim>::FFTPoissonSolverPeriodic(ChargedParticles<T,Dim> *beam):
  layout_m(&(beam->getFieldLayout())),
  mesh_m(&(beam->getMesh()))
{
  gDomain_m = layout_m->getDomain();        // global domain
  lDomain_m = layout_m->getLocalNDIndex(); // local domain
  // create additional objects for the use in the FFT's
  green_m.initialize(*mesh_m, *layout_m);
  rho_m.initialize(*mesh_m, *layout_m);
  // create the FFT object
  bool compressTemps = true;
  fft_m = new FFT_t(layout_m->getDomain(), compressTemps);
  fft_m->setDirectionName(+1, "forward");
  fft_m->setDirectionName(-1, "inverse");
  greensFunction();
  saveGreensHatFunction();
}


////////////////////////////////////////////////////////////////////////////
// destructor
template <class T, unsigned int Dim>
FFTPoissonSolverPeriodic<T,Dim>::~FFTPoissonSolverPeriodic()
{
  // delete the FFT object
  delete fft_m;
}


////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation
template<class T, unsigned int Dim>
void FFTPoissonSolverPeriodic<T,Dim>::computeForceField(ChargedParticles<T,Dim> *beam)
{
  // deposit charge/mass on grid.  We first zero out the rho field, then
  // scatter particles using our selected interpolater type.  In this
  // scatter routine, we cache the values used to calculate where the
  // particles go in the field, and use these values later on in this
  // spaceCharge function to gather elements back from the field to the
  // particle positions.  We can do that because the particles do not
  // change their position between the scatter and gather phases.
  
  beam->rho_m = 0.0;
  beam->F     = Vector_t(0.0);

  beam->M.scatter(beam->rho_m, beam->rtmp, IntrplCIC_t());

  // rho_m is 
  rho_m = beam->rho_m;

  // FFT density field
  fft_m->transform( "forward", rho_m);
  
  // multiply transformed charge density
  // and transformed Green function
  rho_m *= green_m;
  
  // inverse FFT
  fft_m->transform("invers", rho_m);
  
  // rho is now phi, the scalar potential
  beam->rho_m = real(rho_m);

  // finally, compute the gradient of the scalar field 
  beam->gphi_m = -Grad(beam->rho_m, beam->gphi_m);
  
  // interpolate the field at particle positions.  We reuse the
  // cached information about where the particles are relative to the
  // field, since the particles have not moved since this the most recent
  // scatter operation.
  beam->F.gather(beam->gphi_m, beam->rtmp, IntrplCIC_t() );
}

template<class T, unsigned int Dim>
void FFTPoissonSolverPeriodic<T,Dim>::testConvolution(ChargedParticles<T,Dim> *beam)
{
  // assume field rho_m is valid
  Inform msg("testConvolution ");

  rho_m = beam->rho_m;
  
  // FFT density field
  fft_m->transform( "forward", rho_m);
  
  // multiply transformed charge density
  // and transformed Green function
  rho_m *= green_m;
  
  double scale = beam->getGridSize()(0)*beam->getGridSize()(1)*beam->getGridSize()(2);
  msg << "max(real(rho))= " << max(real(rho_m))*scale << " min(real(rho))= " << min(real(rho_m))*scale << endl; 
  msg << "max(imag(rho))= " << max(imag(rho_m))*scale << " min(imag(rho))= " << min(imag(rho_m))*scale << endl; 
}




///////////////////////////////////////////////////////////////////////////
// calculate the FFT of the Green's function for the given field
template<class T, unsigned int Dim>
void FFTPoissonSolverPeriodic<T,Dim>::greensFunction()
{
  /*
    Create Greens function for periodic 1/r potential
  */  
  NDIndex<3> elem;
  double pi = acos(-1.0);
  double tpinx = 2.0*pi/(1.0*gDomain_m[0].max()+1);
  double tpiny = 2.0*pi/(1.0*gDomain_m[1].max()+1);
  double tpinz = 2.0*pi/(1.0*gDomain_m[2].max()+1);
  
  green_m = dcomplex(0.0);

  for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
    elem[0]=Index(i,i);
    for (int j=lDomain_m[1].min(); j<=lDomain_m[1].max(); ++j) {
      elem[1]=Index(j,j);
      for (int k=lDomain_m[2].min(); k<=lDomain_m[2].max(); ++k) {
        elem[2]=Index(k,k);
        if (i+j+k != 0) {
	  double val = 0.5 / (cos(i*tpinx)+cos(j*tpiny)+cos(k*tpinz)-3.0);
          green_m.localElement(elem) = dcomplex(val);
	}
      }
    }
  }
}

template<class T, unsigned int Dim>
void FFTPoissonSolverPeriodic<T,Dim>::saveGreensHatFunction()
{
  /*
    Saves  Greens function works only serial
  */
  ofstream o;
  NDIndex<3> elem;
  string fn("Ghat.dat");
  if (Ippl::myNode() == 0) {
    o.precision(20);
    o.open(fn.c_str());
    o.setf(ios_base::scientific);
    
    for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
      elem[0]=Index(i,i);
      for (int j=lDomain_m[1].min(); j<=lDomain_m[1].max(); ++j) {
	elem[1]=Index(j,j);
	for (int k=lDomain_m[2].min(); k<=lDomain_m[2].max(); ++k) {
	  elem[2]=Index(k,k);
	  o << i << "  " << j << "  " << k << " \t " <<  real(green_m.localElement(elem)) << endl;
	}
      }
    }
    o.close();
  }
}









/***************************************************************************
 * $RCSfile: FFTPoissonSolverPeriodic.cc,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2001/08/16 09:36:09 $
 ***************************************************************************/




/*
  if ( Pooma::getNodes() == 1) {
  // Calculate self energy of charge distribution
  NDIndex<Dim> elem;
  for (unsigned int i = 0 ; i < domain_m[0].length(); i++) {
  elem[0] = Index(i,i);
  for (unsigned int j = 0 ; j < domain_m[1].length(); j++) {
  elem[1] = Index(j,j);
  for (unsigned int k = 0 ; k < domain_m[2].length(); k++) {
  elem[2] = Index(k,k);
  beam->U_m += 0.5*beam->getcoupling()*beam->rho_m.localElement(elem)*rho2_m.localElement(elem);
  }
  }
  }
  }
  else
  beam->U_m = -1.0;
*/
  
