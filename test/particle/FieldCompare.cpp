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
 * FieldCompare -
 *   Compare the performance of expression evaluation for a 1D Field and for
 * particle objects with the same number of elements.
 ***************************************************************************/

#include "Ippl.h"

// A simple particle object using an InteractLayout
template<class T, unsigned Dim>
class InteractAtoms : public ParticleBase< ParticleInteractLayout<T, Dim> > {

public:
  // the attributes for this set of particles (atoms).
  ParticleInteractAttrib<T> A;
  ParticleInteractAttrib<T> B;

  // the constructor: we need FieldLayout for the particles
  InteractAtoms(FieldLayout<Dim>& fl)
    : ParticleBase< ParticleInteractLayout<T,Dim> >(
			    new ParticleInteractLayout<T,Dim>(fl)) {
      // register our attributes with the base class
      addAttribute(A);
      addAttribute(B);
  }
};


// A simple particle object using a SpatialLayout
template<class T, unsigned Dim>
class SpatialAtoms : public ParticleBase< ParticleSpatialLayout<T, Dim> > {

public:
  // the attributes for this set of particles (atoms).
  ParticleAttrib<T> A;
  ParticleAttrib<T> B;

  // the constructor: we need FieldLayout for the particles
  SpatialAtoms(FieldLayout<Dim>& fl)
    : ParticleBase< ParticleSpatialLayout<T,Dim> >(
			    new ParticleSpatialLayout<T,Dim>(fl)) {
      // register our attributes with the base class
      addAttribute(A);
      addAttribute(B);
  }
};


// the main program
int main(int argc, char *argv[]) {
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);

  // find out how many elements we should use ... we will make sure to only
  // have one vnode per processor
  int elements = 0;
  if (argc != 2 || (elements = atoi(argv[1])) < 1) {
    msg << "Usage: " << argv[0] << " <elements>" << endl;
    exit(1);
  }

  // create a 1D Field with no bc's or gc's, and initialize all the values
  Index I(elements);
  FieldLayout<1U> layout(I);
  Field<double,1U> F(layout);
  Field<double,1U> F2(layout);
  Field<double,1U> F3(layout);
  F[I] = I;
  F2[I] = I*I;
  F3[I] = I*I*I;

  // create particle objects
  SpatialAtoms<double,1U>  SA(layout);
  InteractAtoms<double,1U> IA(layout);

  // put in particles with same number on each processor as size of local
  // lfield
  BareField<double,1U>::iterator_if lf_i = F.begin_if();
  LField<double,1U> &lf = *((*lf_i).second);
  const NDIndex<1U>& lo = lf.getOwned();
  msg << "Creating " << lo.size() << " local particles." << endl;
  SA.create(lo.size());
  IA.create(lo.size());

  // initialize values of particles and do an update to get all particles
  // in the right place
  SA.A = (double)(lo[0].first());
  SA.B = (double)(lo[0].last());
  IA.A = (double)(lo[0].first());
  IA.B = (double)(lo[0].last());
  SA.update();
  IA.update();
  msg << "After particle update: localsize(SA)=" << SA.getLocalNum();
  msg << ", locasize(IA)=" << IA.getLocalNum() << endl;

  // now time operations involving the three quantities
  double clock1;
  int i;
  Timer timer;

  timer.clear(); timer.start();
  F[I] = F2[I] - F3[I];
  timer.stop(); clock1 = timer.cpu_time();
  msg << "Time for Field operation = " << clock1 << endl;

  timer.clear(); timer.start();
  SA.A = SA.A - SA.B;
  timer.stop(); clock1 = timer.cpu_time();
  msg << "Time for Spatial Particle operation = " << clock1 << endl;

  timer.clear(); timer.start();
  IA.A = IA.A - IA.B;
  timer.stop(); clock1 = timer.cpu_time();
  msg << "Time for Interact Particle operation = " << clock1 << endl;

  timer.clear(); timer.start();
  for (i=0; i < SA.getLocalNum(); ++i)
    SA.A[i] = SA.A[i] - SA.B[i];
  timer.stop(); clock1 = timer.cpu_time();
  msg << "Time for Spatial Particle loop = " << clock1 << endl;

  timer.clear(); timer.start();
  for (i=0; i < IA.getLocalNum(); ++i)
    IA.A[i] = IA.A[i] - IA.B[i];
  timer.stop(); clock1 = timer.cpu_time();
  msg << "Time for Interact Particle loop = " << clock1 << endl;

  return 0;
}
  

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

