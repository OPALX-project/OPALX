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
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

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
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************

 This test program sets up a simple sine-wave electric field in 3D, 
  creates a population of particles with random q/m values (charge-to-mass
 ratio) and velocities, and then tracks their motions in the static
   electric field using nearest-grid-point interpolation. 
***************************************************************************/

#include "Ippl.h"
//#include "Particle/IntNGP.h"

// dimension of our positions
const unsigned Dim = 3;

// typedef our particle layout type
typedef ParticleSpatialLayout<double,Dim> playout_t;
typedef ParticleSpatialLayout<double,Dim>::SingleParticlePos_t Vector_t;
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;

const int nx=32, ny=32, nz=32;
const unsigned int totalP = 1000;
const int nt = 100;

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {
public:
  ParticleAttrib<double>     qm; // charge-to-mass ratio
  typename PL::ParticlePos_t V;  // particle velocity
  typename PL::ParticlePos_t E;  // electric field at particle position

  Vector_t rmin_m,rmax_m;
  Vector_t hr_m;
  Vektor<int,Dim> nr_m;

  ChargedParticles(PL* pl) : IpplParticleBase<PL>(pl) {
    // register the particle attributes
    addAttribute(qm);
    addAttribute(V);
    addAttribute(E);
  }

  inline const Mesh_t& getMesh() const { return getLayout().getLayout().getMesh(); }
  inline Mesh_t& getMesh() { return getLayout().getLayout().getMesh(); }
  inline const FieldLayout_t& getFieldLayout() const {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }
  inline FieldLayout_t& getFieldLayout() {
    return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
  }
  void boundp() {
    
      bounds(this->R, rmin_m, rmax_m);
    rmin_m = -rmin_m;
    
    for(int i=0; i<Dim; i++)  
      rmax_m[i] = max(rmax_m[i], rmin_m[i]);
    rmin_m = -rmax_m;
    
    Vector_t dr = (rmax_m - rmin_m)*(1.05);
    
    rmin_m -= dr;
    rmax_m += dr;
    hr_m = (rmax_m - rmin_m) / (32.0 - 1.0); 
    
    INFOMSG(" hr=" <<  hr_m );
    INFOMSG(" rmax=" << rmax_m);
    INFOMSG(" rmin=" << rmin_m << endl);
    // rescale mesh 
    //getMesh().set_meshSpacing( &(hr_m[0])  );
    //getMesh().set_origin( rmin_m ); 
  }
};

int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  msg << "Particle test PIC3d: Begin." << endl;

  // potential phi = phi0 * sin(2*pi*x/Lx) * cos(4*pi*y/Ly)
  double phi0 = 0.1*nx;             // electric potential amplitude

  // create interpolater object (nearest-grid-point method)
  IntNGP myinterp;

  // create layout objects
  Index I(nx), J(ny), K(nz);
  Index I1(nx+1), J1(ny+1), K1(nz+1);
  UniformCartesian<Dim> mymesh(I1,J1,K1);
  FieldLayout<Dim> FL(I,J,K);

  // initialize static electric field
  Field<Vektor<double,Dim>,Dim> EFD(mymesh,FL);
  
  assign(EFD[I][J][K](0), -2.0*pi*phi0/nx * cos(2.0*pi*(I+0.5)/nx) * 
                    cos(4.0*pi*(J+0.5)/ny) * cos(pi*(K+0.5)/nz));

  assign(EFD[I][J][K](1),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) *
                     sin(4.0*pi*(J+0.5)/ny));
  
  // create an empty ChargedParticles object, setting it to use periodic BC's
  playout_t* PL = new playout_t(FL);
  ChargedParticles<playout_t> P(PL);
  P.getBConds()[0] = ParticlePeriodicBCond;
  P.getBConds()[1] = ParticlePeriodicBCond;
  P.getBConds()[2] = ParticlePeriodicBCond;
  P.getBConds()[3] = ParticlePeriodicBCond;
  P.getBConds()[4] = ParticlePeriodicBCond;
  P.getBConds()[5] = ParticlePeriodicBCond;

  // initialize the particle object: do all initialization on one node,
  // and distribute to others
  P.create(totalP / Ippl::getNodes());

  // quiet start for particle positions
  assign(P.R(0), IpplRandom * nx);
  assign(P.R(1), IpplRandom * ny);
  assign(P.R(2), IpplRandom * nz);

  // random initialization for charge-to-mass ratio
  assign(P.qm, (2 * qmmax * IpplRandom) - qmmax);

  P.boundp();
  // redistribute particles based on spatial layout
  P.update(); 
  
  // begin main timestep loop
  
  msg << "Starting iterations ..." << endl;
  for (unsigned int it=0; it<nt; it++) {
    // advance the particle positions
    // basic leapfrogging timestep scheme.  velocities are offset
    // by half a timestep from the positions.
    
    assign(P.R, P.R + dt * P.V);
    
    // update particle distribution across processors
    P.update();

    // gather the local value of the E field
    P.E.gather(EFD, P.R, myinterp);

    // advance the particle velocities
    assign(P.V, P.V + dt * P.qm * P.E);
    Vektor<double,Dim> ek = sum(P.V);
    msg << "Finished iteration " << it << " ek= " << ek << endl;
  }

  Ippl::Comm->barrier();
  msg << "Particle test PIC3d: End." << endl;
  return 0;
}


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

