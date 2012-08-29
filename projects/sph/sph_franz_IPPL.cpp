#include "Ippl.h"
#include "mpi.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <vector>

// dimension of our positions
const unsigned Dim = 3;

// some typedefs
//Particle Layout Interaction 
typedef ParticleInteractLayout<double,Dim> playout_inter;

//Fieldlayout
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>       Field_t;

#define GUARDCELL 1

template<class T, unsigned Dim>
class SPHparticles : public ParticleBase< ParticleInteractLayout<T, Dim> > {

public:

  double nmax;

    // the attributes for this set of particles.
  ParticleInteractAttrib<T> A;
  ParticleInteractAttrib<T> B;

 // ParticleInteractAttrib<T> re;
  typename playout_inter::ParticlePos_t acc;				//particle acceleration	
  typename playout_inter::ParticlePos_t v;				//particle velocity
  typename playout_inter::ParticlePos_t vold;				//particle velocity at old timestep
  typename playout_inter::ParticlePos_t Rold;				//particle position at old timestep	


  ParticleAttrib<T> n;							//particle density
  ParticleAttrib<T> m;							//particle mass
  ParticleInteractAttrib<Vektor<double,Dim*Dim> > dvdx;			//particle velocity divergence

  SPHparticles(FieldLayout_t& fl)
    : ParticleBase< ParticleInteractLayout<T, Dim, Mesh_t> >(
			    new ParticleInteractLayout<T,Dim, Mesh_t>(fl)) {

      // register our attributes with the base class
      addAttribute(A);
      addAttribute(B);
 //     addAttribute(re);
      addAttribute(acc);
      addAttribute(v);
      addAttribute(n);
      addAttribute(m);
      addAttribute(dvdx);
      addAttribute(vold);
      addAttribute(Rold);
  }
  inline Vektor<T,Dim> getR(unsigned i) { return this->R[i]; }
};


int main(int argc, char *argv[]){
    // init IPPL
    Ippl ippl(argc, argv);
    Inform message(argv[0]);
    Inform message2all(argv[0],INFORM_ALL_NODES);
    Inform messageNode1(argv[0],1);

    Vektor<int,Dim> nr(2,2,2);

    unsigned int np = 8;

    bool gCells = 1;

    if (gCells)
        message << "Using guard cells" << endl;
    else
        message << "Not using guard cells" << endl;

    e_dim_tag decomp[Dim];
    int serialDim = 5;
    message << "Serial dimension is " << serialDim  << endl;

    Mesh_t *mesh;
    FieldLayout_t *FL;

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

    decomp[0] = PARALLEL;
    decomp[1] = PARALLEL;
    decomp[2] = PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);

    SPHparticles<double,Dim> mp(*FL);

    double dx = 1.0;//sqrt(2.0);//1.0;
    Vektor<double,Dim> pos;
    int currlocalpart = mp.getLocalNum();
    message<< " Number of particle on processor " << currlocalpart<<endl;
    int currtotalpart = mp.getTotalNum();
    message<< " Number of particle ALL processor " << currtotalpart<<endl;
    if(mp.singleInitNode()){
      mp.create(np);
      for(int i=0; i<2; i++){
	for(int j=0; j<2; j++){
	  for(int k=0; k<2; k++){
	    mp.R[currlocalpart++] = Vektor<double,3>(i*dx, j*dx, k*dx);	    
	  }
	}
      }
    }
   
    // Set interaction radius for NN-Search with constant re
    mp.getLayout().setInteractionRadius(3.1 * 1.0/sqrt(3.0)*dx);
    //mp.getLayout().setInteractionRadius(dx);
    // Set interaction radius for NN-Search with constant re
    //mp.re = dx;
    //mp.getLayout().setInteractionRadius(mp.re);

    mp.update();

    //Check particle distribution on nodes
    currlocalpart = mp.getLocalNum();
    message2all<< " Get local num " << currlocalpart<<endl;
    currtotalpart = mp.getTotalNum();
    message<< " Number of particle ALL processor " << currtotalpart<<endl;  


    //Find Nearest Neigbours   
    double re;
    re = mp.getLayout().getInteractionRadius(re);
    message<< " Interaction Radius " << re << endl;
    int j;
    double rr2;   
    SPHparticles<double,Dim>::pair_iterator jlist, jlist_end;

    // Finding neighbors for the first time...
    double clock1;
    Timer timer;
    timer.clear(); timer.start();
    for(unsigned int i=0; i<mp.getLocalNum(); i++){
      double testtesttes =mp.getLayout().getMaxInteractionRadius();
      mp.getLayout().getPairlist(i,jlist,jlist_end,mp);
      message2all << " Find for 1st time: Particle ID "<< mp.ID[i]<< ",  Num of Neighbors " << jlist_end -jlist <<endl; 
    }
    timer.stop(); clock1 = timer.cpu_time();
    message << "Time for Nbssearch 1 = " << clock1 << endl;

    // Finding neighbors for the SECOND time...
    timer.clear(); timer.start();
    for(unsigned int i=0; i<mp.getLocalNum(); i++){
      mp.getLayout().getPairlist(i,jlist,jlist_end,mp);   
    }
    timer.stop(); clock1 = timer.cpu_time();
    message << "Time for Nbssearch  2 = " << clock1 << endl;

    for(unsigned int i=0; i<mp.getLocalNum(); i++){
      mp.m[i] = 1.0;
    }
  
    double m, n, rr;
    for(unsigned int i=0; i<mp.getLocalNum(); i++){
      m = mp.m[i];
      n = m;
      std::cout<<" ---------------------------------i " << i << " m " << m << " n " << n << endl;
      mp.getLayout().getPairlist(i,jlist,jlist_end,mp);
      for ( ; jlist != jlist_end ; jlist++ ) {
	SPHparticles<double,Dim>::pair_t jdata = *jlist;
	j = jdata.first;
	rr = jdata.second;
	rr = sqrt(rr);	
	m = mp.m[j];
	n += m ;
	std::cout<<" j " << j << " mj " << m << " nj " << n << " rr " << rr << endl; 
      }
      mp.n[i]=n;
    }
     
    Ippl::Comm->barrier();
    message << "my_first_IPPL has ENDED." << endl;
    return 0;
}


