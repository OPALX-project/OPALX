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
 *
 ***************************************************************************/

/***************************************************************************

 This test program to check the getNeighborNode

 Usage:
 
 mpirun -np 2 neighNodeTest1 128 128 128 10000 --commlib mpi --info 0 | grep neighN | sort

neighNodeTest1{0}> dim=0 NeighburPID= 4
neighNodeTest1{0}> dim=1 NeighburPID= 2
neighNodeTest1{0}> dim=2 NeighburPID= 1
neighNodeTest1{1}> dim=0 NeighburPID= 5
neighNodeTest1{1}> dim=1 NeighburPID= 3
neighNodeTest1{1}> dim=2 NeighburPID= 0
neighNodeTest1{2}> dim=0 NeighburPID= 6
neighNodeTest1{2}> dim=1 NeighburPID= 0
neighNodeTest1{2}> dim=2 NeighburPID= 3
neighNodeTest1{3}> dim=0 NeighburPID= 7
neighNodeTest1{3}> dim=1 NeighburPID= 1
neighNodeTest1{3}> dim=2 NeighburPID= 2
neighNodeTest1{4}> dim=0 NeighburPID= 0
neighNodeTest1{4}> dim=1 NeighburPID= 6
neighNodeTest1{4}> dim=2 NeighburPID= 5
neighNodeTest1{5}> dim=0 NeighburPID= 1
neighNodeTest1{5}> dim=1 NeighburPID= 7
neighNodeTest1{5}> dim=2 NeighburPID= 4
neighNodeTest1{6}> dim=0 NeighburPID= 2
neighNodeTest1{6}> dim=1 NeighburPID= 4
neighNodeTest1{6}> dim=2 NeighburPID= 7
neighNodeTest1{7}> dim=0 NeighburPID= 3
neighNodeTest1{7}> dim=1 NeighburPID= 5
neighNodeTest1{7}> dim=2 NeighburPID= 6

***************************************************************************/

#include "Ippl.h"
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

template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {
public:

    ChargedParticles(PL* pl, e_dim_tag decomp[Dim]) :
        IpplParticleBase<PL>(pl)
    {
	setBCAllOpen();
        for(int i=0; i<Dim; i++)
            decomp_m[i]=decomp[i];
    }

    inline const Mesh_t& getMesh() const { return this->getLayout().getLayout().getMesh(); }

    inline Mesh_t& getMesh() { return this->getLayout().getLayout().getMesh(); }

    inline const FieldLayout_t& getFieldLayout() const {
        return dynamic_cast<FieldLayout_t&>( this->getLayout().getLayout().getFieldLayout());
    }

    inline FieldLayout_t& getFieldLayout() {
        return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
    }

    void myUpdate() {

        double hz   = hr_m[2];
        double zmin = rmin_m[2];
        double zmax = rmax_m[2];

	bounds(this->R, rmin_m, rmax_m);

	NDIndex<Dim> domain = this->getFieldLayout().getDomain();
	
	for(int i=0; i<Dim; i++)
	    nr_m[i] = domain[i].length();

	for(int i=0; i<Dim; i++)
	    hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1.0);

	getMesh().set_meshSpacing(&(hr_m[0]));
	getMesh().set_origin(rmin_m);

        this->update();
    }

private:

    inline void setBCAllOpen() {
        for (int i=0; i < 2*Dim; i++) {
            this->getBConds()[i] = ParticleNoBCond;
        }
    }

    Vektor<int,Dim> nr_m;

    Vector_t hr_m;
    Vector_t rmin_m;
    Vector_t rmax_m;
    e_dim_tag decomp_m[Dim];

};

int main(int argc, char *argv[]){
    Ippl ippl(argc, argv);
    Inform msg(argv[0]);
    Inform msg2all(argv[0],INFORM_ALL_NODES);

    Vektor<int,Dim> nr(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

    const unsigned int totalP = atoi(argv[4]);

    e_dim_tag decomp[Dim];

    NDIndex<Dim> domain;
    
    for(int d=0; d<Dim; d++) {
	domain[d] = domain[d] = Index(nr[d] + 1);
        decomp[d] = PARALLEL;
    }

    // create mesh and layout objects for this problem domain
    Mesh_t*        mesh            = new Mesh_t(domain);
    FieldLayout_t* FL              = new FieldLayout_t(*mesh, decomp);
    playout_t*     PL              = new playout_t(*FL, *mesh);
    ChargedParticles<playout_t>* P = new ChargedParticles<playout_t>(PL,decomp);
    
    P->create(totalP);

    for (unsigned long int i = 0; i< totalP; i++) {
        for (int d = 0; d<3; d++)
            P->R[i](d) =  IpplRandom() * nr[d];
    }

    P->myUpdate();
    
    msg << P->getFieldLayout() << endl;

    int npid = 0;
    for (unsigned d=0; d<Dim; d++) {
	for (unsigned n=0; n<Ippl::getNodes(); n++) {
	    if ((npid=PL->getNeighborNode(d,n)) != -1) 
		msg2all << "dim=" << d << " NeighburPID= " << npid << endl;
	}
    }

    msg << argv[0] << " test End." << endl;
    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/

