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

This test program sets up a simple sine-wave electric field in 3D,
  creates a population of particles with random q/m values (charge-to-mass
  ratio) and velocities, and then tracks their motions in the static
  electric field using nearest-grid-point interpolation.

Usage:

 mpirun -np 2 testAmrPartBunch 128 128 128 10000 10

***************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

#include "PartBunch.h"

const double dt = 1.0;          // size of timestep


int main(int argc, char *argv[]){
    Ippl ippl(argc, argv);
    Inform msg(argv[0]);
    Inform msg2all(argv[0],INFORM_ALL_NODES);

    Vektor<int,Dim> nr(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

    const unsigned int totalP = atoi(argv[4]);
    const unsigned int nt     = atoi(argv[5]);

    msg << "Particle test testAmrPartBunch " << endl;
    msg << "nt " << nt << " Np= " << totalP << " grid = " << nr <<endl;

    e_dim_tag decomp[Dim];
    unsigned serialDim = 2;

    msg << "Serial dimension is " << serialDim  << endl;

    Mesh_t *mesh;
    FieldLayout_t *FL;
//     PartBunch<playout_t>  *bunch;

    NDIndex<Dim> domain;
    for (unsigned i=0; i<Dim; i++)
        domain[i] = domain[i] = Index(nr[i] + 1);

    for (unsigned d=0; d < Dim; ++d)
        decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);

    /*
     * In case of periodic BC's define
     * the domain with hr and rmin
     */
    Vector_t hr(1.0);
    Vector_t rmin(0.0);
    Vector_t rmax(nr);
    
    PartBunchBase* bunch = new PartBunch<playout_t>(PL,hr,rmin,rmax,decomp);

    // initialize the particle object: do all initialization on one node,
    // and distribute to others

    unsigned long int nloc = totalP / Ippl::getNodes();

    bunch->create(nloc);
    for (unsigned long int i = 0; i< nloc; i++) {
        for (int d = 0; d<3; d++)
            bunch->getR(i)(d) =  IpplRandom() * nr[d];
    }

    double q = 1.0/totalP;

    // random initialization for charge-to-mass ratio
    assign(bunch->getQM(),q);

    msg << "particles created and initial conditions assigned " << endl;

    // redistribute particles based on spatial layout
    bunch->myUpdate();

    msg << "initial update and initial mesh done .... Q= " << sum(bunch->getQM()) << endl;
    msg << bunch->getMesh() << endl;
    msg << bunch->getFieldLayout() << endl;

    msg << "scatter test done delta= " <<  bunch->scatter() << endl;

    bunch->initFields();
    msg << "bunch->initField() done " << endl;

    // begin main timestep loop
    msg << "Starting iterations ..." << endl;
    for (unsigned int it=0; it<nt; it++) {
        bunch->gatherStatistics();
        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.
        assign(bunch->getR(), bunch->getR() + dt * bunch->getP());

        // update particle distribution across processors
        bunch->myUpdate();

        // gather the local value of the E field
        bunch->gatherCIC();

        // advance the particle velocities
        assign(bunch->getP(), bunch->getP() + dt * bunch->getQM() * bunch->getE());
        msg << "Finished iteration " << it << " - min/max r and h " << bunch->getRMin()
            << bunch->getRMax() << bunch->getHr() << endl;
    }
    Ippl::Comm->barrier();
    msg << "Particle test testAmrPartBunch: End." << endl;
    return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
