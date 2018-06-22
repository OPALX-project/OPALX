/*!
 * @file testInitialBox.cpp
 * @author Matthias Frey
 * 
 * @details Initialize particles first in a default box
 * and then initialize AMR. After that a tagging is performed.
 */

#include <iostream>

#include "Ippl.h"

#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "AmrParticleBase1.h"
#include "ParticleAmrLayout.h"
#include "PartBunchAmr.h"

#include "../AmrOpal.h"

#define Dim 3

using namespace amrex;

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase1<amrplayout_t> amrbase_t;
typedef PartBunchAmr<amrplayout_t> amrbunch_t;


void print(amrbunch_t* bunch) {
    std::cout << "Start print" << std::endl;
    // single core
    unsigned int lev = 0, nLocParticles = 0;
    for (unsigned int ip = 0; ip < bunch->getLocalNum(); ++ip) {
        while ( lev != bunch->m_lev[ip] ) {
            std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
            nLocParticles = 0;
            lev++;
        }
        nLocParticles++;
    }
    
    std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
}

void createRandomParticles(amrbunch_t *bunch, int N, int myNode, int seed = 1) {

  srand(seed);
  for (int i = 0; i < N; ++i) {
    bunch->createWithID(myNode * N + i + 1);
    bunch->qm[i] = 1.0; //(double)rand() / RAND_MAX;

    bunch->R[i][0] = 0.015625 * 0.5; // (double)rand() / RAND_MAX;
    bunch->R[i][1] = 0.015625 * 0.5; //(double)rand() / RAND_MAX;
    bunch->R[i][2] = 0.015625 * 0.5; //(double)rand() / RAND_MAX;  
  }
    
}

void initBunch(amrbunch_t* &bunch, int N) {
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n,  0.0);
        real_box.setHi(n,  1.0);
    }
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(31, 31, 31);
    
    const Box bx(domain_lo, domain_hi);
    
    int coord = 0;
    
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 0;
    
    
    Geometry geom;
    geom.define(bx, &real_box, 0, is_per);
    
    BoxArray ba(bx);
    ba.maxSize(16);
    
    DistributionMapping dmap;
    dmap.define(ba, ParallelDescriptor::NProcs());
    
    amrplayout_t* PL = new amrplayout_t(geom, dmap, ba);
    
    std::cout << "initBunch" << std::endl;
    std::cout << geom.CellSize(0) << " "
              << geom.CellSize(1) << " "
              << geom.CellSize(2) << std::endl;
    
    bunch = new PartBunchAmr<amrplayout_t>();
    bunch->initialize(PL);
    bunch->initializeAmr();
    
    createRandomParticles(bunch, N, Ippl::myNode(), 42);
    
    bunch->update();
    
    bunch->gatherStatistics();
    
    print(bunch);
}

void initAmr(AmrOpal* &myAmrOpal) {
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    
    int maxLevel = 3;
    Vector<int> nCells(3);
    nCells[0] = 64;
    nCells[1] = 64;
    nCells[2] = 64;
    
//     Geometry::setProbDomain(real_box);
    
    myAmrOpal = new AmrOpal(&real_box, maxLevel, nCells, 0 /* cartesian */);
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(nCells[0] - 1, nCells[1] - 1, nCells[2] - 1);
    
    const Box bx(domain_lo, domain_hi);
    BoxArray ba(bx);
    ba.maxSize(16);
    
    myAmrOpal->SetBoxArray(0, ba);
    
    DistributionMapping dmap;
    dmap.define(ba, ParallelDescriptor::NProcs());
    
    myAmrOpal->SetDistributionMap(0, dmap);
    
    std::cout << "initAmr" << std::endl;
    const Geometry& geom = myAmrOpal->Geom(0);
    std::cout << geom.CellSize(0) << " "
              << geom.CellSize(1) << " "
              << geom.CellSize(2) << std::endl;
}

void updateBunch(amrbunch_t* bunch, AmrOpal* myAmrOpal) {
    
    const Vector<Geometry>& geom = myAmrOpal->Geom();
    const Vector<DistributionMapping>& dmap = myAmrOpal->DistributionMap();
    const Vector<BoxArray>& ba = myAmrOpal->boxArray();
    const Vector<IntVect>& ref_rato = myAmrOpal->refRatio ();
    
    std::cout << ba[0] << std::endl;
    
    
    Vector<int> rr( ref_rato.size() );
    for (unsigned int i = 0; i < rr.size(); ++i) {
        rr[i] = ref_rato[i][0];
    }
    
    amrplayout_t* PLayout = &bunch->getLayout();
    PLayout->Define(geom, dmap, ba, rr);
    
    std::cout << geom[0].CellSize(0) << " "
              << geom[0].CellSize(1) << " "
              << geom[0].CellSize(2) << std::endl;
    
    bunch->update();
    
    bunch->gatherStatistics();
    
    print(bunch);
    
    myAmrOpal->setBunch(bunch);
}


int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    
    amrex::Initialize(argc,argv, false);
    
    PartBunchAmr<amrplayout_t>* bunch = nullptr;
    
    initBunch(bunch, 2);
    
    AmrOpal* myAmrOpal = nullptr;
    
    initAmr(myAmrOpal);
    
    updateBunch(bunch, myAmrOpal);
    
    myAmrOpal->setCharge( 0.125 );
    
    for (int i = 0; i <= myAmrOpal->finestLevel() && i < myAmrOpal->maxLevel(); ++i)
        myAmrOpal->regrid(i /*lbase*/, 0.0 /*time*/);
    
    const Vector<Geometry>& geom = myAmrOpal->Geom();
    for (unsigned int i = 0; i < geom.size(); ++i) {
        std::cout << geom[i].CellSize(0) << " "
                  << geom[i].CellSize(1) << " "
                  << geom[i].CellSize(2) << std::endl;
    }
    
    std::cout << "Finest level: " << myAmrOpal->finestLevel() << std::endl;
    
    bunch->gatherStatistics();
    
    bunch->python_format(0);
    
    print(bunch);
    
    delete myAmrOpal;
    delete bunch;
    
    return 0;
}