/*!
 * @file testRestrictionMatrix.cpp
 * @author Matthias Frey
 * 
 */
#include <iostream>

#include "Ippl.h"

#include <vector>

#include "build.h"

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nlevs;
  bool verbose;
};


void myUpdate(Teuchos::RCP<Epetra_Vector>& y,
          Teuchos::RCP<Epetra_Vector>& x)
{
    int localnum = y->MyLength();
    for (int i = 0; i < localnum; ++i) {
        if ( (*x)[i] != 0 )
            (*y)[i] = (*x)[i];
    }
    
}


void test(TestParams& parms)
{
    int nlevs = parms.nlevs + 1;
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }
    
    IntVect domain_lo(D_DECL(0 , 0, 0)); 
    
    IntVect domain_hi(D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1)); 
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, coord, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    Vector<DistributionMapping> dmap(nlevs);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        BoxList bl;
        Box b1(IntVect(D_DECL(0, 4, 4)), IntVect(D_DECL(3, 7, 7)));
        
        bl.push_back(b1);
        
        Box b2(IntVect(D_DECL(4, 4, 4)), IntVect(D_DECL(11, 11, 11)));
        
        bl.push_back(b2);
        
        Box b3(IntVect(D_DECL(14, 6, 6)), IntVect(D_DECL(15, 9, 9)));
        
        bl.push_back(b3);
        
        
        ba[1].define(bl);//define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
        
        std::cout << ba[lev] << std::endl;
        
        dmap[lev].define(ba[lev]);
    }
    
    
    
    
    Epetra_MpiComm epetra_comm = Ippl::getComm();
    
    std::vector<Teuchos::RCP<Epetra_Map> > maps(nlevs);
    
    
    
    
    for (int i = 0; i < nlevs; ++i)
        buildMap(maps[i], ba[i], dmap[i],  geom[i], epetra_comm, i);
    
    
    
    
    Teuchos::RCP<Epetra_CrsMatrix> R = Teuchos::null;
    
    IntVect rv(D_DECL(2, 2, 2));
    
    buildRestrictionMatrix(R, maps, ba, dmap, geom, rv, epetra_comm, 1);
    
    
    // fine
    Teuchos::RCP<Epetra_Vector> x = Teuchos::null;
    buildVector(x, maps[1], 2.0);
    
    // coarse
    Teuchos::RCP<Epetra_Vector> y = Teuchos::null;
    buildVector(y, maps[0], 1.0);
    
    Teuchos::RCP<Epetra_Vector> z = Teuchos::null;
    buildVector(z, maps[0], 0.0);
    
    // 7 = R * x
    R->Multiply(false, *x, *z);
    
    // y += z
    myUpdate(y, z);
    
    std::cout << *y << std::endl;
    
    
    container_t rhs(1);
    container_t phi(1);
    container_t efield(1);
    
    if ( parms.verbose ) {
        for (int lev = 0; lev < 1; ++lev) {
            //                                                                       # component # ghost cells                                                                                                                                          
            rhs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 0));
            rhs[lev]->setVal(0.0);
            
            // not used (only for plotting)
            phi[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev],    1          , 1));
            efield[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba[lev], dmap[lev], AMREX_SPACEDIM, 1));
                
            phi[lev]->setVal(0.0, 1);
            efield[lev]->setVal(0.0, 1);
        }
    }
    
    trilinos2amrex(*rhs[0], y);
    
    if ( parms.verbose )
        writeYt(rhs, phi, efield, geom, rr, 1.0, "testRestrictionMatrix");
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nlevs", parms.nlevs);
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << "Num levels: ";
    std::cout << parms.nlevs << std::endl;
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }
  
  test(parms);
  
  amrex::Finalize();
}
