#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Ippl.h"

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <vector>
#include <map>
#include <chrono>

typedef UniformCartesian<3, double> Mesh_t;

typedef Cell Center_t;

typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
// typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::RCP;
using Teuchos::rcp;

typedef Tpetra::Vector<>::scalar_type scalar_type;
typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
typedef Tpetra::Map<> map_type;
typedef Tpetra::Vector<> vector_type;


std::unique_ptr<Mesh_t> mesh;
std::unique_ptr<FieldLayout_t> fl;

std::map<int, amrex::IntVect> mapper;

#define nx 32
#define ny 32
#define nz 32

void setupField() {
    
    using amrex::RealBox;
    using amrex::IntVect;
    using amrex::Geometry;
    using amrex::BoxArray;
    using amrex::Box;
    using amrex::DistributionMapping;
    
    NDIndex<3> domain;
    domain[0] = Index(0, nx);
    domain[1] = Index(0, ny);
    domain[2] = Index(0, nz);
    
    // create prototype mesh and layout objects for this problem domain
    mesh.reset(new Mesh_t(domain));
    
    
    // ================
    int max_grid_size = nx / 4;
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 0.5);
    }

    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(nx - 1, ny - 1, nz - 1); 
    const Box amr_domain(domain_lo, domain_hi);

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 0; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Geometry geom;
    geom.define(amr_domain, &real_box, coord, is_per);

    BoxArray ba;
    ba.define(amr_domain);
    ba.maxSize(max_grid_size);
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    DistributionMapping dmap;
    
    dmap.define(ba);
    
    // ================
    
    
    
    const auto& pmap = dmap.ProcessorMap();
    
    std::vector< NDIndex<3> > regions;
    std::vector< int > nodes;
    for (uint i = 0; i < pmap.size(); ++i) {
        Box bx = ba[i];
        
        NDIndex<3> tmp;
        for (int j = 0; j < 3; ++j)
            tmp[j] = Index(bx.smallEnd(j), bx.bigEnd(j));
        
        regions.push_back( tmp );
        nodes.push_back( pmap[i] );
    }
    
    
    
    
//     for (uint i = 0; i < regions.size(); ++i)
//         std::cout << regions[i] << " " << nodes[i] << std::endl;

    fl.reset(new FieldLayout_t(*mesh, &regions[0], &regions[0] + regions.size(), &nodes[0], &nodes[0] + nodes.size()));
}


void buildMap(const amrex::BoxArray& ba,
              const amrex::DistributionMapping& dmap)
{
//     int localid = 0;
    for (amrex::MFIter mfi(ba, dmap, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box&       tbx = mfi.tilebox();
        
        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    amrex::IntVect iv(D_DECL(i, j, k));
                    int gidx = iv[0] + (iv[1] + ny * iv[2]) * nx;
                    mapper[gidx/*localid++*/] = iv;
                        
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
}

void test(const Teuchos::RCP<Teuchos::MpiComm<int> >& comm)
{
    if ( (nx * ny * nz) % comm->getSize() ) {
        std::cout << "Error: " << nx * ny * nz << " % " << comm->getSize() << " != 0" << std::endl;
        return;
    }
    
    const size_t numLocalEntries = nx * ny * nz / comm->getSize();
    const Tpetra::global_size_t numGlobalEntries = nx * ny * nz;
        
    const global_ordinal_type indexBase = 0;
    
    RCP<const map_type> contigMap =
        rcp (new map_type (numGlobalEntries, numLocalEntries, indexBase, comm));
    
    vector_type x (contigMap);
    
    x.putScalar(1.0);
    
    setupField();
    
    Field_t rho;
    BConds<double, 3, Mesh_t, Center_t> bc_m;
    
    for (int i = 0; i < 2 * 3; ++i) {

        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, 3, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, 3, Mesh_t, Center_t>(i);
        }
    }
    
    rho.initialize(*mesh,
                   *fl,
                   GuardCellSizes<3>(1),
                   bc_m);
  
    
    auto start = std::chrono::high_resolution_clock::now();
    
    double globalsum = 0;
    for (size_t i = 0; i < x.getLocalLength(); ++i) {
            std::size_t gidx = contigMap->getGlobalElement(i);
            rho[mapper[gidx][0]][mapper[gidx][1]][mapper[gidx][2]] = x.getData()[i];
            
            globalsum += rho[mapper[gidx][0]][mapper[gidx][1]][mapper[gidx][2]].get();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> diff = end-start;
    
    std::cout << "Time: " << diff.count() << std::endl;
    
    double pete_sum = sum(rho);
    
    allreduce(globalsum, 1, std::plus<double>());
    
    if ( Ippl::myNode() == 0) {
        std::cout << "reference: " << nx * ny * nz << std::endl
                  << "computed:  " << globalsum << std::endl
                  << "pete sum:  " << pete_sum << std::endl;
    }
}


int main (int argc, char *argv[])
{
    Ippl ippl(argc, argv);
    
    amrex::Initialize(argc, argv, false, Ippl::getComm());
    
    using Teuchos::RCP;
    
    RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp( new Teuchos::MpiComm<int>( Teuchos::opaqueWrapper(Ippl::getComm()) ) );
    
    test(comm);
    
    return 0;
}
