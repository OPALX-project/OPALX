#include "FFTBottomSolver.h"

template <class Level>
FFTBottomSolver<Level>::FFTBottomSolver(Mesh_t *mesh,
                                 FieldLayout_t *fl,
                                 std::string greensFunction,
                                 std::string bcz)
    : FFTPoissonSolver(mesh, fl, greensFunction, bcz)
{
    for (lo_t i = 0; i < 6; ++i) {
        if (Ippl::getNodes()>1)
            bc_m[i] = new ParallelInterpolationFace<double, 3, Mesh_t, Center_t>(i);
        else
            bc_m[i] = new InterpolationFace<double, 3, Mesh_t, Center_t>(i);
    }
    
    rho_m.initialize(*mesh,
                     *fl,
                     GuardCellSizes<3>(1),
                     bc_m);
}


template <class Level>
FieldLayout_t* FFTBottomSolver<Level>::initMesh(Mesh_t *mesh,
                                                       AmrOpal* amrobject_p)
{
    if ( amrobject_p == nullptr )
        throw OpalException("FFTBottomSolver::initMesh()",
                            "No AMR object initialized.");
    
    // find the number of grid points in each dimension
    const AmrGeometry_t& geom = amrobject_p->Geom(0);
    const AmrBox_t& bx = geom.Domain();
    AmrIntVect_t iv = bx.size();
    
    NDIndex<3> domain;
    domain[0] = Index(0, iv[0]);
    domain[1] = Index(0, iv[1]);
    domain[2] = Index(0, iv[2]);
    
    return new Mesh_t(domain);
}


template <class Level>
FieldLayout_t* FFTBottomSolver<Level>::initFieldLayout(Mesh_t *mesh,
                                                       AmrOpal* amrobject_p)
{
    if ( amrobject_p == nullptr )
        throw OpalException("FFTBottomSolver::initFieldLayout()",
                            "No AMR object initialized.");
    
    const AmrGrid_t& ba = amrobject_p->boxArray(0);
    const AmrProcMap_t& dmap = amrobject_p->DistributionMap(0);
    auto pmap = dmap.ProcessorMap();
    
    std::vector< NDIndex<3> > regions;
    std::vector< int > nodes;
    for (uint i = 0; i < pmap.size(); ++i) {
        AmrBox_t bx = ba[i];
        
        NDIndex<3> range;
        for (int j = 0; j < 3; ++j)
            range[j] = Index(bx.smallEnd(j), bx.bigEnd(j));
        
        regions.push_back( range );
        nodes.push_back( pmap[i] );
    }
    
    return new FieldLayout_t(*mesh.get(),
                             &regions[0],
                             &regions[0] + regions.size(),
                             &nodes[0],
                             &nodes[0] + nodes.size());
}


template <class Level>
void FFTBottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                            const Teuchos::RCP<mv_t>& b )
{
    Vector_t hr;
    
    this->vector2field_m(b);
    
    FFTPoissonSolver::computePotential(rho_m, hr);
    
    this->field2vector_m(x);
}


template <class Level>
void FFTBottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A)
{
    // do nothing here
};


template <class Level>
void FFTBottomSolver<Level>::fillMap_m(const AmrGrid_t& ba,
                                const AmrProcMap_t& dmap)
{
    map_m.clear();
    
    lo_t localidx = 0;
    for (amrex::MFIter mfi(ba, dmap, true); mfi.isValid(); ++mfi) {
        const amrex::Box&       tbx = mfi.tilebox();
        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    map_m[localidx++] = iv;
//                     go_t gidx = iv[0] + (iv[1] + ny * iv[2]) * nx;
//                     map_m[gidx] = iv;
                        
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
}


template <class Level>
void FFTBottomSolver<Level>::field2vector_m(const Teuchos::RCP<mv_t>& vector)
{
    for (size_t i = 0; i < vector->getLocalLength(); ++i) {
        AmrIntVect_t iv = map_m[i];
        vector->replaceLocalValue(i, 0, rho_m[iv[0]][iv[1]][iv[2]].get());
    }
}
    

template <class Level>
void FFTBottomSolver<Level>::vector2field_m(const Teuchos::RCP<mv_t>& vector)
{
    for (size_t i = 0; i < vector->getLocalLength(); ++i) {
        AmrIntVect_t iv = map_m[i];
        rho_m[iv[0]][iv[1]][iv[2]] = vector->getData(0)[i];
    }
}
