#include "MueLuPreconditioner.h"

#include <MueLu_CreateTpetraPreconditioner.hpp>



MueLuPreconditioner::MueLuPreconditioner(const AmrIntVect_t& grid,
                                         const bool& rebalance)
    : prec_mp(Teuchos::null),
      coords_mp(Teuchos::null),
      grid_m(grid),
      rebalance_m(rebalance)
{
    this->init_m();
}


void MueLuPreconditioner::create(const Teuchos::RCP<amr::matrix_t>& A) {

    coords_mp = Teuchos::null;

    if ( rebalance_m ) {
        coords_mp = Teuchos::rcp( new amr::multivector_t(A->getDomainMap(),
                                                         AMREX_SPACEDIM, false) );
    
        const Teuchos::RCP<const amr::dmap_t>& map_r = A->getMap();
    
        for (std::size_t i = 0; i < map_r->getNodeNumElements(); ++i) {
            AmrIntVect_t iv = deserialize_m(map_r->getGlobalElement(i));
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                coords_mp->replaceLocalValue(i, 0, iv[d]);
        }
    }

    prec_mp = MueLu::CreateTpetraPreconditioner(A, params_m, coords_mp);

    coords_mp = Teuchos::null;
}


Teuchos::RCP<amr::operator_t> MueLuPreconditioner::get() {
    return prec_mp;
}


void MueLuPreconditioner::init_m() {
    params_m.set("problem: type", "Poisson-3D");
//    params_m.set("problem: symmetric", false);
    params_m.set("verbosity", "extreme");
    params_m.set("number of equations", 1);
    params_m.set("max levels", 8);
    params_m.set("cycle type", "V");

    params_m.set("coarse: max size", grid_m[0]);
    params_m.set("multigrid algorithm", "sa");
    
    params_m.set("repartition: enable", rebalance_m);
    params_m.set("repartition: rebalance P and R", rebalance_m);
    params_m.set("repartition: partitioner", "zoltan2");
    params_m.set("repartition: min rows per proc", grid_m[0]);
    params_m.set("repartition: start level", 1);
    
    params_m.set("smoother: type", "CHEBYSHEV");
    params_m.set("smoother: pre or post", "both");

    params_m.set("coarse: type", "KLU2");

    params_m.set("aggregation: type", "uncoupled");

    params_m.set("transpose: use implicit", true);

    params_m.set("reuse: type", "none");
}


MueLuPreconditioner::AmrIntVect_t
MueLuPreconditioner::deserialize_m(const std::size_t& gidx)
{
    int nxny = grid_m[0] * grid_m[1];
    int z = gidx / nxny;
    int y = (gidx %  nxny) / grid_m[0];
    int x = gidx - z * nxny - y * grid_m[0];
    return AmrIntVect_t(x, y, z);
}
