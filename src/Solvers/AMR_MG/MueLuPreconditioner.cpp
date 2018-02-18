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
                                                         1, false) );
    
        const Teuchos::RCP<const amr::dmap_t>& map_r = A->getMap();
    
        for (std::size_t i = 0; i < map_r->getNodeNumElements(); ++i)
            coords_mp->replaceLocalValue(i, 0, map_r->getGlobalElement(i));
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
    params_m.set("max levels", 4);
    params_m.set("cycle type", "V");

    params_m.set("coarse: max size", 2000);
    params_m.set("multigrid algorithm", "sa");
    
    params_m.set("repartition: enable", rebalance_m);
    params_m.set("repartition: rebalance P and R", rebalance_m);
    params_m.set("repartition: partitioner", "zoltan2");
    params_m.set("repartition: min rows per proc", 800);
    params_m.set("repartition: start level", 2);

    Teuchos::ParameterList reparms;
    reparms.set("algorithm", "rcb");
    //    reparms.set("partitioning_approach", "partition");

    params_m.set("repartition: params", reparms);
    
    params_m.set("smoother: type", "CHEBYSHEV");
    params_m.set("smoother: pre or post", "both");

    params_m.set("coarse: type", "KLU2");

    params_m.set("aggregation: type", "uncoupled");
    params_m.set("aggregation: min agg size", 3);
    params_m.set("aggregation: max agg size", 27);

    params_m.set("transpose: use implicit", true);

    params_m.set("reuse: type", "none");
}
