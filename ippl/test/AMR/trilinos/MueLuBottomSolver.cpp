#include "MueLuBottomSolver.h"

/* https://trilinos.org/wordpress/wp-content/uploads/2015/03/MueLu_tutorial.pdf
 * http://prod.sandia.gov/techlib/access-control.cgi/2014/1418624r.pdf
 */

MueLuBottomSolver::MueLuBottomSolver()
    : hierarchy_mp(Teuchos::null),
      finest_mp(Teuchos::null),
      A_mp(Teuchos::null),
      tolerance_m(1.0e-4)
{
    this->initMueLuList_m();
}


void MueLuBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                              const Teuchos::RCP<mv_t>& b)
{
    // MueLu requires Xpetra multivectors (wrap them)
    Teuchos::RCP<xmv_t> xx = MueLu::TpetraMultiVector_To_XpetraMultiVector(x);
    Teuchos::RCP<xmv_t> xb = MueLu::TpetraMultiVector_To_XpetraMultiVector(b);

    // InitialGuessIsZero = true
    // startLevel = 0
    hierarchy_mp->Iterate(*xb, *xx, 4/*tolerance_m*/, true, 0);
    
    // put multivector back
    x->assign(*util_t::MV2NonConstTpetraMV2(*xx));
}


void MueLuBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    
    A_mp = MueLu::TpetraCrs_To_XpetraMatrix<scalar_t, lo_t, go_t, node_t>(A);
    A_mp->SetFixedBlockSize(1); // only 1 DOF per node (pure Laplace problem)

    Teuchos::RCP<mv_t> coords_mp = Teuchos::rcp( new amr::multivector_t(A->getDomainMap(),
                                                                        1/*AMREX_SPACEDIM*/, false) );

    const Teuchos::RCP<const amr::dmap_t>& map_r = A->getMap();

    for (std::size_t i = 0; i < map_r->getNodeNumElements(); ++i) {
//            AmrIntVect_t iv = deserialize_m(map_r->getGlobalElement(i));
//            for (int d = 0; d < AMREX_SPACEDIM; ++d) 
        coords_mp->replaceLocalValue(i, 0, map_r->getGlobalElement(i)); //iv[d]);                                 
    }


    Teuchos::RCP<xmv_t> coordinates = MueLu::TpetraMultiVector_To_XpetraMultiVector(coords_mp);


    Teuchos::RCP<manager_t> mueluFactory = Teuchos::rcp(
        new pListInterpreter_t(mueluList_m)
	);

    // empty multigrid hierarchy with a finest level only
    hierarchy_mp = mueluFactory->CreateHierarchy();

    hierarchy_mp->GetLevel(0)->Set("A", A_mp);

    Teuchos::RCP<mv_t> nullspace = Teuchos::rcp(new mv_t(A->getRowMap(), 1));
    Teuchos::RCP<xmv_t> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector(nullspace);
    xnullspace->putScalar(1.0);

    hierarchy_mp->GetLevel(0)->Set("Nullspace", xnullspace);
    hierarchy_mp->GetLevel(0)->Set("Coordinates", coordinates);
    hierarchy_mp->IsPreconditioner(false);
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
   
    finest_mp = hierarchy_mp->GetLevel();
    finest_mp->Set("A", A_mp);
    finest_mp->Set("Coordinates", coordinates);
    finest_mp->Set("Nullspace", xnullspace);
    
    mueluFactory->SetupHierarchy(*hierarchy_mp);
}


std::size_t MueLuBottomSolver::getNumIters() {
    return 1;
}


void MueLuBottomSolver::initMueLuList_m() {
    mueluList_m.set("problem: type", "Poisson-3D");
//    mueluList_m.set("problem: symmetric", false);
    mueluList_m.set("verbosity", "extreme");
    mueluList_m.set("number of equations", 1);
    mueluList_m.set("max levels", 8);
    mueluList_m.set("cycle type", "V");

    mueluList_m.set("coarse: max size", 200);
    mueluList_m.set("multigrid algorithm", "sa");
    
    mueluList_m.set("repartition: enable", true);
    mueluList_m.set("repartition: rebalance P and R", true);
    mueluList_m.set("repartition: partitioner", "zoltan2");
    mueluList_m.set("repartition: min rows per proc", 800);
    mueluList_m.set("repartition: start level", 2);

    Teuchos::ParameterList reparms;
    reparms.set("algorithm", "rcb");
    //    reparms.set("partitioning_approach", "partition");

    mueluList_m.set("repartition: params", reparms);
    
    mueluList_m.set("smoother: type", "CHEBYSHEV");
    mueluList_m.set("smoother: pre or post", "both");

    mueluList_m.set("coarse: type", "KLU2");

    mueluList_m.set("aggregation: type", "uncoupled");
    mueluList_m.set("aggregation: min agg size", 3);
    mueluList_m.set("aggregation: max agg size", 27);

    mueluList_m.set("transpose: use implicit", true);

    mueluList_m.set("reuse: type", "none");
}


// void MueLuBottomSolver::setupSmoother_m(fManager_t& manager) {
//     std::string ifpackType = "RELAXATION";
//     Teuchos::ParameterList ifpackList;
//     ifpackList.set("relaxation: sweeps", (lo_t) 4);
//     ifpackList.set("relaxation: damping factor", (scalar_t) 1.0);
//     Teuchos::RCP<smootherPrototype_t> smootherPrototype =
//         Teuchos::rcp(new smoother_t(ifpackType, ifpackList));
//     
//     manager.SetFactory("Smoother", Teuchos::rcp(new smootherFactory_t(smootherPrototype)));
// }
// 
// 
// void MueLuBottomSolver::setupCoarseSolver_m(fManager_t& manager) {
//     Teuchos::RCP<smootherFactory_t> coarseSolverPrototype =
//         Teuchos::rcp( new solver_t() );
//     
//     Teuchos::RCP<smootherFactory_t> coarseSolverFact =
//         Teuchos::rcp( new smootherFactory_t(coarseSolverPrototype, Teuchos::null) );
//     
//     manager.SetFactory("CoarseSolver", coarseSolverFact);
// }
