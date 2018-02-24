#include "MueLuBottomSolver.h"

MueLuBottomSolver::MueLuBottomSolver()
    : hierarchy_mp(Teuchos::null),
      finest_mp(Teuchos::null),
      A_mp(Teuchos::null),
      tolerance_m(1.0e-4)
{
/*
    hierarchy_mp = Teuchos::rcp(new hierarchy_t());
    
    // WCYCLE also supported
    hierarchy_mp->SetCycle(MueLu::CycleType::VCYCLE);
    hierarchy_mp->SetMaxCoarseSize(200);
    hierarchy_mp->SetPRrebalance(false);
    hierarchy_mp->SetImplicitTranspose(true);
    hierarchy_mp->IsPreconditioner(false);
    
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    
    finest_mp = hierarchy_mp->GetLevel();
    
    finest_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
*/
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
    Teuchos::RCP<tentPFactory_t> TentativePFact = Teuchos::rcp( new tentPFactory_t() );
    Teuchos::RCP<saPFactory_t> SaPFact = Teuchos::rcp( new saPFactory_t() );
    Teuchos::RCP<transPFactory_t> RFact = Teuchos::rcp( new transPFactory_t());
    Teuchos::RCP<rAPFactory_t> AFact = rcp(new rAPFactory_t());
    AFact->setVerbLevel(Teuchos::VERB_HIGH);
    
    fManager_t M;
    M.SetFactory("Ptent", TentativePFact);
    //M.SetFactory("P",     SaPFact);
    //M.SetFactory("R",     RFact);

    M.SetFactory("Smoother", Teuchos::null);      //skips smoother setup
    M.SetFactory("CoarseSolver", Teuchos::null);  //skips coarsest solve setup

    int startLevel = 0;
    int maxLevels = 10;


    // Indicate which Hierarchy operators we want to keep
    hierarchy_mp = Teuchos::rcp(new hierarchy_t());  
    hierarchy_mp->Keep("P", SaPFact.get());  //SaPFact is the generating factory for P.
    hierarchy_mp->Keep("R", RFact.get());    //RFact is the generating factory for R.
    hierarchy_mp->Keep("Ptent", TentativePFact.get());  //SaPFact is the generating factory for P.

    hierarchy_mp->SetCycle(MueLu::CycleType::VCYCLE);
    hierarchy_mp->SetMaxCoarseSize(200);
    hierarchy_mp->SetPRrebalance(true);
    hierarchy_mp->SetImplicitTranspose(true);

    Teuchos::ParameterList Aclist = *(AFact->GetValidParameterList());
    Aclist.set("transpose: use implicit", true);
    AFact->SetParameterList(Aclist);
    
    // --------------------
    // The Factory Manager will be configured to return the rebalanced versions of P, R, A by default.
    // Everytime we want to use the non-rebalanced versions, we need to explicitly define the generating factory.
    RFact->SetFactory("P", SaPFact);
    //
    AFact->SetFactory("P", SaPFact);
    AFact->SetFactory("R", RFact);

    // Transfer coordinates
    Teuchos::RCP<coordTransferFactory_t> TransferCoordinatesFact = Teuchos::rcp(new coordTransferFactory_t());
    AFact->AddTransferFactory(TransferCoordinatesFact); // FIXME REMOVE

    // Repartitioning (creates "Importer" from "Partition")
    Teuchos::RCP<repartFactory_t> RepartitionFact = Teuchos::rcp(new repartFactory_t());
    {
        Teuchos::ParameterList paramList;
        paramList.set("repartition: start level", 2);
        paramList.set("repartition: min rows per proc", 200);
        paramList.set("repartition: max imbalance", 1.3);
        RepartitionFact->SetParameterList(paramList);
    }
    RepartitionFact->SetFactory("A", AFact);


    Teuchos::RCP<MueLu::Factory> ZoltanFact = Teuchos::rcp(new zoltan_t());
    ZoltanFact->SetFactory("A", AFact);
    ZoltanFact->SetFactory("Coordinates", TransferCoordinatesFact);
    ZoltanFact->SetFactory("number of partitions", RepartitionFact);

    // Reordering of the transfer operators
    Teuchos::RCP<rebalTransferFactory_t> RebalancedPFact = Teuchos::rcp(new rebalTransferFactory_t());
    RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
    RebalancedPFact->SetFactory("P", SaPFact);
    RebalancedPFact->SetFactory("Coordinates", TransferCoordinatesFact);
    RebalancedPFact->SetFactory("Nullspace", M.GetFactory("Ptent")); // TODO

    Teuchos::RCP<rebalTransferFactory_t> RebalancedRFact = Teuchos::rcp(new rebalTransferFactory_t());
    RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedRFact->SetFactory("R", RFact);

    // Compute Ac from rebalanced P and R
    Teuchos::RCP<rebalAcFactory_t> RebalancedAFact = Teuchos::rcp(new rebalAcFactory_t());
    RebalancedAFact->SetFactory("A", AFact);

    // Configure FactoryManager
    M.SetFactory("A", RebalancedAFact);
    M.SetFactory("P", RebalancedPFact);
    M.SetFactory("R", RebalancedRFact);
    M.SetFactory("Nullspace",   RebalancedPFact);
    M.SetFactory("Coordinates", RebalancedPFact);
    //M.SetFactory("Importer",    RepartitionFact);

    // -------------------

    //M.SetFactory("A", AFact);

    hierarchy_mp->IsPreconditioner(false);
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    // Create Gauss-Seidel smoother.                                                                                  
    std::string ifpackType = "RELAXATION";
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (lo_t) 4);
    ifpackList.set("relaxation: damping factor", (scalar_t) 1.0);
    Teuchos::RCP<smootherPrototype_t> smootherPrototype = Teuchos::rcp(new smoother_t(ifpackType, ifpackList));

    M.SetFactory("Smoother", Teuchos::rcp(new smootherFactory_t(smootherPrototype)));

    // Create coarsest solver.
    Teuchos::RCP<smootherFactory_t> coarseSolverPrototype = Teuchos::rcp( new solver_t() );
    Teuchos::RCP<smootherFactory_t>   coarseSolverFact = Teuchos::rcp( new smootherFactory_t(coarseSolverPrototype, Teuchos::null) );
    M.SetFactory("CoarseSolver", coarseSolverFact);

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
   
    finest_mp = hierarchy_mp->GetLevel();
    finest_mp->Set("A", A_mp);
    finest_mp->Set("Coordinates", coordinates);
    
    Teuchos::RCP<mv_t> nullspace = Teuchos::rcp(new mv_t(A->getRowMap(), 1));
    Teuchos::RCP<xmv_t> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector(nullspace);
    nullspace->putScalar(1.0);
    finest_mp->Set("Nullspace", xnullspace);
    //finest_mp->Set("Coordinates", coordinates)

    hierarchy_mp->Setup(M,startLevel,maxLevels);
}


std::size_t MueLuBottomSolver::getNumIters() {
    return 1;
}


void MueLuBottomSolver::setupAggregation_m() {

}


void MueLuBottomSolver::setupSmoother_m() {


}


void MueLuBottomSolver::setupCoarser_m() {

}
