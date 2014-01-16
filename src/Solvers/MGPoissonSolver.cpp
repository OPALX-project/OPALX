#ifdef HAVE_ML_SOLVER
#define DBG_STENCIL
#include "Algorithms/PartBunch.h"
#include "MGPoissonSolver.h"
#include "Physics/Physics.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include <algorithm>

using Physics::c;

MGPoissonSolver::MGPoissonSolver(PartBunch &beam,Mesh_t *mesh, FieldLayout_t *fl, std::vector<BoundaryGeometry *> geometries, std::string itsolver, std::string interpl, double tol, int maxiters, std::string precmode):
    itsBunch_m(&beam),
    mesh_m(mesh), 
    layout_m(fl), 
    geometries_m(geometries), 
    tol_m(tol), 
    maxiters_m(maxiters), 
    Comm(Ippl::getComm()) {
    domain_m = layout_m->getDomain();

    for(int i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        orig_nr_m[i] = domain_m[i].length();
    }

    if(itsolver == "CG") itsolver_m = AZ_cg;
    else if(itsolver == "BiCGStab") itsolver_m = AZ_bicgstab;
    else if(itsolver == "GMRES") itsolver_m = AZ_gmres;
    else throw OpalException("MGPoissonSolver", "No valid iterative solver selected!");

    precmode_m = STD_PREC;
    if(precmode == "STD") precmode_m = STD_PREC;
    else if(precmode == "HIERARCHY") precmode_m = REUSE_HIERARCHY;
    else if(precmode == "REUSE") precmode_m = REUSE_PREC;

    tolerableIterationsCount_m = 2;
    numIter_m = -1;
    forcePreconditionerRecomputation_m = false;

    hasParallelDecompositionChanged_m = true;
    useRCB_m = false;
    if(Ippl::Info->getOutputLevel() > 0)
        verbose_m = true;
    else
        verbose_m = false;

    // Find CURRENT geometry
    currentGeometry = geometries_m[0];
    if(currentGeometry->getFilename() == "") {
        if(currentGeometry->getTopology() == "ELLIPTIC")
            bp = new EllipticDomain(currentGeometry->getA(), currentGeometry->getB(), orig_nr_m, hr_m, interpl);
        else if(currentGeometry->getTopology() == "BOXCORNER") {
            bp = new BoxCornerDomain(currentGeometry->getA(), currentGeometry->getB(), currentGeometry->getC(), currentGeometry->getLength(),
                                     currentGeometry->getL1(), currentGeometry->getL2(), orig_nr_m, hr_m, interpl);
            bp->Compute(itsBunch_m->get_hr());
        } else {
            ERRORMSG("Geometry not known" << endl);
            exit(1);
        }
    } else {
	cout << "hr_m = " << hr_m << endl;
        localidx_m =  layout_m->getLocalNDIndex();
	//bp = new ArbitraryDomain(currentGeometry, currentGeometry->getmincoords(), currentGeometry->getmaxcoords(), orig_nr_m, hr_m, localidx_m, interpl);
    }
#if 0
    Map = 0;
    A = Teuchos::null;
    LHS = Teuchos::null;
    RHS = Teuchos::null;
#endif
    MLPrec = Teuchos::null;
    numBlocks_m = Options::numBlocks;
    recycleBlocks_m = Options::recycleBlocks;
    nLHS = Options::nLHS;
    SetupMLList();
    SetupBelosList();

    // setup Belos solver
    cout<<"numBlocks_m "<<numBlocks_m<<" " <<"recycleBlocks_m "<<recycleBlocks_m<<endl;

    if(numBlocks_m == 0 || recycleBlocks_m == 0)
        solver = rcp(new Belos::BlockCGSolMgr<double, MV, OP>());
    else
        solver = rcp(new Belos::RCGSolMgr<double, MV, OP>());
    convStatusTest = rcp(new Belos::StatusTestGenResNorm<ST, MV, OP> (tol));
    convStatusTest->defineScaleForm(Belos::NormOfRHS, Belos::TwoNorm);
#if defined(setUserConv)
    solver->setUserConvStatusTest(convStatusTest);
#endif

    //all timers used here
    FunctionTimer1_m = IpplTimings::getTimer("Create Stencil");
    FunctionTimer2_m = IpplTimings::getTimer("Prepare RHS");
    FunctionTimer3_m = IpplTimings::getTimer("CG");
    FunctionTimer4_m = IpplTimings::getTimer("LHS to Rho");
    FunctionTimer5_m = IpplTimings::getTimer("recalculate Map");
    FunctionTimer6_m = IpplTimings::getTimer("ML");
}


MGPoissonSolver::MGPoissonSolver(PartBunch &beam):
    layout_m(&beam.getFieldLayout()), mesh_m(&beam.getMesh()), itsBunch_m(&beam), LHS(0), Comm(Ippl::getComm()) {

    throw OpalException("MGPoissonSolver", "Copy Constructor not implemented yet");
}

MGPoissonSolver::~MGPoissonSolver() {
    if(Map) delete Map;
}

void MGPoissonSolver::computePotential(Field_t &rho, Vector_t hr, double zshift) {
    throw OpalException("MGPoissonSolver", "not implemented yet");
}

void MGPoissonSolver::computeMap() {

    IpplTimings::startTimer(FunctionTimer5_m);
    NDIndex<3> localidx =  layout_m->getLocalNDIndex();

    //FIXME: we have to set this to true after constructor has been called
    //(otherwise Map 'size' is 0)
    //if(hasParallelDecompositionChanged_m == true) {
    if(true) {
        if(Map != 0)
            delete Map;
        if(useRCB_m)
            redistributeWithRCB(localidx);
        else
            IPPLToMap3D(localidx);
#if 0
        //we also have to redistribute LHS
        if(LHS == Teuchos::null) {
		std::cout << "COMPUTEMAP LHS->PutScalar" <<std::endl;
            LHS = rcp(new Epetra_Vector(*Map));
            LHS->PutScalar(1.0);
        } else {
		std::cout << "COMPUTEMAP LHS->PutScalar NO" <<std::endl;
            RCP<Epetra_Vector> tmplhs = rcp(new Epetra_Vector(*Map));
            Epetra_Import importer(*Map, LHS->Map());
            tmplhs->Import(*LHS, importer, Add);
            LHS = tmplhs;
        }

        //...and all previously saved LHS
        std::deque< Epetra_Vector >::iterator it = OldLHS.begin();
        if(OldLHS.size() > 0) {
            int n = OldLHS.size();
            for(int i = 0; i < n; ++i) {
                Epetra_Vector tmplhs = Epetra_Vector(*Map);
                Epetra_Import importer(*Map, it->Map());
                tmplhs.Import(*it, importer, Add);
                *it = tmplhs;
                ++it;
            }
        }
#endif
    }
    IpplTimings::stopTimer(FunctionTimer5_m);
}

void MGPoissonSolver::extrapolateLHS() {
// Aitken-Neville
// Pi0 (x) := yi , i = 0 : n
// Pik (x) := (x − xi ) Pi+1,k−1(x) − (x − xi+k ) Pi,k−1(x) /(xi+k − xi )
// k = 1, . . . , n, i = 0, . . . , n − k.

    std::deque< Epetra_Vector >::iterator it = OldLHS.begin();

    if(nLHS == 0)
        LHS->PutScalar(1.0);
    else if(OldLHS.size() == 1)
        *LHS = *it;
    else if(OldLHS.size() == 2){
        LHS->Update (2.0, *it++, -1.0, *it, 0.0);
    }
    else if(OldLHS.size() > 0)
    {
        int n = OldLHS.size();
        for(int i=0; i<n; ++i){
            *(*P)(i) = *it++;
        }
        for(int k = 1; k < n; ++k){// x==0 & xi==i+1
            for(int i = 0; i < n-k; ++i){
                (*P)(i)->Update(-(i+1)/(float)k, *(*P)(i+1), (i+k+1)/(float)k);//TODO test
            }
        }
        *LHS = *(*P)(0);
    }
    else
    {
        std::cout << "Invalid number of old LHS: " + OldLHS.size() << std::endl;
    }
}


// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation
// XXX: use matrix stencil in computation directly (no Epetra, define operators
// on IPPL GRID)
void MGPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {
    //TODO: needed?
    nr_m[0] = orig_nr_m[0];
    nr_m[1] = orig_nr_m[1];
    nr_m[2] = orig_nr_m[2];

    const double scaleFactor = itsBunch_m->getdT() * Physics::c;
//    hr *= Vector_t(scaleFactor); //m  
    bp->setMinMaxZ(itsBunch_m->get_origin()[2]*scaleFactor, itsBunch_m->get_maxExtend()[2]*scaleFactor);
    bp->setNr(nr_m);
    bp->Compute(hr);

    // Define the Map
    computeMap();

    // Allocate the RHS and LHS with the new Epetra Map
    RHS = rcp(new Epetra_Vector(*Map));
    LHS = rcp(new Epetra_Vector(*Map));
    
    // get charge densities from IPPL field and store in Epetra vector (RHS)
    NDIndex<3> localidx = layout_m->getLocalNDIndex();

    IpplTimings::startTimer(FunctionTimer2_m);
    int idx = 0;
    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

                if(bp->isInside(x, y, z)) {
                    RHS->Values()[idx++] = rho[x][y][z].get();
                }
            }
        }
    }
    IpplTimings::stopTimer(FunctionTimer2_m);

    extrapolateLHS();
    std::cout << "Done with extrapolateLHS " << std::endl; 
    std::cout << " .... and starting to solve " << std::endl; 


    // build discretization matrix
    IpplTimings::startTimer(FunctionTimer1_m);
    A = rcp(new Epetra_CrsMatrix(Copy, *Map,  7, true));
    ComputeStencil(hr, RHS);
    IpplTimings::stopTimer(FunctionTimer1_m);

#ifdef DBG_STENCIL
    EpetraExt::RowMatrixToMatlabFile("DiscrStencil.dat", *A);
#endif
    
    IpplTimings::startTimer(FunctionTimer6_m);

//    if(MLPrec == Teuchos::null) // first repetition we need to create a new preconditioner
        MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
	/*
    else{	
		std::cout << "HERE I AM " << std::endl;
	    switch(precmode_m) {
       		case REUSE_PREC:
               		MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
            		break;

	        case REUSE_HIERARCHY:
        	        MLPrec->ReComputePreconditioner();
            		break;

        	case STD_PREC:
	   		delete MLPrec.get();
            		MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
            		break;
	    }
    }
	*/
    IpplTimings::stopTimer(FunctionTimer6_m);

    // setup preconditioned iterative solver
    // use old LHS solution as initial guess
    problem.setOperator(A);
    problem.setLHS(LHS);
    problem.setRHS(RHS);
    prec = rcp(new Belos::EpetraPrecOp(MLPrec));
    problem.setLeftPrec(prec);
    solver->setParameters(rcp(&belosList, false));
    solver->setProblem(rcp(&problem, false));
    if(!problem.isProblemSet()) {
        if(problem.setProblem() == false)
            throw OpalException("MGPoissonSolver", "Belos::LinearProblem failed to set up correctly!");
    }

    std::ofstream timings;
    char filename[50];
    sprintf(filename, "timing_MX%d_MY%d_MZ%d_nProc%d_recB%d_numB%d_nLHS%d", orig_nr_m[0], orig_nr_m[1], orig_nr_m[2], Comm.NumProc(), recycleBlocks_m, numBlocks_m, nLHS);
    if(Comm.MyPID() == 0) timings.open(filename, std::ios::app);
    double time = MPI_Wtime();

    IpplTimings::startTimer(FunctionTimer3_m);
    	solver->solve();
    IpplTimings::stopTimer(FunctionTimer3_m);

    time = MPI_Wtime() - time;
    double minTime, maxTime, avgTime;
    Comm.MinAll(&time, &minTime, 1);
    Comm.MaxAll(&time, &maxTime, 1);
    Comm.SumAll(&time, &avgTime, 1);
    avgTime /= Comm.NumProc();
    if(Comm.MyPID() == 0) timings <<
                                      solver->getNumIters() << "\t" <<
                                      //time <<" "<<
                                      minTime << "\t" <<
                                      maxTime << "\t" <<
                                      avgTime << "\t" <<
                                      numBlocks_m << "\t" <<
                                      recycleBlocks_m << "\t" <<
                                      nLHS << "\t" <<
                                      //OldLHS.size() <<"\t"<<
                                      std::endl;
    if(Comm.MyPID() == 0) timings.close();


    //FIXME: what happens when iteration count changes due to geometry?
    //int numiter = solver.NumIters();
    //if(numIter_m != -1) numIter_m = numiter;
    //if(numiter > numIter_m + tolerableIterationsCount_m) {
    //numIter_m = numiter;
    //forcePreconditionerRecomputation_m == true;
    //delete MLPrec;
    //}

    // Store new LHS in OldLHS
	#if 0
    if(nLHS > 1) 
	#endif
    OldLHS.push_front(*(LHS.get()));
    if(OldLHS.size() > nLHS) OldLHS.pop_back();

    //now transfer solution back to IPPL grid
    IpplTimings::startTimer(FunctionTimer4_m);
    idx = 0;
    rho = 0.0;
    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

                NDIndex<3> l(Index(x, x), Index(y, y), Index(z, z));
                if(bp->isInside(x, y, z))
                    rho.localElement(l) = LHS->Values()[idx++];
            }
        }
    }
    IpplTimings::stopTimer(FunctionTimer4_m);

}


void MGPoissonSolver::redistributeWithRCB(NDIndex<3> localidx) {

    int numMyGridPoints = 0;

    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
                if(bp->isInside(x, y, z))
                    numMyGridPoints++;
            }
        }
    }

    Epetra_BlockMap bmap(-1, numMyGridPoints, 1, 0, Comm);
    Teuchos::RCP<const Epetra_MultiVector> coords = Teuchos::rcp(new Epetra_MultiVector(bmap, 3, false));

    double *v;
    int stride, stride2;

    coords->ExtractView(&v, &stride);
    stride2 = 2 * stride;

    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
                if(bp->isInside(x, y, z)) {
                    v[0] = (double)x;
                    v[stride] = (double)y;
                    v[stride2] = (double)z;
                    v++;
                }
            }
        }
    }

    Teuchos::ParameterList paramlist;
    paramlist.set("Partitioning Method", "RCB");
    Teuchos::ParameterList &sublist = paramlist.sublist("ZOLTAN");
    sublist.set("RCB_RECTILINEAR_BLOCKS", "1");
    sublist.set("DEBUG_LEVEL", "1");

    Teuchos::RCP<Isorropia::Epetra::Partitioner> part = Teuchos::rcp(new Isorropia::Epetra::Partitioner(coords, paramlist));
    Isorropia::Epetra::Redistributor rd(part);
    Teuchos::RCP<Epetra_MultiVector> newcoords = rd.redistribute(*coords);

    newcoords->ExtractView(&v, &stride);
    stride2 = 2 * stride;
    numMyGridPoints = 0;
    std::vector<int> MyGlobalElements;
    for(int i = 0; i < newcoords->MyLength(); i++) {
        MyGlobalElements.push_back(bp->getIdx(v[0], v[stride], v[stride2]));
        v++;
        numMyGridPoints++;
    }

    Map = new Epetra_Map(-1, numMyGridPoints, &MyGlobalElements[0], 0, Comm);
}

void MGPoissonSolver::IPPLToMap3D(NDIndex<3> localidx) {

    int NumMyElements = 0;
    vector<int> MyGlobalElements;

    for(int x = localidx[0].first(); x <= localidx[0].last(); x++)
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++)
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
                if(bp->isInside(x, y, z)) {
                    MyGlobalElements.push_back(bp->getIdx(x, y, z));
                    NumMyElements++;
                }
            }
	std::cout <<"NumMyElements" << NumMyElements << std::endl;
    Map = new Epetra_Map(-1, NumMyElements, &MyGlobalElements[0], 0, Comm);
}

void MGPoissonSolver::ComputeStencil(Vector_t hr, Teuchos::RCP<Epetra_Vector> RHS) {

    int NumMyElements = Map->NumMyElements();
    int *MyGlobalElements = Map->MyGlobalElements();

    vector<double> Values(6);
    vector<int> Indices(6);

    //XXX: can we use x,y,z and remove idx->coord mapping?
    for(int i = 0 ; i < NumMyElements ; ++i) {

        int NumEntries = 0;

        double WV, EV, SV, NV, FV, BV, CV, scaleFactor; 
        int W, E, S, N, F, B;

        bp->getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, FV, BV, CV, scaleFactor);
        RHS->Values()[i] *= scaleFactor; 

        bp->getNeighbours(MyGlobalElements[i], W, E, S, N, F, B);

        if(E != -1) {
            Indices[NumEntries] = E;
            Values[NumEntries++] = EV;
        }
        if(W != -1) {
            Indices[NumEntries] = W;
            Values[NumEntries++] = WV;
        }
        if(S != -1) {
            Indices[NumEntries] = S;
            Values[NumEntries++] = SV;
        }
        if(N != -1) {
            Indices[NumEntries] = N;
            Values[NumEntries++] = NV;
        }
        if(F != -1) {
            Indices[NumEntries] = F;
            Values[NumEntries++] = FV;
        }
        if(B != -1) {
            Indices[NumEntries] = B;
            Values[NumEntries++] = BV;
        }

        // if matrix has already been filled (FillComplete()) we can only
        // replace entries

            // off-diagonal entries
            A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
            // diagonal entry
            A->InsertGlobalValues(MyGlobalElements[i], 1, &CV, &MyGlobalElements[i]);

    }

    A->FillComplete();
    A->OptimizeStorage();
}

void MGPoissonSolver::printLoadBalanceStats() {

    //compute some load balance statistics
    size_t myNumPart = Map->NumMyElements();
    size_t NumPart = Map->NumGlobalElements() * 1.0 / Comm.NumProc();
    double imbalance = 1.0;
    if(myNumPart >= NumPart)
        imbalance += (myNumPart - NumPart) / NumPart;
    else
        imbalance += (NumPart - myNumPart) / NumPart;

    double max = 0.0, min = 0.0, avg = 0.0;
    int minn = 0, maxn = 0;
    MPI_Reduce(&imbalance, &min, 1, MPI_DOUBLE, MPI_MIN, 0, Ippl::getComm());
    MPI_Reduce(&imbalance, &max, 1, MPI_DOUBLE, MPI_MAX, 0, Ippl::getComm());
    MPI_Reduce(&imbalance, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, Ippl::getComm());
    MPI_Reduce(&myNumPart, &minn, 1, MPI_INT, MPI_MIN, 0, Ippl::getComm());
    MPI_Reduce(&myNumPart, &maxn, 1, MPI_INT, MPI_MAX, 0, Ippl::getComm());

    avg /= Comm.NumProc();
    if(Comm.MyPID() == 0) cout << "LBAL min = " << min << ", max = " << max << ", avg = " << avg << endl;
    if(Comm.MyPID() == 0) cout << "min nr gridpoints = " << minn << ", max nr gridpoints = " << maxn << endl;
}

Inform &MGPoissonSolver::print(Inform &os) const {
    os << "* *************** M G P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
    return os;	
}

#endif /* HAVE_ML_SOLVER */
