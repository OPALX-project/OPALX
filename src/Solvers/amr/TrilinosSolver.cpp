#include "TrilinosSolver.h"

#include <cmath>
#include <ctime>

extern Inform *gmsg;

Solver::Solver(Epetra_MpiComm& comm,
               std::string itsolver,
               bool verbose, 
               double tol,
               int maxit,
               int nBlocks,
               int nRecycleBlocks, 
               unsigned int nLhs,
               int nLevels,
               const Array<MultiFab*>& rhs,
               const Array<MultiFab*>& soln,
               const Real* hr_in,
               const Real* prob_lo_in,
               BoundaryPointList& xlo,
               BoundaryPointList& xhi,
               BoundaryPointList& ylo,
               BoundaryPointList& yhi):
    comm_m(comm),
    itsolver_m(itsolver),
    isReusingHierarchy_m(true),
    isReusingPreconditioner_m(false),
    nBlocks_m(nBlocks),
    nRecycleBlocks_m(nRecycleBlocks),
    nLhs_m(nLhs),
    MLList_m(),
    xlo_m(xlo),
    xhi_m(xhi),
    ylo_m(ylo),
    yhi_m(yhi),
    idxMap_m(1)
{
    idxMap_m[0].reset(new iMultiFab(rhs[0]->boxArray(), 1 , 0));
    verbose_m = 0;
    
    
    for (int i = 0; i < 3; i++) {
        hr[i] = hr_in[i];
        prob_lo[i] = prob_lo_in[i];
    }
    
    // setup ml preconditioner parameters
    setupMultiLevelList();
    
    // get the problem dimension ( just single level support )
    BoxArray ba(rhs[0]->boxArray());
    domain_m = ba.minimalBox();
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        // +1 is needed since for e.g. for a 32x32x32 mesh: loVect = 0 and hiVect = 31
        nGridPoints_m[i] = domain_m.hiVect()[i] - domain_m.loVect()[i] + 1;
    }
    
    SetupProblem(rhs,soln);
    
    MLPrec = Teuchos::null;
    
    // setup Belos parameters
    setupBelosList_m(maxit, tol);
    
    
    // setup Belos solver
    if (nBlocks_m == 0 || nRecycleBlocks_m == 0) {
	if (itsolver_m == "CG")
	    solver = Teuchos::rcp( new Belos::BlockCGSolMgr<double,MV,OP>() );
	else if (itsolver_m == "GMRES")
	    solver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>() );
    } else
	solver = Teuchos::rcp( new Belos::RCGSolMgr<double,MV,OP>() );
    
    convStatusTest = Teuchos::rcp( new Belos::StatusTestGenResNorm<ST,MV,OP> (tol) );
    convStatusTest->defineScaleForm(Belos::NormOfRHS, Belos::TwoNorm);
    
#ifdef UserConv
    solver->setUserConvStatusTest(convStatusTest);
#endif
}

Solver::~Solver() {}


void
Solver::SetupProblem(const Array<MultiFab*>& rhs , const Array<MultiFab*>& soln) {
    // only single level support
    
    // get my global indices
    std::vector<int> gidx;
    fillGlobalIndex_m(gidx, *rhs[0 /*level*/]);
    
    /*
     * construct the mapping for the linear system of equations
     */
    // total number of elements
    int N = 1;
    for (int i = 0; i < BL_SPACEDIM; ++i)
        N *= nGridPoints_m[i];
    
    // Epetra_Map(int NumGlobalElements, int NumMyElements, const int* MyGlobalElements,
    //            int IndexBase, const Epetra_Comm& Comm)
    Epetra_Map map(N, gidx.size(), &gidx[0], 0, comm_m);
    
    
    if (nLhs_m > 0)
	P = Teuchos::rcp(new Epetra_MultiVector(/**Map*/map, nLhs_m, false));
    
    
    /*
     * fill right-hand side vector (false: do not initialize with default value zero)
     */
    
    RHS = Teuchos::rcp( new Epetra_Vector(map, false) );
    
    fillRHS_m(*rhs[0/*level*/]);
    
    /*
     * fill left-hand side vector (initial guess is 1.0 everywhere)
     */
    LHS = Teuchos::rcp( new Epetra_Vector(map, false) );
    
    /*
     * fill system matrix (7 elements per row)
     */
    fillSysMatrix_m(map);
}

void Solver::extrapolateLHS() {
// Aitken-Neville
// Pi0 (x) := yi , i = 0 : n
// Pik (x) := (x − xi ) Pi+1,k−1(x) − (x − xi+k ) Pi,k−1(x) /(xi+k − xi )
// k = 1, . . . , n, i = 0, . . . , n − k.

    if ( !nLhs_m ) {
	// default initial guess
	LHS->PutScalar(1.0);
	return;
    } else if ( nLhs_m < 0 )
	throw OpalException("Error in TrilinosSolver::extrapolateLHS",
			    "Invalid number of previous lhs solutions");

    std::deque< Epetra_Vector >::iterator it = OldLHS.begin();

    switch ( OldLHS.size() ) {

    case 1:
	*LHS = *it;
	break;
    case 2:
	LHS->Update (2.0, *it++, -1.0, *it, 0.0);
	break;
    default:
	// > 2 previous solutions
	int n = OldLHS.size();
        for (int i = 0; i < n; ++i)
            *(*P)(i) = *it++;
	
        for (int k = 1; k < n; ++k) // x==0 & xi==i+1
            for (int i = 0; i < n - k; ++i)
                (*P)(i)->Update( -(i + 1) / (float)k, *(*P)(i + 1), (i + k + 1) / (float)k);//TODO test

        *LHS = *(*P)(0);
	break;
    }
}

void Solver::CopySolution(const Array<MultiFab*>& soln)
{
    int lidx = 0;
    for (MFIter mfi(*soln[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox& lhs = (*soln[0])[mfi];
        IArrayBox& idx = (*idxMap_m[0])[mfi];
        
        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    lhs({i, j, k}) = LHS->Values()[ /*idx( {i, j, k} )*/lidx++ ];
                }
            }
        }
#ifdef DBG_SCALARFIELD
        /* Write the potential of all levels to file in the format: x, y, z, \phi
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = "data/amr-phi_scalar-level-" + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                if ( !out.is_open() )
                    throw OpalException("Error in TrilinosSolver::CopySolution",
                                        "Couldn't open the file: " + outfile);
                    
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " "
                                << lhs(ivec, 0) << " " << proc << std::endl;
                        }
                    }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
#endif
    }
}

int 
Solver::getNumIters() {
    return solver->getNumIters();
}

void 
Solver::Compute() {
    // estimate a good initial guess
    extrapolateLHS();

    // create the preconditioner object and compute hierarchy
    // true -> create the multilevel hirarchy
    // ML allows the user to cheaply recompute the preconditioner. You can
    // simply uncomment the following line:
    //
    // MLPrec->ReComputePreconditioner();
    //
    // It is supposed that the linear system matrix has different values, but
    // **exactly** the same structure and layout. The code re-built the
    // hierarchy and re-setup the smoothers and the coarse solver using
    // already available information on the hierarchy. A particular
    // care is required to use ReComputePreconditioner() with nonzero
    // threshold.
    
    if(MLPrec == Teuchos::null) { // first repetition we need to create a new preconditioner
        MLPrec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    } else if(isReusingHierarchy_m) {
        MLPrec->ReComputePreconditioner();
    } else if(isReusingPreconditioner_m) {
        // do nothing since we are reusing old preconditioner
    } else { // create a new preconditioner in every repetition
        delete MLPrec.get();//MLPrec now RCP => delete??? TODO
        MLPrec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    }
    
    // Setup problem
    problem.setOperator(A);
    problem.setLHS(LHS);
    problem.setRHS(RHS);
    prec = Teuchos::rcp(new Belos::EpetraPrecOp(MLPrec));
    problem.setLeftPrec(prec);
    solver->setParameters(Teuchos::rcp(&belosList_m, false));
    solver->setProblem(Teuchos::rcp(&problem,false));
    if(!problem.isProblemSet()){
        if ( !problem.setProblem() )
	  throw OpalException("Error in TrilinosSolver::Compute",
				"Belos::LinearProblem failed to set up correctly.");
    }

    // Solve problem
    solver->solve();

    // Store new LHS in OldLHS
    OldLHS.push_front(*(LHS.get()));
    if(OldLHS.size() > nLhs_m) OldLHS.pop_back();
}

// ------------------------------------------------------------------------------------------------------------
// PROTECTED MEMBER VARIABLES
// ------------------------------------------------------------------------------------------------------------

inline void Solver::setupMultiLevelList() {
    ML_Epetra::SetDefaults("SA", MLList_m);
    MLList_m.set("max levels", 8);
    MLList_m.set("increasing or decreasing", "increasing");

    // we use a V-cycle
    MLList_m.set("prec type", "MGV");

    // uncoupled aggregation is used (every processor aggregates
    // only local data)
    MLList_m.set("aggregation: type", "Uncoupled");

    // smoother related parameters
    MLList_m.set("smoother: type","Chebyshev");
    MLList_m.set("smoother: sweeps", 3);
    MLList_m.set("smoother: pre or post", "both");

    // on the coarsest level we solve with  Tim Davis' implementation of
    // Gilbert-Peierl's left-looking sparse partial pivoting algorithm,
    // with Eisenstat & Liu's symmetric pruning. Gilbert's version appears
    // as \c [L,U,P]=lu(A) in MATLAB. It doesn't exploit dense matrix
    // kernels, but it is the only sparse LU factorization algorithm known to be
    // asymptotically optimal, in the sense that it takes time proportional to the
    // number of floating-point operations.
    MLList_m.set("coarse: type", "Amesos-KLU");

    //XXX: or use Chebyshev coarse level solver
    // SEE PAPER FOR EVALUATION KLU vs. Chebyshev
    //MLList.set("coarse: sweeps", 10);
    //MLList.set("coarse: type", "Chebyshev");

    // turn on all output
    if(verbose_m)
        MLList_m.set("ML output", 101);
    else
        MLList_m.set("ML output", -1);

    // heuristic for max coarse size depending on number of processors
    int coarsest_size = std::max(comm_m.NumProc() * 10, 1024);
    MLList_m.set("coarse: max size", coarsest_size);
}


// -------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -------------------------------------------------------------------------------------------------------------

void Solver::setupBelosList_m(int maxit, double tol) {
    belosList_m.set("Maximum Iterations", maxit);    // Maximum number of iterations allowed
    belosList_m.set("Convergence Tolerance", tol);     // Relative convergence tolerance requested

    if(nBlocks_m != 0 && nRecycleBlocks_m != 0) {    // only set if solver==RCGSolMgr
	belosList_m.set("Num Blocks", nBlocks_m);      // Maximum number of blocks in Krylov space
	belosList_m.set("Num Recycled Blocks", nRecycleBlocks_m); // Number of vectors in recycle space
    }

    if(verbose_m) {
	belosList_m.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);
	belosList_m.set("Output Frequency", 1);
    } else
	belosList_m.set("Verbosity", Belos::Errors + Belos::FinalSummary);
}


int Solver::coordBox2globalIdx_m(const IntVect& coord) {
    return coord[0] + nGridPoints_m[0] * coord[1] + nGridPoints_m[0] * nGridPoints_m[1] * coord[2];
}


void Solver::fillGlobalIndex_m(std::vector<int>& gidx, const MultiFab& rhs) {
    
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        IArrayBox& ifab = (*idxMap_m[0])[mfi];
        const Box& bx = mfi.validbox();
        
        // iterate over local indices given by the box boundary
        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    IntVect coord(i, j, k);
                    ifab(coord) = coordBox2globalIdx_m( coord );
                    gidx.push_back( ifab(coord) );
                }
            }
        }
    }
}

void Solver::fillRHS_m(const MultiFab& rhs) {
    
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        const FArrayBox& fab = rhs[mfi];
        const Box& bx = mfi.validbox();
        
        // local index
        int lidx = 0;
        
        // iterate over local indices given by the box boundary
        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    RHS->Values()[lidx++] = fab({i, j, k});
                }
            }
        }
    }
}

void Solver::fillSysMatrix_m(const Epetra_Map& map) {
    // not nice coding
    
    A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, map, 7) );
    
    double values[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int indices[7] = { 0, 0, 0, 0, 0, 0, 0};
    int nEntries = 0;
    
    for (MFIter mfi(*idxMap_m[0]); mfi.isValid(); ++mfi) {
        IArrayBox& ifab = (*idxMap_m[0])[mfi];
        const Box& bx = mfi.validbox();
        
        // iterate over local indices given by the box boundary
        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
            for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    // ordering: center, west, east, north, south, front, back
                    getStencil({i, j, k}, &values[0], &indices[0], nEntries, ifab);
                    int row = ifab( {i, j, k} );
                    A->InsertGlobalValues(row , nEntries, &values[0], &indices[0]);
                }
            }
        }
    }
    
    A->FillComplete();
    
    *gmsg << "Number of non-zeros: " << A->NumGlobalNonzeros () << endl;
}

void Solver::getStencil(const IntVect& coord, double* values, int* indices, int& nEntries, const IArrayBox& ifab) {
    // not nice coding
    
    
    nEntries = 0;
    
    
    int center = ifab( coord );
    
    values[nEntries] = - 2.0 / ( hr[0] * hr[0]) - 2.0 / ( hr[1] * hr[1] ) - 2.0 / ( hr[2] * hr[2] );
    indices[nEntries++] = ifab( coord );
    
    
    int nx = nGridPoints_m[0];
    int ny = nGridPoints_m[1];
    int nz = nGridPoints_m[2];
    
    
    // west
    if ( center % nx != 0 ) {
        values[nEntries] = 1.0 / (hr[0] * hr[0]);
        indices[nEntries++] = center - 1;
    }
    
    // east
    if ( (center + 1) % nx != 0 && center + 1 < nx * ny * nz ) {
        values[nEntries] = 1.0 / (hr[0] * hr[0]);
        indices[nEntries++] = center + 1;
    }
    
    // south
    if ( center - nx >= 0 && center % (nx * ny) >= nx ) {
        values[nEntries] = 1.0 / (hr[1] * hr[1]);
        indices[nEntries++] = center - nx;
    }
    
    // north
    if ( (center + nx) % (nx * ny) > ny - 1 && (center + nx) < nx * ny * nz ) {
        values[nEntries] = 1.0 / (hr[1] * hr[1]);
        indices[nEntries++] = center + nx;
    }
    
    // front
    if ( center - nx * ny >= 0 ) {
        values[nEntries] = 1.0 / (hr[2] * hr[2]);
        indices[nEntries++] = center - nx * ny;
    }
    
    // back
    if ( center + nx * ny < nx * ny * nz ) {
        values[nEntries] = 1.0 / (hr[2] * hr[2]);
        indices[nEntries++] = center + nx * ny;
    }
}
