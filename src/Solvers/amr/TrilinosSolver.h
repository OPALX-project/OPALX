#ifndef _TRILINOS_SOLVER_H
#define _TRILINOS_SOLVER_H

#include "ml_include.h"

#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include <iMultiFab.H>

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"
#include <Epetra_Import.h>

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosStatusTestGenResNorm.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_ParameterList.hpp>

#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

#include <vector>
#include <deque>
#include <fstream>
#include <sstream>
#include <iostream>

/// enumeration corresponding to different interpolation methods at the boundary
enum {
    CONSTANT,
    LINEAR,
    QUADRATIC
};

enum {
    INSIDE,
    OUTSIDE,
    COARSER,
    COVERED
};

typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

//using namespace Teuchos;
// using Teuchos::RCP;
// using Teuchos::rcp;
// using namespace ML_Epetra;
// using namespace Isorropia;

class Solver {

public:
    /** Constructor for the SAAMGPCG solver.
     * \param Comm the communicator object used in the solver routines.
     * \param verbose flag if output should be verbose.
     * \sa ~Solver(), SetupMLList() and SetupProblem()
     */
    Solver(Epetra_MpiComm& Comm, std::string& itsolver, std::string& interpl, bool verbose,
           double tol, int maxIterations, int numBlocks, int recycleBlocks,
           int maxOldLHS, int nlevs, PArray<MultiFab>& rhs, PArray<MultiFab>& soln,
           const Real* hr_in,  const Real* prob_lo_in,
           BoundaryPointList& xlo_in, BoundaryPointList& xhi_in,
           BoundaryPointList& ylo_in, BoundaryPointList& yhi_in) :
          Comm_m(Comm), mask(nlevs), imap(nlevs)
          {
            Diag = 0;

//          verbose_m = verbose;
            verbose_m = 0;

            for (int i = 0; i < 3; i++)
            {
                hr[i] = hr_in[i];
                prob_lo[i] = prob_lo_in[i];
            }

            xlo = xlo_in;
            xhi = xhi_in;
            ylo = ylo_in;
            yhi = yhi_in;

            if (interpl == "constant")
                interpolationMethod = CONSTANT;
            else if (interpl == "linear")
                interpolationMethod = LINEAR;
            else if (interpl == "quadratic")
                interpolationMethod = QUADRATIC;

	    itsolver_m = itsolver;

            isReusingHierarchy_m = false;
            isReusingPreconditioner_m = false;
            //XXX: currently set reusing to hierarchy
            isReusingHierarchy_m = true;

            numBlocks_m = numBlocks;
            recycleBlocks_m = recycleBlocks;
            nLHS_m = maxOldLHS;

            // setup ml preconditioner parameters
            SetupMLList();

            // First construct the map and mask
            for (int i = 0; i < nlevs; i++)
            {
               // Create new imultifabs with 1 component and 1 ghost cell
               mask.set(i,new iMultiFab(rhs[i].boxArray(),1,1));
               imap.set(i,new iMultiFab(rhs[i].boxArray(),1,1));
            }

            construct_mask_and_imap();

            SetupProblem(rhs,soln);

            // setup extrapolation
	    if(nLHS_m > 0)
		    P = rcp(new Epetra_MultiVector(*Map, nLHS_m, false));

	    MLPrec = Teuchos::null;

	    // setup Belos parameters
	    belosList.set( "Maximum Iterations", maxIterations );  // Maximum number of iterations allowed
	    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
	    if(numBlocks_m != 0 && recycleBlocks_m != 0){          // only set if solver==RCGSolMgr
		    belosList.set("Num Blocks", numBlocks_m);          // Maximum number of blocks in Krylov space
		    belosList.set("Num Recycled Blocks", recycleBlocks_m); // Number of vectors in recycle space
	    }
	    if(verbose_m) {
		    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);
		    belosList.set("Output Frequency", 1);
	    } else
		   // belosList.set("Verbosity", Belos::Errors + Belos::Warnings);
		    belosList.set("Verbosity", Belos::Errors + Belos::FinalSummary);

	    // setup Belos solver
	    if(numBlocks_m == 0 || recycleBlocks_m == 0)
	    {
		if(itsolver_m == "CG")
		    solver = rcp( new Belos::BlockCGSolMgr<double,MV,OP>() );
		else if (itsolver_m == "GMRES")
		    solver = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>() );
	    }
	    else
		    solver = rcp( new Belos::RCGSolMgr<double,MV,OP>() );
	    convStatusTest = rcp( new Belos::StatusTestGenResNorm<ST,MV,OP> (tol) );
	    convStatusTest->defineScaleForm(Belos::NormOfRHS, Belos::TwoNorm);
#ifdef UserConv
	    solver->setUserConvStatusTest(convStatusTest);
#endif
    }
      /** Destructor for the SAAMGPCG solver.
     */
    ~Solver() {
	    delete Diag;
	    delete Map;
    }

    /** Actual computation of self field.
     * \param repetitions the number of calls to the solver.
     * \param tol the convergence criterion for the PCG solver.
     * \param maxIterations the maximal number of iterations performed by the PCG solver.  */
    void Compute();

    /** Performs setup of data distribution and problem (system matrix, RHS, LHS).  */
    void SetupProblem(PArray<MultiFab>& rhs, PArray<MultiFab>& soln);

    /** Copy LHS->Values() into MultiFab soln */
    void CopySolution(PArray<MultiFab>& soln);

    int getNumIters();

private:

    Real hr[3];
    Real prob_lo[3];

    BoundaryPointList xlo;
    BoundaryPointList xhi;
    BoundaryPointList ylo;
    BoundaryPointList yhi;

    int* proc_offset;

    /// Map corresponding to data distribution
    Epetra_Map* Map;

    /// preconditioner object
    RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;//now RCP TODO

    /// parameter list used for the preconditioner (filled in SetupMLList)
    Teuchos::ParameterList MLList_m;

    /// system matrix
    RCP<Epetra_CrsMatrix> A;

    /// diagonal matrix of A
    Epetra_CrsMatrix* Diag;

    /// right hand side of our problem
    RCP<Epetra_Vector> RHS;

    /// left hand side
    RCP<Epetra_Vector> LHS;
    RCP<Epetra_MultiVector> lhssol;

    /// last N LHS's for extrapolating the new LHS as starting vector
    //uint nLHS_m;
    int nLHS_m;
    RCP<Epetra_MultiVector> P;
    std::deque< Epetra_Vector > OldLHS;

    /// Solver (Belos BlockCG or RCG)
    /// if numBlocks_m OR recycleBlocks_m == 0 BlockCG is used, else RCG
    /// maximum number of blocks in Krylov space
    int numBlocks_m;
    /// number of vectors in recycle space
    int recycleBlocks_m;

    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST,MV,OP> OPT;// TODO
    typedef Belos::MultiVecTraits<ST,MV>    MVT;// TODO

    Teuchos::ParameterList belosList;
    Belos::LinearProblem<ST,MV,OP> problem;
    RCP< Belos::EpetraPrecOp > prec;//TODO just wrapper for MLPrec->LinearProblem
    RCP< Belos::StatusTestGenResNorm< ST, MV, OP > > convStatusTest;
	RCP< Belos::SolverManager<ST,MV,OP> > solver;

    /// communicator
    Epetra_MpiComm Comm_m;

    /// iterativer solver type
    std::string	 itsolver_m;
    /// interpolation type
    int interpolationMethod;

    /// verbosity of the code
    bool verbose_m;

    /// flag specifying if hierarchy is reused
    bool isReusingHierarchy_m;

    /// flag specifying if whole preconditioner is reued
    bool isReusingPreconditioner_m;

    /// rectangular data distribution
    int my_slide_size[3], my_start_offset[3], my_end_offset[3];

    /** Performs a redistribution of the data with RCB to avoid idle processors.
     */
    void redistributeWithRCB();

    void getBoundaryStencil(int x, int y, int z, double& WV, double& EV, double& SV, double& NV,
                            double& FV, double& BV, double& CV, double& scaleFactor, IArrayBox& fab_mask)
    {
    // Determine which interpolation method we use for points near the boundary
    if (interpolationMethod == CONSTANT)
        ConstantInterpolation(x,y,z,WV,EV,SV,NV,FV,BV,CV,scaleFactor,fab_mask);
#if 0
    else if (interpolationMethod == LINEAR)
        LinearInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
    else if (interpolationMethod == QUADRATIC)
        QuadraticInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
#endif
    // stencil center value has to be positive!
    assert(CV > 0);
    };

    void ConstantInterpolation(int x, int y, int z, double& WV, double& EV, double& SV, double& NV,
                               double& FV, double& BV, double& CV, double& scaleFactor,
                               IArrayBox& fab_mask)
    {

    scaleFactor = 1.0;

    WV = -1/(hr[0]*hr[0]);
    EV = -1/(hr[0]*hr[0]);
    NV = -1/(hr[1]*hr[1]);
    SV = -1/(hr[1]*hr[1]);
    FV = -1/(hr[2]*hr[2]);
    BV = -1/(hr[2]*hr[2]);
    CV = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    IntVect iv(x,y,z);

    // we are a left boundary point
    if (fab_mask(IntVect(x-1,y,z),0) == OUTSIDE)
        WV = 0.0;

    // we are a right boundary point
    if (fab_mask(IntVect(x+1,y,z),0) == OUTSIDE)
        EV = 0.0;

    // we are a lower boundary point
    if (fab_mask(IntVect(x,y-1,z),0) == OUTSIDE)
        SV = 0.0;

    // we are a upper boundary point
    if (fab_mask(IntVect(x,y+1,z),0) == OUTSIDE)
        NV = 0.0;

    // we are at the low z end
    if (fab_mask(IntVect(x,y,z-1),0) == OUTSIDE)
        FV = 0.0;

    // we are at the high z end
    if (fab_mask(IntVect(x,y,z+1),0) == OUTSIDE)
        BV = 0.0;
};

    void printLoadBalanceStats();

    void extrapolateLHS();

    PArray<iMultiFab> mask;
    PArray<iMultiFab> imap;

    void construct_mask_and_imap();

    void avgDownImap(iMultiFab& crse_imap, iMultiFab& fine_imap);

    void FillImapGhostCells(iMultiFab& crse_imap, iMultiFab& fine_imap, iMultiFab& fine_mask);

    void ComputeStencil();

protected:

    /** Setup the parameters for the SAAMG preconditioner.
     */
    inline void SetupMLList() {
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
        int coarsest_size = std::max(Comm_m.NumProc() * 10, 1024);
        MLList_m.set("coarse: max size", coarsest_size);
    }

};

#endif
