#ifndef _TRILINOS_SOLVER_H
#define _TRILINOS_SOLVER_H

// ------------------------------------------------------------------------
// $RCSfile: TrilinosSolver.h,v $
// ------------------------------------------------------------------------
// $Revision: 0.0.0.0 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrilinosSolver
//
// ------------------------------------------------------------------------
//
// $Date: unknown $
// $Author: Tulin Kaman $
//
// ------------------------------------------------------------------------


// BoxLib headers
#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include <iMultiFab.H>

// Trilinos headers
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Operator.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MatrixMatrix.h>
#include <Epetra_Import.h>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosRCGSolMgr.hpp>
#include <BelosStatusTestGenResNorm.hpp>

#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_MultiLevelOperator.h>
#include <ml_epetra_utils.h>

#include <vector>
#include <deque>
#include <fstream>
#include <sstream>

// OPAL headers
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"


typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

using namespace ML_Epetra;

// Class TrilinosSolver
// ------------------------------------------------------------------------
/// This class solves a linear system of equations coming from AMR using SAAMGPCG.
class Solver {
    
private:
    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST,MV,OP> OPT;// TODO
    typedef Belos::MultiVecTraits<ST,MV>    MVT;// TODO

public:
    /// The different interpolation methods at the boundary
    enum {
	CONSTANT,
	LINEAR,
	QUADRATIC
    };

    /// The cell property of AMR (used to convert between Trilinos and BoxLib)
    enum {
	INSIDE,
	OUTSIDE,
	COARSER,
	COVERED
    };

    /*! Constructor for the SAAMGPCG solver.
     * @param comm is the communicator object used in the solver routines.
     * @param itsolver specifies the type of solver (e.g. CG, GMRES)
     * @param interpl is the interpolation scheme at the boundary.
     * @param vrebose flat if output should be verbose.
     * @param tol is the tolerance of the solver.
     * @param maxit is the maximum number of allowed iterations for the solver.
     * @param nRecycleBlocks is the number of vectors in recycle space.
     * @param nLhs is the max. number of lhs solutions that are kept to compute a new initial guess.
     * @param nLevels is the max. number of refinement levels.
     * @param rhs stores the data of the density for all levels
     * @param soln will be filled with the solution.
     * @param hr_in is the spacing in all directions
     * @param prob_lo_in is the lower boundary of the problem domain
     * @param xlo
     * @param xhi
     * @param ylo
     * @param yhi
     */
    Solver(Epetra_MpiComm& comm,
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
	   BoundaryPointList& yhi);
	
    /*! 
     * Destructor for the SAAMGPCG solver.
     */
    ~Solver();

    /** Actual computation of self field.
     * \param repetitions the number of calls to the solver.
     * \param tol the convergence criterion for the PCG solver.
     * \param maxit the maximal number of iterations performed by the PCG solver.  */
    void Compute();

    /** Performs setup of data distribution and problem (system matrix, RHS, LHS).  */
    void SetupProblem(const Array<MultiFab*>& rhs, const Array<MultiFab*>& soln); 

    /*!
     * Copy the solution of the linear system of equations (done by Trilinos)
     * back to BoxLib. This function is called in Electrostatic::solve_for_phi.
     * If OPAL is compiled with DBG_SCALARFIELD enabled, it also prints the potential
     * on all levels. At the moment just single-core supported.
     * @param soln will be filled with the solution of Trilinos.
     */
    void CopySolution(const Array<MultiFab*>& soln);

    /*!
     * Number of iterations
     * @returns the number of iterations used for the finding the solution.
     */
    inline int getNumIters();
    

private:
    /*! Solver (Belos BlockCG or RCG)
     * if nBlocks_m OR nRecycleBlocks_m == 0 BlockCG is used, else RCG
     */
    int nBlocks_m;                      ///< Maximum number of blocks in Krylov space
    int nRecycleBlocks_m;               ///< Number of vectors in recycle space
    unsigned int nLhs_m;                ///< Number of left-hand side solutions for extrapolating new guess.

    Epetra_MpiComm comm_m;              ///< Communicator for the solver.

    bool isReusingHierarchy_m;          ///< Flag specifying if hierarchy is reused
    bool isReusingPreconditioner_m;     ///< Flag specifying if whole preconditioner is reused

    std::string	 itsolver_m;            ///< Iterativer solver type


    Real hr[3];
    Real prob_lo[3];

    BoundaryPointList xlo_m;
    BoundaryPointList xhi_m;
    BoundaryPointList ylo_m;
    BoundaryPointList yhi_m;
    
    Box domain_m;                       ///< Dimension of the coarsest level --> problem dimension
    int nGridPoints_m[BL_SPACEDIM];     ///< Number of grid points in each direction.
    Array<std::unique_ptr<iMultiFab>>  idxMap_m;         ///< Mapping between Trilinos and BoxLib

    /// preconditioner object
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;

    /// parameter list used for the preconditioner (filled in SetupMLList)
    Teuchos::ParameterList MLList_m;

    /// system matrix
    Teuchos::RCP<Epetra_CrsMatrix> A;

    /// right hand side of our problem
    Teuchos::RCP<Epetra_Vector> RHS;

    /// left hand side
    Teuchos::RCP<Epetra_Vector> LHS;

    Teuchos::RCP<Epetra_MultiVector> P;
    std::deque< Epetra_Vector > OldLHS;

    Teuchos::ParameterList belosList_m;
    Belos::LinearProblem<ST,MV,OP> problem;
    Teuchos::RCP< Belos::EpetraPrecOp > prec;//TODO just wrapper for MLPrec->LinearProblem
    Teuchos::RCP< Belos::StatusTestGenResNorm< ST, MV, OP > > convStatusTest;
    Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > solver;

    /// verbosity of the code
    bool verbose_m;

    // -------------------------------------------------------------------------------------------------------------
    // PRIVATE MEMBER FUNCTIONS
    // -------------------------------------------------------------------------------------------------------------

    /*!
     * Fill in the list used for the solver
     * @param maxit is the maximum number of allowed iterations
     * @param tol is the tolerance
     */
    void setupBelosList_m(int maxit, double tol);

    /*!
     * Computes a new initial guess of the left-hand side (i.e. the potential)
     * with the Aitken-Neville scheme using the previous solutions.
     */
    void extrapolateLHS();
    
    /*!
     * Compute the global index used in construcing the system matrix from the
     * index of the box.
     * @param coord is the index in thre dimensions coming from the MultiFab.
     */
    inline int coordBox2globalIdx_m(const IntVect& coord);
    
    /*!
     * Get the global index of the Multifab and add it to the Epetra::Epetra_Map.
     * @param gidx contains all global index of this process
     * @param rhs is the MultiFab from which we compute the global indices.
     */
    void fillGlobalIndex_m(std::vector<int>& gidx, const MultiFab& rhs);
    
    /*!
     * Fill the right-hand side vector with the densities on the grid coming
     * from the MultiFab
     * @param rhs is the right-hand side
     */
    void fillRHS_m(const MultiFab& rhs);
    
    /*!
     * Fill the system matrix A of the linear system of equations: Ax = b
     * @param map is the element mapping to the processors
     */
    void fillSysMatrix_m(const Epetra_Map& map);
    
    /*!
     * Compute the stencil values.
     * @param coord is the coordinate (i, j, k)
     * @param values contains the values of the center, west, east, north, south, front, back in this ordering.
     * @param indices contains the global indices where the values are inserted in the row.
     * @param nEntries is the number of non-zero entries in the matrix row (number of valid values).
     * @param ifab does the mapping.
     */
    void getStencil(const IntVect& coord, double* values, int* indices, int& nEntries, const IArrayBox& ifab);
    
protected:
    /// Setup the parameters for the SAAMG preconditioner.
    inline void setupMultiLevelList();

};

#endif
