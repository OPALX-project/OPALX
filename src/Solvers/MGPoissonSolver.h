////////////////////////////////////////////////////////////////////////////
// This class contains methods for solving Poisson's equation for the
// space charge portion of the calculation.
////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ML_SOLVER

#ifndef MG_POISSON_SOLVER_H_
#define MG_POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Structure/BoundaryGeometry.h"
#include "PoissonSolver.h"
#include "IrregularDomain.h"
#include "EllipticDomain.h"
#include "BoxCornerDomain.h"
#include "RectangularDomain.h"
//////////////////////////////////////////////////////////////
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"
#include <Epetra_Import.h>

#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_ParameterList.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using namespace ML_Epetra;
using namespace Isorropia;
//////////////////////////////////////////////////////////////

typedef UniformCartesian<3, double> Mesh_t;
typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

enum {
    STD_PREC,
    REUSE_PREC,
    REUSE_HIERARCHY
};

/**
 * \class MGPoissonSolver
 * \brief A smoothed aggregation based AMG preconditioned iterative solver for space charge
 * \see FFTPoissonSolver
 * \warning This solver is in an EXPERIMENTAL STAGE. For reliable simulations use the FFTPoissonSolver
 *
 */
class PartBunch;
class BoundaryGeometry;

class MGPoissonSolver : public PoissonSolver {

public:
    MGPoissonSolver(PartBunch &beam,Mesh_t *mesh, FieldLayout_t *fl, std::vector<BoundaryGeometry *> geometries, std::string itsolver, std::string interpl, double tol, int maxiters, std::string precmode);
    MGPoissonSolver(PartBunch &bunch);
    ~MGPoissonSolver();

    /// given a charge-density field rho and a set of mesh spacings hr, compute
    /// the scalar potential in 'open space'
    /// \param rho (inout) scalar field of the potential
    /// \param hr mesh spacings in each direction
    void computePotential(Field_t &rho, Vector_t hr);
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    /// set a geometry
    void setGeometry(std::vector<BoundaryGeometry *> geometries);

    /// force Solver to recompute Epetra_Map
    void recomputeMap() { hasParallelDecompositionChanged_m = true; }

    double getXRangeMin() { return bp->getXRangeMin(); }
    double getXRangeMax() { return bp->getXRangeMax(); }
    double getYRangeMin() { return bp->getYRangeMin(); }
    double getYRangeMax() { return bp->getYRangeMax(); }

    /// useful load balance information
    void printLoadBalanceStats();

    Inform &print(Inform &os) const;

private:

    //TODO: we need to update this an maybe change attached
    //solver!
    /// holding the currently active geometry
    BoundaryGeometry *currentGeometry;
    /// container for multiple geometries
    std::vector<BoundaryGeometry *> geometries_m;


    /// flag notifying us that the geometry (discretization) has changed
    bool hasGeometryChanged_m;
    /// flag is set when OPAL changed decomposition of mesh
    bool hasParallelDecompositionChanged_m;
    /// flag specifying if problem is redistributed with RCB
    bool useRCB_m;
    /// flag specifying if we are verbose
    bool verbose_m;

    /// tolerance for the iterative solver
    double tol_m;
    /// maximal number of iterations for the iterative solver
    int maxiters_m;
    /// iterative solver we are applying: CG, BiCGStab or GMRES
    int itsolver_m;
    /// preconditioner mode
    int precmode_m;
    /// number of iterations in the solve of the previous time step
    int numIter_m;
    /// percentage the iteration count can increase before recomputing the preconditioner
    int tolerableIterationsCount_m;
    /// force the solver to recompute the preconditioner
    bool forcePreconditionerRecomputation_m;
    /// maximum number of blocks in Krylov space
    int numBlocks_m;
    /// number of vectors in recycle space
    int recycleBlocks_m;

    /// structure that holds boundary points
    IrregularDomain *bp;

    /// right hand side of our problem
    RCP<Epetra_Vector> RHS;
    /// left hand side of the linear system of equations we solve
    RCP<Epetra_Vector> LHS;
    /// matrix used in the linear system of equations
    RCP<Epetra_CrsMatrix> A;
    /// ML preconditioner object
    RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;
    /// Epetra_Map holding the processor distribution of data
    Epetra_Map *Map;
    /// communicator used by Trilinos
    Epetra_MpiComm Comm;

    /// last N LHS's for extrapolating the new LHS as starting vector
    uint nLHS;
    RCP<Epetra_MultiVector> P;
    std::deque< Epetra_Vector > OldLHS;

    /// Solver (Belos RCG)
    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST, MV, OP> OPT;
    typedef Belos::MultiVecTraits<ST, MV>    MVT;
    Belos::LinearProblem<double, MV, OP> problem;
    RCP< Belos::EpetraPrecOp > prec;
    RCP< Belos::StatusTestGenResNorm< ST, MV, OP > > convStatusTest;
    RCP< Belos::SolverManager<double, MV, OP> > solver;

    /// parameter list for the ML solver
    Teuchos::ParameterList MLList_m;
    /// parameter list for the iterative solver (Belos)
    Teuchos::ParameterList belosList;

    /// PartBunch object
    PartBunch *itsBunch_m;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    // domains for the various fields
    NDIndex<3> domain_m;

    /// mesh spacings in each direction
    Vector_t hr_m;
    /// current number of mesh points in each direction
    Vektor<int, 3> nr_m;
    /// global number of mesh points in each direction
    Vektor<int, 3> orig_nr_m;

    // timers
    IpplTimings::TimerRef FunctionTimer1_m;
    IpplTimings::TimerRef FunctionTimer2_m;
    IpplTimings::TimerRef FunctionTimer3_m;
    IpplTimings::TimerRef FunctionTimer4_m;
    IpplTimings::TimerRef FunctionTimer5_m;
    IpplTimings::TimerRef FunctionTimer6_m;

    /// recomputes the Epetra_Map
    void computeMap();

    /// redistributes Map with RCB
    /// \param localidx local IPPL grid node indices
    void redistributeWithRCB(NDIndex<3> localidx);

    /// converts IPPL grid to a 3D Epetra_Map
    /// \param localidx local IPPL grid node indices
    void IPPLToMap3D(NDIndex<3> localidx);

    /** returns a discretized stencil that has Neumann BC in z direction and
     * Dirichlet BC on the surface of a specified geometry
     * \param hr gridspacings in each direction
     * \param Epetra_CrsMatrix holding the stencil
     * \param RHS right hand side might be scaled
     */
    void ComputeStencil(Vector_t hr, Teuchos::RCP<Epetra_Vector> RHS);

protected:

    /// Setup the parameters for the Belos iterative solver.
    inline void SetupBelosList() {
        belosList.set("Maximum Iterations", maxiters_m);
        belosList.set("Convergence Tolerance", tol_m);

        if(numBlocks_m != 0 && recycleBlocks_m != 0){//only set if solver==RCGSolMgr
            belosList.set("Num Blocks", numBlocks_m);               // Maximum number of blocks in Krylov space
            belosList.set("Num Recycled Blocks", recycleBlocks_m); // Number of vectors in recycle space
        }
        if(verbose_m) {
            belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);
            belosList.set("Output Frequency", 1);
        } else
            belosList.set("Verbosity", Belos::Errors + Belos::Warnings);
    }

    /// Setup the parameters for the SAAMG preconditioner.
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
        MLList_m.set("smoother: type", "Chebyshev");
        MLList_m.set("smoother: sweeps", 3);
        MLList_m.set("smoother: pre or post", "both");

        // on the coarsest level we solve with  Tim Davis' implementation of
        // Gilbert-Peierl's left-looking sparse partial pivoting algorithm,
        // with Eisenstat & Liu's symmetric pruning. Gilbert's version appears
        // as \c [L,U,P]=lu(A) in MATLAB. It doesn't exploit dense matrix
        // kernels, but it is the only sparse LU factorization algorithm known to be
        // asymptotically optimal, in the sense that it takes time proportional to the
        // number of floating-point operations.
      //  MLList_m.set("coarse: type", "Amesos-KLU");

        //FIXME: CHEBY COARSE LEVEL SOLVER
        // SEE PAPER FOR EVALUATION KLU vs. Chebyshev
        MLList_m.set("coarse: sweeps", 20);
        MLList_m.set("coarse: type", "Chebyshev");

        // turn on all output
        if(verbose_m)
            MLList_m.set("ML output", 101);
        else
            MLList_m.set("ML output", 0);

        // try to optimize mem for xt3
        //MLList_m.set("low memory usage", true);

        // heuristic for max coarse size depending on number of processors
        int coarsest_size = std::max(Comm.NumProc() * 10, 1024);
        MLList_m.set("coarse: max size", coarsest_size);
    }

};

inline Inform &operator<<(Inform &os, const MGPoissonSolver &fs) {
    return fs.print(os);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[]) {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    puts("Please configure ML with:");
    puts("--enable-epetra");
    puts("--enable-teuchos");
    puts("--enable-aztecoo");

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO) */

#endif /* #ifndef MG_POISSON_SOLVER_H_ */

#endif /* #ifdef HAVE_ML_SOLVER */
