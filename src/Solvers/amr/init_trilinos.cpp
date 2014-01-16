
#include "ml_include.h"

//the following code cannot be compiled without these ML Trilinos packages
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#include <TrilinosSolver.h>
#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Electrostatic.H>

#include "Epetra_MpiComm.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include <vector>
#include <string>
#include <iostream>

using namespace Teuchos;

Solver*
Electrostatic::init_trilinos(PArray<MultiFab>& rhs, PArray<MultiFab>& soln, const Real* dx)
{
    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    // default values

    int maxiters = 1000;
    int numBlocks = 1;
    int recycleBlocks = 0;
    int maxOldLHS = 1;
    double tol = 1e-8;
    bool verbose = false;

    // Default to constant interpolation here.
    std::string interpl = "constant";

    /// A block CG iteration for SPD linear problems."CG"
    /// A block GMRES iteration for non-Hermitian linear problems. "GMRES"
    std::string itsolver = "GMRES"; 

    int nlevs = rhs.size();

    BoundaryPointList& xlo = parent->getIntersectLoX();
    BoundaryPointList& xhi = parent->getIntersectHiX();
    BoundaryPointList& ylo = parent->getIntersectLoY();
    BoundaryPointList& yhi = parent->getIntersectHiY();

    const Real* prob_lo = Geometry::ProbLo();

    Solver* s = new Solver(Comm, itsolver, interpl, verbose,
                           tol, maxiters, numBlocks, recycleBlocks, maxOldLHS,
                           nlevs, rhs, soln, dx, prob_lo, xlo, xhi, ylo, yhi);

   return s;
}


#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO) */
