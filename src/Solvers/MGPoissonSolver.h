////////////////////////////////////////////////////////////////////////////
// This class contains methods for solving Poisson's equation for the 
// space charge portion of the calculation.
////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ML_SOLVER

#ifndef MG_POISSON_SOLVER_H_
#define MG_POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Ippl.h"
//class PartBunch;
class MGPoissonSolver;
#include "Algorithms/PartBunch.h"
#include "PoissonSolver.h"
#include "IrregularDomain.h"
#include "EllipticalDomain.h"
#include "RectangularDomain.h"
#include "../Structure/BoundaryGeometry.h"
//////////////////////////////////////////////////////////////
#include "ml_include.h"
 
// the following code cannot be compiled without these Trilinos
// packages. Note that Galeri is required in the examples only (to
// generate the linear system), not by the ML library
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_Utils.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include "ml_MatrixFreePreconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

using namespace Teuchos;
using namespace Galeri;
using namespace ML_Epetra;
//////////////////////////////////////////////////////////////

typedef UniformCartesian<3,double> Mesh_t;
typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

enum {
  CG,
  BiCGStab,
  GMRES,
  NONE
};

/**
 * \class MGPoissonSolver
 * \brief a Multigrid solver for space charges
 * \see FFTPoissonSolver
 * \warning this solver is in an EXPERIMENTAL STAGE. For reliable simulations use the FFTPoissonSolver
 *
 */
class MGPoissonSolver : public PoissonSolver
{
public:
  /// constructor
  MGPoissonSolver(PartBunch &bunch);

  /// constructor
  MGPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::vector<BoundaryGeometry*> geometries, std::string itsolver, std::string interpl, double tol, int maxiters);
  
  /// destructor 
  ~MGPoissonSolver();
    
  /// given a charge-density field rho and a set of mesh spacings hr, compute the scalar potential in 'open space'
  /// \param rho (inout) scalar field of the potential
  /// \param hr mesh spacings in each direction
  void computePotential(Field_t &rho, Vector_t hr);
  /// this function is unused at the moment and should NOT be called!
  void computePotential(Field_t &rho, Vector_t hr, double zshift) {};

  void setGeometry(std::vector<BoundaryGeometry*> geometries);

  double getXRangeMin() { return bp->getXRangeMin(); }
  double getXRangeMax() { return bp->getXRangeMax(); }
  double getYRangeMin() { return bp->getYRangeMin(); }
  double getYRangeMax() { return bp->getYRangeMax(); }

  Inform &print(Inform &os) const;

private:

  /// holding the currently active geometry
  //TODO: we need to update this an maybe change attached
  //solver!
  BoundaryGeometry* currentGeometry;

  //TODO: how set value clean?
  bool hasGeometryChanged;
  bool hasDomainDecompositionChanged;

  //TODO: remove!
  /// semi-major of the elliptical beam pipe 
  double semi_major_m;
  /// semi-minor of the elliptical beam pipe 
  double semi_minor_m;

  /// tolerance for the iterative solver
  double tol_m;
  /// maximal number of iterations for the iterative solver
  int maxiters_m;
  /// iterative solver we are applying: CG, BiCGStab or GMRES
  int itsolver_m;

  /// structure that holds boundary points
  IrregularDomain *bp;

  /// left hand side of the linear system of equations we solve
  Epetra_Vector *LHS;
  /// matrix used in the linear system of equations
  Epetra_CrsMatrix *A;
  /// an Epetra_Map holding the processor distribution of data
  Epetra_Map* Map;
  /// communicator used by Trilinos
  Epetra_MpiComm Comm;
  /// ML preconditioner object
  ML_Epetra::MultiLevelPreconditioner* MLPrec;

  /// parameter list for the ML solver
  Teuchos::ParameterList MLList;
  /// parameter list for Galeri usage
  Teuchos::ParameterList GaleriList;

  // mesh and layout objects for rho_m
  Mesh_t *mesh_m;
  FieldLayout_t *layout_m;

  /// store geometries
  std::vector<BoundaryGeometry*> geometries_m;
    
  // domains for the various fields
  NDIndex<3> domain_m;
     
  /// mesh spacings in each direction
  Vector_t hr_m;
  /// global number of nodes in each direction
  Vektor<int,3> nr_m;
  Vektor<int,3> orig_nr_m;
  
  /// PartBunch object
  PartBunch *itsBunch_m;

  // timers
  IpplTimings::TimerRef FunctionTimer1_m;
  IpplTimings::TimerRef FunctionTimer2_m;
  IpplTimings::TimerRef FunctionTimer3_m;
  IpplTimings::TimerRef FunctionTimer4_m;
  IpplTimings::TimerRef FunctionTimer5_m;
  IpplTimings::TimerRef FunctionTimer6_m;
 
  /// converts IPPL grid to a 3D Epetra_Map
  /// \param localidx local IPPL grid node indices
  /// \return Epetra_Map corresponding to the parallel IPPL grid
  inline Epetra_Map* IPPLToMap3D(NDIndex<3> localidx);
  
  /// returns a discretized stencil that has Neumann BC everywhere except in the z=0 plane Dirichlet BC
  /// \param hr gridspacings in each direction
  /// \return Epetra_CrsMatrix holding the stencil
  inline Epetra_CrsMatrix* Stencil3DOneSidedDirichlet(Vector_t hr);
  /// returns a discretized stencil that has Neumann BC in longitudinal direction and Dirichlet BC in transversal direction
  /// \param hr gridspacings in each direction
  /// \return Epetra_CrsMatrix holding the stencil
  inline Epetra_CrsMatrix* Stencil3DLongitudinalNeumann(Vector_t hr);
  /// returns a discretized stencil that has Neumann BC in z direction and Dirichlet BC on the surface of a specified geometry
  /// \param hr gridspacings in each direction
  /// \return Epetra_CrsMatrix holding the stencil
  inline Epetra_CrsMatrix* Stencil3DGeometry(Vector_t hr, Epetra_CrsMatrix*& Matrix, Epetra_Vector& RHS);

  /// helper function to determine position
  /// \param idx index of the current element in the matrix
  /// \param W index of the element in the west of idx (x-1)
  /// \param E index of the element in the east of idx (x+1)
  /// \param S index of the element in the south of idx (y-1)
  /// \param N index of the element in the north of idx (y+1)
  /// \param F index of the element in front of idx (z-1)
  /// \param B index of the element in the back of idx (z+1)
  /// \param numOutOfDomain total number of elements out of bounds
  void GetNeighbours(const int idx, int& W, int& E, int& S, int& N, int& F, int& B, int& numOutOfDomain);

  /*
  /// helper function to check if a given grid node is on the front surface (z=0)
  inline bool isOnFrontSurface(const int idx);
  /// helper function to check if a given grid node is on the back surface (z=max_z)
  inline bool isOnBackSurface(const int idx);
  */

};

inline Inform &operator<<(Inform &os, const MGPoissonSolver &fs)
{
  return fs.print(os);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) */

#endif /* #ifndef MG_POISSON_SOLVER_H_ */

#endif /* #ifdef HAVE_ML_SOLVER */

/***************************************************************************
 * $RCSfile: MGPoissonSolver.hh,v $   $Author: ineichen $
 * $Revision: 1.1.1.1 $   $Date: 2008/19/03 11:21:48 $
 ***************************************************************************/
