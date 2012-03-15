#ifdef HAVE_ML_SOLVER

#include "MGPoissonSolver.h"
#include "Physics/Physics.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include <algorithm>
#include <fstream>

#include "EpetraExt_RowMatrixOut.h"
#include <Epetra_Import.h>

MGPoissonSolver::MGPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl):
  layout_m(fl),
  mesh_m(mesh),
  LHS(0),
  Comm(MPI_COMM_WORLD)
{ 
  domain_m = layout_m->getDomain();
  e_dim_tag decomp[3]; 

  for(int i=0; i<3; i++) {
    decomp[i] = layout_m->getRequestedDistribution(i);
    hr_m[i] = mesh_m->get_meshSpacing(i);
    orig_nr_m[i] = domain_m[i].length();
  }

  if(decomp[0] == 1 || decomp[1] == 1) {
    //*gmsg << "WARNING: MG SOLVER CAN ONLY RUN PARALLEL IN Z DIRECTION!" << endl;
    //exit(-1);
    throw OpalException("MGPoissonSolver","MGPoissonSolver only works parallel in z-direction");
  }

  //TODO: if BC == elliptic: ADD GEOMETRY COMMAND
  //read in semi-major and semi-minor defining the elliptical beam pipe
  semi_major = 1.0;
  semi_minor = 1.0;
  RealVariable *ma = dynamic_cast<RealVariable*>(OPAL.find("SEMIMAJOR"));
  RealVariable *mi = dynamic_cast<RealVariable*>(OPAL.find("SEMIMINOR"));
    
  if (ma) {
    semi_major = ma->getReal();
    *gmsg << "SEMIMAJOR=" << semi_major << endl;
  } else
    *gmsg << "using std SEMIMAJOR=" << semi_major << endl;

  if (mi) {
    semi_minor = mi->getReal();
    *gmsg << "SEMIMINOR=" << semi_minor << endl;
  } else
    *gmsg << "using std SEMIMINOR=" << semi_minor << endl;

  //TEST: ANALYTICAL SOL ON ZYL!
  semi_major = 0.5;
  semi_minor = 0.5;
  
  if(2*semi_minor < orig_nr_m[1]*hr_m[1]  || 2*semi_major < orig_nr_m[0]*hr_m[0]) {
    //*gmsg << "FATAL: semi-minor and semi-major have to be > than current grid" << endl;
    throw OpalException("MGPoissonSolver","semi-minor and semi-major have to be > than current grid");
  }
  
  //TODO: UNIT CONVERT!
  //double factor = 1e-3 * Physics::c;// * getdT();
  //semi_major *= factor;
  //semi_minor *= factor;

  //ANALYTICAL TESTCASE:
  nr_m[0] = 41;
  nr_m[1] = 41;
  nr_m[2] = 40;

  /*
  //IFF! hr_m is HERE = (1.0,1.0,1.0)!!!
  //the new rectangular grid containing the ellipse
  //has the size:
  nr_m[0] = (int) 2*ceil(semi_major/hr_m[0]);
  nr_m[1] = (int) 2*ceil(semi_minor/hr_m[1]);
  nr_m[2] = orig_nr_m[2];
  */

  cout << "new domain with enclosing ellipse: " << nr_m << endl;
  
  bp = new EllipticalDomain(semi_major, semi_minor, nr_m, hr_m);

  Map = 0;
  A = 0;
  LHS = 0;
 
  //all timers used here
  FunctionTimer1_m = IpplTimings::getTimer("Create Stencil");
  FunctionTimer2_m = IpplTimings::getTimer("Prepare RHS");
  FunctionTimer3_m = IpplTimings::getTimer("ML + CG");
  FunctionTimer4_m = IpplTimings::getTimer("LHS to rho");
  FunctionTimer5_m = IpplTimings::getTimer("recalc Map");
}


MGPoissonSolver::MGPoissonSolver(PartBunch &beam):
  layout_m(&beam.getFieldLayout()),
  mesh_m(&beam.getMesh()),
  itsBunch_m(&beam),
  LHS(0),
  Comm(MPI_COMM_WORLD)
{
  domain_m = layout_m->getDomain();
  e_dim_tag decomp[3]; 

  for(int i=0; i<3; i++) {
    decomp[i] = layout_m->getRequestedDistribution(i);
    hr_m[i] = mesh_m->get_meshSpacing(i);
    orig_nr_m[i] = domain_m[i].length();
  }
  
  if(decomp[0] == 1 || decomp[1] == 1) 
    *gmsg << "WARNING: MG SOLVER CAN ONLY BE PARALLEL IN Z DIRECTION!!" << endl;
  
  //TEST: ANALYTICAL SOL ON ZYL!
  semi_major = 0.5;
  semi_minor = 0.5;

  nr_m[0] = (int) 2*ceil(semi_major/hr_m[0]);
  nr_m[1] = (int) 2*ceil(semi_minor/hr_m[1]);
  nr_m[2] = orig_nr_m[2];
  bp = new EllipticalDomain(semi_major, semi_minor, nr_m, hr_m);
  
  Map = 0;
  A = 0;
  LHS = 0;
 
  FunctionTimer1_m = IpplTimings::getTimer("Create Stencil");
  FunctionTimer2_m = IpplTimings::getTimer("Prepare RHS");
  FunctionTimer3_m = IpplTimings::getTimer("ML + CG");
  FunctionTimer4_m = IpplTimings::getTimer("LHS to rho");
  FunctionTimer5_m = IpplTimings::getTimer("recalc Map");
}

////////////////////////////////////////////////////////////////////////////
// destructor
MGPoissonSolver::~MGPoissonSolver()
{
  //cleanup
  delete LHS;
  delete A;
  delete Map;
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation

void MGPoissonSolver::computePotential(Field_t &rho, Vector_t hr)
{

  //IFF: MATRIX FREE: reduce time to sol, increase MEM
  //IFF: use matrix stencil in computation directly (no Epetra, define operators on IPPL GRID)
  
  //IFF: params (sane settings? read from input file?)
  int maxIterations = 100;
  double tol=1e-8;
 
  //ANALYITICAL TEST CASE: in domain [0,1]
  nr_m[0] = 41;
  nr_m[1] = 41;
  nr_m[2] = 40;
  hr[0] = 0.025;
  hr[1] = 0.025;
  hr[2] = 0.025;

  /*
  //the new rectangular grid containing the ellipse
  //has the size:
  nr_m[0] = (int) 2*ceil(semi_major/hr[0]);
  nr_m[1] = (int) 2*ceil(semi_minor/hr[1]);
  nr_m[2] = orig_nr_m[2];
  */
  cout << "new domain with enclosing ellipse: " << nr_m << endl;

  //TODO: recreate mesh (elliptic stuff)
  //if mesh changes, ellipse may change!
  //recompute intersection (hr may have changed)
  cout << "recomputing boundary points" << endl;
  //IFF: may be obsolete with geometry command
  bp->setNr(nr_m);
  bp->Compute(hr);
  cout << "finished recomputing boundary points" << endl;

  IpplTimings::startTimer(FunctionTimer5_m);
  NDIndex<3> localidx =  layout_m->getLocalNDIndex();
  cout << "CPU " << Comm.MyPID() << " " << localidx << endl;
  if(Map != 0)
    delete Map;
  Map = IPPLToMap3D(localidx);
  
  //we also have to redistribute LHS
  if(LHS == 0) {
    LHS = new Epetra_Vector(*Map);
    LHS->Random();
  } else {
    Epetra_Vector *tmp = new Epetra_Vector(*Map);
    Epetra_Import importer(*Map, LHS->Map());
    tmp->Import(*LHS, importer, Add);
    delete LHS;
    LHS = tmp;
  }
  IpplTimings::stopTimer(FunctionTimer5_m);
  cout << "created NEW map" << endl;

  IpplTimings::startTimer(FunctionTimer2_m);

  Epetra_Vector RHS(*Map);
  //Epetra_Vector RHS(A->Map());
  RHS.PutScalar(0.0);
  register int idx = 0;
  register double fact = 1.0;
  int nxy = bp->getNumXY(0);

  ///TEST: ANALYTICAL SOLUTION TEST
//#ifdef ANALY_TEST
  double r = 0.0;
  for(int x = 0; x < nr_m[0]; x++) {
    for(int y = 0; y < nr_m[1]; y++) {
      for(int z = 0; z < nr_m[2]; z++) {
    //for(int xy = 0; xy < bp->getNumXY(0); xy++) {
      //for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

        fact = 1.0;
        if(z == 0 || z == nr_m[2]-1) {fact = 0.5;} 
        if(bp->isInside(x,y,z)) {

          //RHS.Values()[idx++] = idx;

          double xx = (x-nr_m[0]/2)*hr[0];
          double yy = (y-nr_m[1]/2)*hr[1];
          double zz = (z-nr_m[2]/2)*hr[2];
          r = sqrt(xx*xx+yy*yy);
          if(r >= 0.5)
            RHS.Values()[idx++] = 0.0;
          else
            RHS.Values()[idx++] = fact*(96*(1-4*r)+(1-2*r)*(1-2*r)*(6*r+1)*Physics::pi*Physics::pi)*(1-2*r)*sin(Physics::pi*(zz-0.5));

        }
        }
      }
    }
  //}
//#endif

  /*
  //IFF: get charge densities from ippl field -> Epetra vector (RHS)
  idx = 0;
  RHS.PutScalar(0.0);
  for(int x = 0; x < nr_m[0]; x++) {
    for(int y = 0; y < nr_m[1]; y++) {
      for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {


        if(bp->isInside(x,y,z)) {
          RHS.Values()[idx++] = rho[x][y][z].get();
        }

      }
    }
  }*/

  IpplTimings::stopTimer(FunctionTimer2_m);

  IpplTimings::startTimer(FunctionTimer1_m);

  if(A != 0)
    delete A;
  //A = CreateCrsMatrix("Cross3D", Map, GaleriList);
  //A = CreateCrsMatrix("Laplace3D", Map, GaleriList);
  //A = Stencil3DLongitudinalNeumann(hr);
  A = Stencil3DGeometry(hr, RHS);
  IpplTimings::stopTimer(FunctionTimer1_m);

  //debug out
#ifdef DBG_STENCIL 
  EpetraExt::RowMatrixToMatlabFile("A.dat", *A);
#endif

  //use old LHS solution as start vector
  //Epetra_LinearProblem Problem(A, LHS, &RHS);
  AztecOO solver(A, LHS, &RHS);

  /////////////////////////////////////////
  //PROBLEM SETUP FINISHED -- PRECOND BELOW
  /////////////////////////////////////////
  
  IpplTimings::startTimer(FunctionTimer3_m);
  // set defaults
  ML_Epetra::SetDefaults("SA", MLList);
  MLList.set("max levels", 5);
  MLList.set("increasing or decreasing", "decreasing");
  MLList.set("prec type", "MGW");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: type","Chebyshev");
  MLList.set("smoother: sweeps",3);
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: type", "Amesos-KLU");

  // try to optimize mem for xt3
  //MLList.set("low memory usage", true);
  //maybe also helps
  //MLList.set("coarse: max size", 1024);

  // create the preconditioner object and compute hierarchy
  // true -> create the multilevel hirarchy
  // IFF: maybe we can use one precond obj and then use .setMatrix() and .Construct() (or something
  // similar)
  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
  
  ////////////////////////////////////////
  //PRECOND SETUP FINISHED -- SOLVER BELOW
  ////////////////////////////////////////

  //CG SOLVER 
  solver.SetPrecOperator(MLPrec);

  //solver.SetAllAztecOptions(options);
  //solver.SetAllAztecParams(params);
  //solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
 
  //IFF: need to decide when to use CG and when non-symmetric solver
  //if(semi_major == semi_minor)
  //  solver.SetAztecOption(AZ_solver, AZ_cg);
  //else
    solver.SetAztecOption(AZ_solver, AZ_bicgstab);

  solver.SetAztecOption(AZ_conv, AZ_rhs);
  solver.SetAztecParam(AZ_kspace, maxIterations);
  solver.SetAztecOption(AZ_output, AZ_all);
  //solver.SetAztecOption(AZ_output, AZ_none);
  solver.Iterate(maxIterations, tol);

  if( Comm.MyPID()==0 ) 
    cout << "\t\t||b-Ax||_2 = " << solver.TrueResidual() << endl;
  IpplTimings::stopTimer(FunctionTimer3_m);

  //now transfer solution back to grid
  //only original grid back to rho

  IpplTimings::startTimer(FunctionTimer4_m);

  idx = 0;
 
  ///TEST: DUMP SOLUTION
//#ifdef ANALY_TEST
  *gmsg << "*** START DUMPING SOLUTION ***" << endl;
  ostringstream oss;
  MPI_File file;
  MPI_Status status;
  MPI_Info fileinfo;
  MPI_File_open(MPI_COMM_WORLD, "rho_solution", MPI_MODE_WRONLY | MPI_MODE_CREATE, fileinfo, &file);
  for(int x=localidx[0].first(); x<=localidx[0].last(); x++) {
    for(int y=localidx[1].first(); y<=localidx[1].last(); y++) {
      for(int z=localidx[2].first(); z<=localidx[2].last(); z++) {
         if(bp->isInside(x,y,z))
           oss << x+1 << " " << y+1 << " " << z+1 << " " <<  LHS->Values()[idx++]/*(hr[0]*hr[1]*hr[2])*/ << endl;
         else
           oss << x+1 << " " << y+1 << " " << z+1 << " " <<  "0.0" << endl;
      }
    }
  }

  char* tmpch = new char [oss.str().size()+1];
  strcpy(tmpch, oss.str().c_str());
  MPI_File_write_shared(file, tmpch, oss.str().size(), MPI_CHAR, &status);
  //MPI_File_write_shared(file, (char*)oss.str().c_str(), oss.str().size(), MPI_CHAR, &status);
  MPI_File_close(&file);
  delete tmpch;
  *gmsg << "*** FINISHED DUMPING SOLUTION ***" << endl;
//#endif

  /*
  //IFF: put solution to ippl grid
  idx = 0;
  for(int x = 0; x < nr_m[0]; x++) {
    for(int y = 0; y < nr_m[1]; y++) {
      for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

      NDIndex<3> l(Index(x,x),Index(y,y),Index(z,z));
        if(bp->isInside(x,y,z)) 
          rho.localElement(l) = LHS->Values()[idx++];
        else
          rho.localElement(l) = 0.0;

      }
    }
  }
  */

  IpplTimings::stopTimer(FunctionTimer4_m);

  //cleanup
  delete MLPrec;

}

/*inline bool MGPoissonSolver::isOnFrontSurface(const int idx) 
{

  int nxy = nr_m[0]*nr_m[1];
  int ixy = idx % nxy;
  return ((idx - ixy) / nxy) == 0;

}

inline bool MGPoissonSolver::isOnBackSurface(const int idx) 
{

  int nxy = nr_m[0]*nr_m[1];
  int ixy = idx % nxy;
  return ((idx - ixy) / nxy) == nr_m[2]-1;

}*/

//NOTE: the map is only PARALLEL in Z direction!
//IFF FIXME: think i should move this to the appropriate BoundaryPoints class
inline Epetra_Map* MGPoissonSolver::IPPLToMap3D(NDIndex<3> localidx)
{
 
  //Version for Elliptic Boundary Points
  int idx = 0;
  int xy = bp->getNumXY(0);
  cout << "my num xy = " << xy << endl;
  int NumMyElements = xy*localidx[2].length();
  cout << "my num eles = " << NumMyElements << endl;
  vector<int> MyGlobalElements(NumMyElements);
 
  for(int pl = 0; pl < xy; pl++) {
    //for (int z = 0; z < nr_m[2]; z++)
    for (int z = localidx[2].first(); z <= localidx[2].last(); z++)
      MyGlobalElements[idx++] = pl + z * xy;
  }

  /*
  // Version ? 
  for(int x=0; x < nr_m[0]; x++) {
    for(int y=0; y < nr_m[1]; y++) {
      for(int z=0; z < nr_m[2]; z++) {
  
        if(bp->isInside(x,y,z))
          MyGlobalElements[idx++] = 

      }
    }
  }*/

  return(new Epetra_Map (-1, NumMyElements, &MyGlobalElements[0], 0, Comm));

}

///////////////////////////////////// STENCIL CREATION //////////////////////////////////////////

//at the moment this stencil has neumann BC in z-direction and dirichlet on the ellipse surface in x,y-direction
//this code is experimental/not-optimized/not-final at the moment
//IFF: since we are using lists for intersection with X and Y direction gridlines we can reuse this stencil 
//even with complex geometries. This method only depends on the two lists.
inline Epetra_CrsMatrix* MGPoissonSolver::Stencil3DGeometry(Vector_t hr, Epetra_Vector& RHS)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  7);

  //this stencil will be unsymmetric when we have an ellipse
  //when we have a circle we get a symmetric discr. matrix
  //the first and the last xy-plane (z=0, z=max) have a differenet discretization! (neuman in z dir)
  
  //ALGORITHM IDEA:
  //in this method we go over all intersections with gridlines in y-dir to number all interior points
  //then we calculate the shortely-weller approximation for boundary points

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();
  int i=0;

  vector<double> Values(6);
  vector<int> Indices(6);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {

    int NumEntries = 0;

    //if(bp->isInside(x,y,z)) {

    //create stencil
    double WV, EV, SV, NV, FV, BV, CV, scaleFactor;
    //bp->getBoundaryStencil(x, y, z, WV, EV, SV, NV, FV, BV, CV);
    bp->getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, FV, BV, CV, scaleFactor);

    //scale RHS:
    RHS[i] *= scaleFactor;

    double W, E, S, N, F, B;
    //bp->getNeighbours(x, y, z, W, E, S, N, F, B);
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

    // put the off-diagonal entries
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);

    // Put in the diagonal entry
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);

  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);

}

void MGPoissonSolver::GetNeighbours(const int idx, int& W, int& E, int& S, int& N, int& F, int& B, int& numOutOfDomain)
{ 

  int ixy, iz, iy, ix; 
  int nx = nr_m[0];
  int ny = nr_m[1];
  int nz = nr_m[2];
  int nxy = nx*ny;
  ixy = idx % nxy;

  ix = ixy % nx;
  iy = (ixy - ix) / nx;
  iz = (idx - ixy) / nxy;
  numOutOfDomain = 0;

  if (iz == 0)
    {F = -1; numOutOfDomain++;}
  else              
    F = idx - nxy; 
  if (iz == nz - 1) 
    {B = -1; numOutOfDomain++;}
  else              
    B = idx + nxy;

  if (ix == 0)      
    {W = -1; numOutOfDomain++;}
  else              
    W = ixy - 1;
  if (ix == nx - 1) 
    {E = -1; numOutOfDomain++;}
  else              
    E = ixy + 1;

  if (iy == 0)      
    {S = -1; numOutOfDomain++;}
  else              
    S = ixy - nx;
  if (iy == ny - 1) 
    {N = -1; numOutOfDomain++;}
  else              
    N = ixy + nx;

  int cor = iz * nxy;
  switch(W) {
    case -1: break;
    default: W += cor;
  }
  switch(E) {
    case -1: break;
    default: E += cor;
  }
  switch(S) {
    case -1: break;
    default: S += cor;
  }
  switch(N) {
    case -1: break;
    default: N += cor;
  }

}

//at the moment this stencil has neumann everywhere except on the z=0 plane (Dirichlet)
//this code is experimental/not-optimized/not-final at the moment
inline Epetra_CrsMatrix* MGPoissonSolver::Stencil3DOneSidedDirichlet(Vector_t hr)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  7);
    
  Vector_t hr2 = hr*hr;
  double a = 2*(1/hr2[0]+1/hr2[1]+1/hr2[2]);
  double b = -1/hr2[0];
  double d = -1/hr2[1];
  double f = -1/hr2[2]; 

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int W, E, S, N, F, B, numout;
  vector<double> Values(6);
  vector<int> Indices(6);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    GetNeighbours(MyGlobalElements[i], W, E, S, N, F, B, numout);

    int NumEntries = 0;
    double diag = a;

    if(F != -1) {
    
    switch(numout) {

       case 3: {
        if(W == -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.25*b;
        }

        if(E == -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.25*b;
        }

        if(N == -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.25*d;
        }

        if(S == -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.25*d;
        }

        if(B == -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = 0.25*f;
        }

        if(F == -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = 0.25*f;
        }

        diag *= 0.125;
        
        if(NumEntries != 3)
          cout << "ERROR: while calc corners in stencil" << endl;

        break;

      } 
    case 2: { 

        if(W != -1 && E != -1) { 
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.25*b;
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.25*b;
        }
        if(E == -1 && W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
        }
        if(W == -1 && E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }

        
        if(N != -1 && S != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.25*d;
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.25*d;
        }
        if(S == -1 && N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }
        if(N == -1 && S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }

        
        if(B != -1 && F != -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = 0.25*f;
          Indices[NumEntries] = F;
          Values[NumEntries++] = 0.25*f;
        }
        if(B == -1 && F != -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = 0.5*f;
        }
        if(F == -1 && B != -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = 0.5*f;
        }

        diag *= 0.25;

        if(NumEntries != 4)
          cout << "ERROR: while calc edge stencil elements" << endl;

        break;
      } 
    case 1: {

        
        if(W != -1 && E != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        } 
        if(E == -1 && W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = b;
        }
        if(W == -1 && E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = b;
        }

        
        if(N != -1 && S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }
        if(S == -1 && N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = d;
        }
        if(N == -1 && S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = d;
        }

        
        if(B != -1 && F != -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = 0.5*f;
          Indices[NumEntries] = B;
          Values[NumEntries++] = 0.5*f;
        }
        if(B == -1 && F != -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        }
        if(F == -1 && B != -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        diag *= 0.5;

        if(NumEntries != 5)
          cout << "ERROR: calc surface elements of stencil" << endl;

        break;
      }

    case 0: { 
        if(W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = b;
        } if(E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = b;
        } if(S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = d;
        } if(N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = d;
        } if(F != -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        } if(B != -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        break;
      }
    default: {
        cout << "ERROR CREATING THE STENCIL: points out of domain (" << numout << ") is >3 OR <0" << endl;
      }
    }
    } else {
      //F == -1 --> this is our dirichlet boundary surface
      switch(numout) {
        case 3: {
          
        if(W == -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }

        if(E == -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
        }

        if(N == -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }

        if(S == -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }

        Indices[NumEntries] = B;
        Values[NumEntries++] = 0.25*f;
        
        if(NumEntries != 3)
          cout << "ERROR: calc corner (below) elements of stencil" << endl;

        diag *= 0.25;
        break;

        }
        case 2:{

        if(W != -1 && E != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }
        if(W == -1 && E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = b;
        }

        if(E == -1 && W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = b;
        }

        if(N != -1 && S != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }
        if(N == -1 && S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = d;
        }

        if(S == -1 && N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = d;
        }

        Indices[NumEntries] = B;
        Values[NumEntries++] = 0.25*f;

        if(NumEntries != 4)
          cout << "ERROR: calc edge (below) elements of stencil" << endl;

        diag *= 0.25;
        break;

        }
        case 1:{

          Indices[NumEntries] = E;
          Values[NumEntries++] = b;

          Indices[NumEntries] = W;
          Values[NumEntries++] = b;

          Indices[NumEntries] = S;
          Values[NumEntries++] = d;

          Indices[NumEntries] = N;
          Values[NumEntries++] = d;

          Indices[NumEntries] = B;
          Values[NumEntries++] = f;

          if(NumEntries != 5)
            cout << "ERROR: calc surface (below) elements of stencil: " << NumEntries << endl;

          break;
        
        }
        default: {
          cout << "error in dirichlet surface!" << endl;
        }
      }
    }
    // put the off-diagonal entries
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
	
    // Put in the diagonal entry
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, 
                               &diag, MyGlobalElements + i);

  }
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

//stencil for longitudinal neumann and transversal dirichlet BC
//this code is experimental/not-optimized/not-final at the moment
inline Epetra_CrsMatrix* MGPoissonSolver::Stencil3DLongitudinalNeumann(Vector_t hr)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  7);
    
  Vector_t hr2 = hr*hr;
  double a = 2*(1/hr2[0]+1/hr2[1]+1/hr2[2]);
  double b = -1/hr2[0];
  double d = -1/hr2[1];
  double f = -1/hr2[2]; 

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int W, E, S, N, F, B, numout;
  vector<double> Values(6);
  vector<int> Indices(6);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    GetNeighbours(MyGlobalElements[i], W, E, S, N, F, B, numout);

    int NumEntries = 0;
    double diag = a;
    //numout = 0;
    
    //transversal direction: if below == -1 || above == -1 ===> NEUMANN
    //else dirichlet:
    if(F != -1 && B != -1)
      numout = 0;
    
    //only on left longitudinal end open bc -- else dirichlet
    //if(F != -1)
    //  numout = 0;

    switch(numout) {

       case 3: {
        if(W == -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }

        if(E == -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
        }

        if(N == -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }

        if(S == -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }

        if(B == -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        }

        if(F == -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        diag *= 0.5;
        
        if(NumEntries != 3)
          cout << "ERROR: while calc corners in stencil" << endl;

        break;

      } 
    case 2: { //edges

        if(E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }
        if(W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
        }

        if(S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }
        if(N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }

        if(B == -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        }
        if(F == -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        diag *= 0.5;

        if(NumEntries != 4)
          cout << "ERROR: while calc edge stencil elements" << endl;

        break;
      } 
    case 1: {

        if(E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = 0.5*b;
        }
        if(W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = 0.5*b;
        }

        if(S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = 0.5*d;
        }
        if(N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = 0.5*d;
        }

        if(B == -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        }
        if(F == -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        diag *= 0.5;

        if(NumEntries != 5)
          cout << "ERROR: calc surface elements of stencil" << endl;

        break;
      }

    case 0: {  //interior points (& dirichlet)
        if(W != -1) {
          Indices[NumEntries] = W;
          Values[NumEntries++] = b;
        } if(E != -1) {
          Indices[NumEntries] = E;
          Values[NumEntries++] = b;
        } if(S != -1) {
          Indices[NumEntries] = S;
          Values[NumEntries++] = d;
        } if(N != -1) {
          Indices[NumEntries] = N;
          Values[NumEntries++] = d;
        } if(F != -1) {
          Indices[NumEntries] = F;
          Values[NumEntries++] = f;
        } if(B != -1) {
          Indices[NumEntries] = B;
          Values[NumEntries++] = f;
        }

        break;
      }
    default: {
        cout << "ERROR CREATING THE STENCIL: points out of domain (" << numout << ") is >3 OR <0" << endl;
      }
    }
    // put the off-diagonal entries
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
	
    // Put in the diagonal entry
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, 
                               &diag, MyGlobalElements + i);

  }
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}


Inform &MGPoissonSolver::print(Inform &os) const
{
  os << "* *************** M G P o i s s o n S o l v e r ************************************ " << endl;
  os << "* h " << hr_m << '\n';
  os << "* ********************************************************************************** " << endl;
}

#endif /* HAVE_ML_SOLVER */
