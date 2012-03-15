#ifdef HAVE_ML_SOLVER

#include "MGPoissonSolver.h"
#include "Physics/Physics.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
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
    nr_m[i] = domain_m[i].length();
  }

  if(decomp[0] == 1 || decomp[1] == 1) 
    *gmsg << "WARNING: MG SOLVER CAN ONLY BE PARALLEL IN Z DIRECTION!!" << endl;

  radius = 1;
  RealVariable *ar = dynamic_cast<RealVariable*>(OPAL.find("RADIUS"));
    
  if (ar) {
    radius = (int)ar->getReal();
    *gmsg << "RADIUS=" << radius << endl;
  } else
    *gmsg << "using std RADIUS=" << radius << endl;


  //new improved super MAP
  NDIndex<3> idx =  layout_m->getLocalNDIndex();
  Map = IPPLToMap3D(idx);
 
  //calculate nrx and nry from radius
  GaleriList.set("nx", radius*nr_m[0]);
  GaleriList.set("ny", radius*nr_m[1]);
  GaleriList.set("nz", nr_m[2]);
  A = CreateCrsMatrix("Laplace3D", Map, GaleriList);
  //if we have a solution from a previous call to this function, 
  //use LHS (containg old solution) else init LHS with 0
  //LHS = new Epetra_Vector(A->Map());
  LHS = new Epetra_Vector(*Map);
  LHS->PutScalar(0.0);
  
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
    nr_m[i] = domain_m[i].length();
  }
  
  if(decomp[0] == 1 || decomp[1] == 1) 
    *gmsg << "WARNING: MG SOLVER CAN ONLY BE PARALLEL IN Z DIRECTION!!" << endl;

  radius = 1;
  RealVariable *ar = dynamic_cast<RealVariable*>(OPAL.find("RADIUS"));
    
  if (ar) {
    radius = (int)ar->getReal();
    *gmsg << "RADIUS=" << radius << endl;
  } else
    *gmsg << "using std RADIUS=" << radius << endl;
  
  //calculate nrx and nry from radius
  GaleriList.set("nx", radius*nr_m[0]);
  GaleriList.set("ny", radius*nr_m[1]);
  GaleriList.set("nz", nr_m[2]);
  /*
  int mx=1,my=1,mz=1;
  NDIndex<3> idx =  layout_m->getLocalNDIndex();
  if(decomp[0] == 1) {
    mx = nr_m[0]/idx[0].length();
    cout << "parallel in x" << idx[0] << ": " << mx << " nodes" << endl;
  }
  if(decomp[1] == 1) {
    my = nr_m[1]/idx[1].length();
    cout << "parallel in y" << idx[1] << ": " << my << " nodes" << endl;
  }
  if(decomp[2] == 1) {
    mz = nr_m[2]/idx[2].length();
    cout << "parallel in z " << idx[2] << ": " << mz << " nodes" << endl;
  }

  GaleriList.set("mx", mx);
  GaleriList.set("my", my);
  GaleriList.set("mz", mz);
  Map = CreateMap("Cartesian3D", Comm, GaleriList);
  */
  NDIndex<3> idx =  layout_m->getLocalNDIndex();
  Map = IPPLToMap3D(idx);
  A = CreateCrsMatrix("Laplace3D", Map, GaleriList);
  //if we have a solution from a previous call to this function, 
  //use LHS (containg old solution) else init LHS with 0
  LHS = new Epetra_Vector(A->Map());
  LHS->PutScalar(0.0);
 
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

  //TODO: MATRIX FREE
  //use matrix stencil in computation directly
  
  //TODO: params (sane settings? read from input file?)
  int maxIterations = 100;
  double tol=1e-8;

  IpplTimings::startTimer(FunctionTimer5_m);
  NDIndex<3> localidx =  layout_m->getLocalNDIndex();
  cout << "CPU " << Comm.MyPID() << " " << localidx << endl;
  delete Map;
  Map = IPPLToMap3D(localidx);
  //we also have to redist LHS
  Epetra_Vector *tmp = new Epetra_Vector(*Map);
  Epetra_Import importer(*Map, LHS->Map());
  tmp->Import(*LHS, importer, Add);
  delete LHS;
  LHS = tmp;
  IpplTimings::stopTimer(FunctionTimer5_m);

  IpplTimings::startTimer(FunctionTimer1_m);
  /*      
  Vector_t hr2 = hr*hr;
  GaleriList.set("a", 2*(1/hr2[0]+1/hr2[1]+1/hr2[2]));
  GaleriList.set("b", -1/hr2[0]);
  GaleriList.set("c", -1/hr2[0]);
  GaleriList.set("d", -1/hr2[1]);
  GaleriList.set("e", -1/hr2[1]);
  GaleriList.set("f", -1/hr2[2]); //below corr to c
  GaleriList.set("g", -1/hr2[2]); //above corr to c
  */

  delete A;
  //A = CreateCrsMatrix("Cross3D", Map, GaleriList);
  //A = CreateCrsMatrix("Laplace3D", Map, GaleriList);
  A = Stencil3DLongitudinalNeumann(hr);
  IpplTimings::stopTimer(FunctionTimer1_m);

  //debug out
#ifdef DBG_STENCIL 
  EpetraExt::RowMatrixToMatlabFile("A.dat", *A);
#endif

  IpplTimings::startTimer(FunctionTimer2_m);

  Epetra_Vector RHS(A->Map());
  RHS.PutScalar(0.0);
  register int idx = 0;
  double fact = 1.0;
  int dx = (radius - 1)*nr_m[0]/2;
  int dy = (radius - 1)*nr_m[1]/2;

  for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
    for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
      idx = (y-localidx[1].first()) * localidx[2].length() + (x-localidx[0].first()) * (radius*nr_m[1] * localidx[2].length()) + dy*localidx[2].length() + dx*radius*nr_m[1]*localidx[2].length();
      for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
          
        fact = 1.0;
        if(isOnFrontSurface(idx) || isOnBackSurface(idx)) {fact = 0.5;} 
        RHS.Values()[idx] = fact*rho[x][y][z].get();

        idx++;

      }
    }
  }

  IpplTimings::stopTimer(FunctionTimer2_m);
  cout << "idx = " << idx << " and RHS mysize = " << RHS.MyLength() << " on CPU " << Comm.MyPID() << endl;

  //use old LHS solution as start vector
  //Epetra_LinearProblem Problem(A, LHS, &RHS);
  AztecOO solver(A, LHS, &RHS);

  /////////////////////////////////////////
  //PROBLEM SETUP FINISHED -- PRECOND BELOW
  
  IpplTimings::startTimer(FunctionTimer3_m);
  // set defaults
  //int *options = new int[AZ_OPTIONS_SIZE];
  //double *params = new double[AZ_PARAMS_SIZE];
  //ML_Epetra::SetDefaults("SA", MLList, options, params);
  ML_Epetra::SetDefaults("SA", MLList);
  // maximum number of levels
  MLList.set("max levels", 5);
  MLList.set("increasing or decreasing", "decreasing");
  // MG cycle type
  MLList.set("prec type", "MGW");

  //set aggregation scheme
  MLList.set("aggregation: type", "Uncoupled");

  // set smoother type
  //MLList.set("smoother: type", "symmetric Gauss-Seidel");
  MLList.set("smoother: type","Chebyshev");
  MLList.set("smoother: sweeps",3);
  MLList.set("smoother: pre or post", "both");

  // set coarse level solver
  MLList.set("coarse: type", "Amesos-KLU");

  // try to optimize mem for xt3
  //MLList.set("low memory usage", true);
  //maybe also helps
  //MLList.set("coarse: max size", 1024);

  // create the preconditioner object and compute hierarchy
  // true -> create the multilevel hirarchy
  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  ////////////////////////////////////////
  //PRECOND SETUP FINISHED -- SOLVER BELOW

  //CG SOLVER 
  solver.SetPrecOperator(MLPrec);

  //solver.SetAllAztecOptions(options);
  //solver.SetAllAztecParams(params);
  //solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_solver, AZ_cg);
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
  for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
    for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
      idx = (y-localidx[1].first()) * localidx[2].length() + (x-localidx[0].first()) * (radius*nr_m[1] * localidx[2].length()) + dy*localidx[2].length() + dx*radius*nr_m[1]*localidx[2].length();
      for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
        
        NDIndex<3> l(Index(x,x),Index(y,y),Index(z,z));
        rho.localElement(l) = LHS->Values()[idx]/(hr[0]*hr[1]*hr[2]);
        
        idx++;

      }
    }
  }

  IpplTimings::stopTimer(FunctionTimer4_m);
  cout << "idx = " << idx << " and LHS mysize = " << LHS->MyLength() << " on CPU " << Comm.MyPID() << endl;

  //cleanup
  delete MLPrec;
  //delete options;
  //delete params;

  //and we're hopefully done
}

inline bool MGPoissonSolver::isOnFrontSurface(const int idx) 
{

  int nxy = radius*nr_m[0]*radius*nr_m[1];
  int ixy = idx % nxy;
  return ((idx - ixy) / nxy) == 0;

}

inline bool MGPoissonSolver::isOnBackSurface(const int idx) 
{

  int nxy = radius*nr_m[0]*radius*nr_m[1];
  int ixy = idx % nxy;
  return ((idx - ixy) / nxy) == nr_m[2]-1;

}

void MGPoissonSolver::GetNeighbours(const int idx, int& W, int& E, int& S, int& N, int& F, int& B, int& numOutOfDomain)
{ 

  int ixy, iz, iy, ix; 
  int nx = radius*nr_m[0];
  int ny = radius*nr_m[1];
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

inline Epetra_Map* MGPoissonSolver::IPPLToMap3D(NDIndex<3> localidx)
{
  
  //TODO: use radius
  //we only parallelize in Z!
  int NumMyElements = radius*localidx[0].length()*radius*localidx[1].length()*localidx[2].length();
  int dx = (radius - 1)*nr_m[0]/2;
  int dy = (radius - 1)*nr_m[1]/2;
  vector<int> MyGlobalElements(NumMyElements);
  int idx = 0;
 
  /*
  for (int x = 0; x < dx; x++)
    for (int y = 0; y < dy; y++)
      for (int z = localidx[2].first(); z <= localidx[2].last(); z++)
        MyGlobalElements[idx++] = x + y * radius*nr_m[0] + z * (radius*nr_m[0] * radius*nr_m[1]);

  for (int x = localidx[0].last(); x < localidx[0].last()+dx; x++)
    for (int y = localidx[1].last(); y < localidx[1].last()+dy; y++)
      for (int z = localidx[2].first(); z <= localidx[2].last(); z++)
        MyGlobalElements[idx++] = x + y * radius*nr_m[0] + z * (radius*nr_m[0] * radius*nr_m[1]);

  for (int x = localidx[0].first(); x <= localidx[0].last(); x++)
    for (int y = localidx[1].first(); y <= localidx[1].last(); y++)
      for (int z = localidx[2].first(); z <= localidx[2].last(); z++)
        MyGlobalElements[idx++] = x + y * radius*nr_m[0] + z * (radius*nr_m[0] * radius*nr_m[1]);
  */

  for (int x = 0; x < radius*localidx[0].length(); x++)
    for (int y = 0; y < radius*localidx[1].length(); y++)
      for (int z = localidx[2].first(); z <= localidx[2].last(); z++)
        MyGlobalElements[idx++] = x + y * radius*nr_m[0] + z * (radius*nr_m[0] * radius*nr_m[1]);

  return(new Epetra_Map (-1, NumMyElements, &MyGlobalElements[0], 0, Comm));
}

Inform &MGPoissonSolver::print(Inform &os) const
{
  os << "* *************** M G P o i s s o n S o l v e r ************************************ " << endl;
  os << "* h " << hr_m << '\n';
  os << "* ********************************************************************************** " << endl;
}

///////////////////////////////////// OLD CODE ////////////////////////////////////////////

// setup a coordinate mapping and set them 
// in the ML parameter list. in a second step
// prepare the RHS of the problem 
// into an epetra vector
/*
void MGPoissonSolver::createMLGrid(Vector_t hr) 
{
  x_coord = new double[nr_m[0]];
  y_coord = new double[nr_m[1]];
  z_coord = new double[nr_m[2]];

  //TODO: the following only correct for CPU 0
  for(int i=0; i < nr_m[0]; i++)
    x_coord[i] = i*hr[0];
  
  for(int i=0; i < nr_m[1]; i++)
    y_coord[i] = i*hr[1];

  for(int i=0; i < nr_m[2]; i++)
    z_coord[i] = i*hr[2];

  MLList.set("x-coordinates", x_coord); 
  MLList.set("y-coordinates", y_coord); 
  MLList.set("z-coordinates", z_coord); 
  MLList.set("nx", nr_m[0]);
  MLList.set("ny", nr_m[1]);
  MLList.set("nz", nr_m[2]);
}
*/
/*
void MGPoissonSolver::createCoordinateGrid(Epetra_VbrMatrix* &Grid, Vector_t hr) {

  const Epetra_BlockMap& RowMap = Grid->RowMap();
  const Epetra_BlockMap& ColMap = Grid->ColMap();

  int NumMyRowElements = RowMap.NumMyElements();
  int NumMyColElements = ColMap.NumMyElements();

  int* MyGlobalRowElements = RowMap.MyGlobalElements();
  int* MyGlobalColElements = ColMap.MyGlobalElements();

  Epetra_Map PtrRowMap(-1, NumMyRowElements, MyGlobalElements, 0, Comm());
  Epetra_Map PtrColMap(-1, NumMyColElements, MyGlobalElements, 0, Comm());

  Epetra_Vector RowX(PtrRowMap);
  Epetra_Vector RowY(PtrRowMap);
  Epetra_Vector RowZ(PtrRowMap);

  //TODO: this only works like this for CPU 0
  Vector_t start_coords;
  start_coords = 0.0;
  for(int i=0; i < Comm.NumProc(); i++) {
    start_coords[0] += x_nodes_proc[i];
    start_coords[1] += y_nodes_proc[i];
    start_coords[2] += z_nodes_proc[i];
  }
  
  start_coords *= hr;
  
  for(int i=0; i < NumMyRowElements; i++) {
    RowX[i] = start_coords[0] + hr[0]*i;
    RowY[i] = start_coords[1] + hr[1]*i;
    RowZ[i] = start_coords[2] + hr[2]*i;
  }

  //non local nodes
  //when using Uncoupled aggregation this may be obsolete
  Epetra_Vector ColX(PtrColMap);
  Epetra_Vector ColY(PtrColMap);
  Epetra_Vector ColZ(PtrColMap);
  Epetra_Import Importer(PtrColMap, PtrRowMap);
  ColX.Import(RowX, Importer, Insert);
  ColY.Import(RowY, Importer, Insert);
  ColZ.Import(RowZ, Importer, Insert);

  int MaxNnz = Grid->MaxNumEntries();
  std::vector<int> colIndex(MaxNnz);
  std::vector<double> colValue(MaxNnz);
  std::vector<double> coord_i(3);
  std::vector<double> coord_j(3);

  Grid->PutScalar(0.0);

  for(int LocalRow = 0; LocalRow < NumMyRowElements; LocalRow++) {
    
    int RowDim, NumBlockEntries;
    int* BlockIndicies;
    Epetra_SerialDenseMatrix** RowValues;

    coord_i[2] = RowZ[LocalRow];
    coord_i[1] = RowY[LocalRow];
    coord_i[0] = RowX[LocalRow];

    Grid->ExtractMyBlockRowView(LocalRow, RowDim, NumBlockEntries, BlockIndicies, RowValues);

    double total = 0.0;

    for(int BlockEntry = 0; BlockEntry < NumBlockEntries; BlockEntry++) {

      int LocalCol = BlockIndicies[BlockEntry];

      //NO diag now
      if(LocalCol != LocalRow) {

        coord_j[2] = ColZ[LocalCol];
        coord_j[1] = ColY[LocalCol];
        coord_j[0] = ColX[LocalCol];

        //save the sqrt
        double d2 = (coord_i[0] - coord_j[0])*(coord_i[0] - coord_j[0]);
        d2 += (coord_i[1] - coord_j[1])*(coord_i[1] - coord_j[1]);
        d2 += (coord_i[2] - coord_j[2])*(coord_i[2] - coord_j[2]);

        if(d2 == 0) {
          cout << "ERROR d2 = 0.0";
          d2 = 1.0;
        }

        for(int k=0; k < RowValues[BlockEntry]->M(); k++) {
          for(int h=0; h, RowValues[BlockEntry]->N(); h++) {
            (*RowValues[j])(k,h) = -1.0/d2;
          }
        }

        total += 1.0/d2;

      }

    }

    //now handle diag: move into if -> else
    bool hasDiagBlock = false;
    int DiagBlock = 0;
    for(int j=0; j < NumBlockEntries; j++) {
      if(BlockIndicies[j] == LocalRow) {
        DiagBlock = j;
        hasDiagBlock = true;
      }
    }

    assert (hasDiagBlock == true);

    for(int k=0; k < RowValues[DiagBlock]->N(); k++)
      (*RowValues[DiagBlock](k,k) = total;

  }

}
*/

#endif /* HAVE_ML_SOLVER */
