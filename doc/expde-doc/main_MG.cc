// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <iostream.h>
#include <fstream.h>
    
#include <math.h>
#include <time.h>
      
#include "source/extpde.h"


          /* Poisson by multigrid with Dirichlet boundary conditions */

         
void Gauss_Seidel(Variable& u, Variable& f, int iter,
		  Res_stencil_boundary& ResS) {
  for(int i=0;i<iter;++i) {
    u = u - (ResS(Laplace_FE(u)) +f) / ResS[Diag_Laplace_FE()];
  }
}


void MG(Variable& u, Variable &f, Variable& r,
	int k, int max_level, int iter,
	Res_stencil_boundary& ResS) {
  if(k > u.Give_grid()->Min_interior_level()) {
    // pre-smoothing
    Gauss_Seidel(u,f,iter,ResS);

    // restriction
    r = f + ResS(Laplace_FE(u));
    f = Restriction_FE(r);
    
    u.Level_down();
    u = 0.0;
    
    // coarse-grid correction
    MG(u,f,r,k-1,max_level,iter,ResS);

    // prolongation
    f.Level_up();
    u = u + Prolongation_FE(u);
    
    // post-smoothing
    Gauss_Seidel(u,f,iter,ResS);
  }
  else {
    // smoothing on coarsest grid
    Gauss_Seidel(u,f,iter,ResS);
  }
}
                      

int main()
{
  cout.precision(12);
  cout.setf(ios::fixed,ios::floatfield);
  ifstream PARAMETER;
  int iteration, i, n, iter_relax;
  double normi;

  D3vector move(0.25,0.25,0.0), scale(0.5,0.5,1.0);
    
  PARAMETER.open("para.dat",ios :: in);
  PARAMETER >> n >> iteration >> iter_relax;
  PARAMETER.close();

  cout << "\n Solving Poisson's equation on a ball \n";
  cout << " with Dirichlet boundary conditions by multigrid: ";
  cout << "\n Depth of grid: " << n;
  cout << "\n Maximal number of iterations: " << iteration;
  cout << "\n Number of smoothing steps: " << iter_relax << endl;     

  Ball ball;
  Periodic_cylinder peri_cylinder;

  //  Domain domain(&ball);                 // Domain is a ball
  Domain domain((peri_cylinder*scale)+move);   

  Grid grid(n,&domain);                  // Construction of grid

  Variable f(&grid);
  Variable r(&grid);
  Variable u(&grid);
  Variable One(&grid);
  Variable u_exakt(&grid);

  Res_stencil_boundary ResS(&grid);  // definition of a stencil 
                                       // near the boundary

  grid.Initialize();                 // initialize storage  

  grid.Print_poi_Info();
         
  u = 1.0 | interior_points;
  normi = product(u,u);
  
  Function Cos(cos);                   // definition of a function
  Function Sin(sin);                   // Definition of functions
  Function Sinh(sinh); 

  cout << "\n Solving equation! \n";
  
  //  u_exakt = Cos(X()*M_PI)*Cos(Y()*M_PI)*Cos(Z()*M_PI); // exact solution

  u_exakt = Sin(2.0*M_PI*Z())*X()*Y();    // Exact solution
  // u_exakt = Sin(X())*Sinh(Y());    // Exact solution   


  u = u_exakt | boundary_points;
  u = 0.0     | interior_points;
  
  /*
  r = Cos(X()*M_PI)*Cos(Y()*M_PI)*Cos(Z()*M_PI) 
    * -3.0 * M_PI * M_PI;                                // right hand side
    */

  // r = 0.0;
  r = - 2.0*M_PI*2.0*M_PI*Sin(2.0*M_PI*Z())*X()*Y();    // Exact solution   

  f = Helm_FE(r)    | interior_points;
  r = 0.0           | all_points;
  r = r             | interior_points;
  u_exakt = u_exakt | interior_points;  
  
  cout << "\n MG: " << endl;
  cout << "-------------------- " << endl;
  
  for(i=1;i<=iteration;++i) {
    MG(u,f,r,n,n,iter_relax,ResS);
    r = u - u_exakt;
    cout << "Iteration: Error: Max:" << L_infty(r);
    cout << ",  L2: " << sqrt(product(r,r)/normi) << endl;        
  }

  
  f = 1.0 | all_points;   
  One = Helm_FE(f)  | all_points;        

  r = DX_FE(u) / One | all_points;      
  // r = DY_FE(u) / One | all_points;      
  // r = DZ_FE(u) / One | all_points;      

  u = r;
  

  ofstream DATEI; 
  // Print FE-solution in file with UCD_file format
  DATEI.open("domain.dat",ios :: out);
  Print_UCD(u,&DATEI);
  DATEI.close();  
             
  // Print FE-solution in file with dx_file format
  DATEI.open("domain.dx",ios :: out);
  Print_Dx(u,&DATEI);
  DATEI.close();               


  cout << "\n end " << endl;
}
