// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <iostream.h>
#include <fstream.h>
    
#include <math.h>
#include <time.h>
      
#include "source/extpde.h"


          /* Poisson on a Ball */


int main()
{
  cout.precision(12);
  cout.setf(ios::fixed,ios::floatfield);
  ifstream PARAMETER;
  int iteration, i, n;
  double delta, delta_prime, beta, tau, eps, normi;
  double h;

  Index1D I, J;

  I = I.son(Ld);
  cout << " ind: "; I.Print(); cout << endl;
  
  J = I.next(Ld,1);

  cout << " ind: "; J.Print(); cout << endl;

  //  return 0;

  PARAMETER.open("para.dat",ios :: in);
  PARAMETER >> n >> iteration >> eps;
  PARAMETER.close();

  cout << "\n Solving Poisson's equation \n";
  cout << " with Dirichlet boundary conditions on on a ball by cg: ";
  cout << "\n Depth of grid: " << n;
  cout << "\n Maximal number of iterations: " << iteration;
  cout << "\n Stopping parameter for cg: " << eps << endl;

  Ball ball;
  Square square;
  Periodic_cylinder peri_cylinder;

  //  Domain domain(&ball);                 // Domain is a ball
  //  Domain domain(&square);                 // Domain is a cube

  D3vector move(0.25,0.25,0.0), scale(0.5,0.5,1.0);

  //  Domain domain(&peri_cylinder);                 // Domain is a cube
  Domain domain((peri_cylinder*scale)+move);                 // Domain is a cube

  /*
  Grid_gen_parameters ggen_parameter;
  //  ggen_parameter.Set_grid_point(D3vector(0.35,0.55,0.45)); 
  ggen_parameter.Set_offset_square(0.3);                    // larger than 0.0
  ggen_parameter.Set_stretch_square(2.2);                   // larger than 1.0
  */

  Grid grid(n,&domain);                  // Construction of grid

  Variable f(&grid);
  Variable g(&grid);
  Variable d(&grid);
  Variable r(&grid);
  Variable e(&grid);
  Variable u(&grid);
  Variable u_exakt(&grid);

  grid.Initialize();                 // Initialize storage

  grid.Print_poi_Info();   


  
  Function Sin(sin);                   // Definition of functions
  Function Sinh(sinh);

  //  u_exakt = Sin(X())*Sinh(Y())*Z();    // Exact solution
  //  u_exakt = Sin(X())*Sinh(Y());    // Exact solution

  u_exakt = Sin(2.0*M_PI*Z())*X()*Y();    // Exact solution

  /*
  // test
  
  u = X()*X() | interior_points;
  f = Laplace_FE(u) | interior_points;
  normi = product(f,f);
  cout << " Normi: " << normi;
  
  // test ende
  */
  /*
  cout << " number_nearb_points[2] :" << grid.number_nearb_points[2]
       << endl;
  */

  
  u = 1.0;
  normi = product(u,u);
  cout << normi;

  u = u_exakt | boundary_points;
  u = 0.0     | interior_points;

  //  f = 0.0;                             // Right hand side
  f = 2.0*M_PI*2.0*M_PI*Sin(2.0*M_PI*Z())*X()*Y();    // Exact solution

  u_exakt = u_exakt | interior_points;

  f = f | interior_points;
  r = r | interior_points;
  e = e | interior_points;
  d = d | interior_points;
  g = g | interior_points;


  cout << "\n cg: " << endl;
  cout << "-------------------- " << endl;

  d = -1.0 * r &
  r == Laplace_FE(u) - Helm_FE(f);


  delta = product(r,r);

  //  cout << "\n delta: " << delta << endl;
  for(i=1;i<=iteration && delta > eps;++i) {
    cout << i << ". ";

    g = Laplace_FE(d);

    tau = delta / product(d,g);

    u =  u + tau * d &
    r == r + tau*g;
    delta_prime = product(r,r);
    beta = delta_prime / delta;    delta = delta_prime;
    d = beta*d - r;
    e = u - u_exakt;

    //    cout << "\n delta: " << delta << endl;

    cout << "Iteration: Error: Max:" << L_infty(e);
    cout << ",  L2: " << sqrt(product(e,e)/normi) << endl;
  }
  

  ofstream DATEI;
  // Print FE-solution in file with UCD_file format
  DATEI.open("domain.dat",ios :: out);
  Print_UCD(u,&DATEI);
  DATEI.close();

  cout << "\n end. " << endl;
}
