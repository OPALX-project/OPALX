// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <iostream>
#include <fstream>
    
#include <math.h>
#include <time.h>
 
#include "source/extpde.h"


void Gauss_Seidel(Variable& u, Variable& f, int iter,
		  Res_stencil_boundary& ResS) {
  for(int i=0;i<iter;++i) {
    u = u - (ResS(Laplace_FE(u)) +f) / ResS[Diag_Laplace_FE()];
  }
}

void MG(Variable& u, Variable &f, Variable& r,
        int k, int iter, Res_stencil_boundary& ResS) {

  if(k > u.Give_grid()->Min_level()) {
    // pre-smoothing
    Gauss_Seidel(u,f,iter,ResS);
    
    // restriction
    r = f + ResS(Laplace_FE(u));
    f = Restriction_FE(r);
    
    u.Level_down();
    u = 0.0;

    // coarse-grid correction
    MG(u,f,r,k-1,iter,ResS);

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
         

int main(int argc, char** argv)
{
  cout.precision(8);
  cout.setf(std::ios::fixed,std::ios::floatfield);
  ifstream PARAMETER;
  int iteration, i, iter_relax;
  double normi;
  double sp;
  int r_parallel;

  double starttime, endtime;
  double start_calc_time, start_MG_time, end_calc_time;
  double o_square, s_square;
  
  double meshsize;

  int my_rank;        /* Rank of process */
  int p;              /* Number of processes */

  double a,b,l;       // defines the cylinder

  PARAMETER.open("para.dat",std::ios :: in);
  PARAMETER >> meshsize  >> iteration >> iter_relax 
            >> r_parallel  >> o_square >> s_square >> a >> b >> l;
  PARAMETER.close();

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  starttime = MPI_Wtime();
  ofstream DATEI;

  if(my_rank==0) {
    cout << "\n Solving Poisson's equation \n";
    cout << "\n Geometry cylinder (" << a << "," << b << "," << l << ")";
    cout << " with Dirichlet boundary conditions by MG: ";
    cout << "\n Meshsize of grid: " << meshsize;
    cout << "\n Number of relaxations: " << iter_relax;
    cout << "\n Maximal number of iterations: " << iteration << endl;
  }    


  D3vector Vs(2.0*a,2.0*b,2.0*l);
  D3vector Vm(-1.0,-1.0,-1.0);
  Square square;
  Cylinder   cylinder;
  Ball ball;
  Domain domain(&square); 
  // Domain domain(&cylinder*Vs+Vm);    // Domain is a cylinder stretched by Vs and translated by Vm
 
  Grid_gen_parameters ggen_parameter;
  ggen_parameter.Set_r_parallel(r_parallel);           // choose larger than 1
  ggen_parameter.Set_offset_square(o_square);
  ggen_parameter.Set_stretch_square(s_square);                    
       
  Grid grid(meshsize,&domain,
	    ggen_parameter,MPI_COMM_WORLD); // Construction of grid

  if(my_rank==0) grid.process_info();
  
  Variable f(&grid);
  Variable r(&grid);
  Variable e(&grid);
  Variable u(&grid);
  Variable u_exakt(&grid);
  
  Res_stencil_boundary ResS(&grid);           // definition of a stencil 
  
  grid.Initialize();                 // Initialize storage


  BrickInfo *brick_info;

  brick_info = grid.getBrickInfo();

  if(my_rank==0) 
    cout << *(brick_info) << endl;

  start_calc_time = MPI_Wtime();
    
    
  grid.Print_poi_Info();
    
  if(my_rank==0) DATEI.open("hex.dat",std::ios::out);
  grid.Print_region_processes_UCD(&DATEI);
  if(my_rank==0) DATEI.close();
 
  ///////////   Test der Funktion mit 4 Parametern
  Function Cos(cos);                   // definition of a function
      
  u_exakt = Cos(X()*M_PI)*Cos(Y()*M_PI)*Cos(Z()*M_PI);   // exact solution
  
  normi = 0.0;
  u = 1.0;
  normi = product(u,u);
  
  if(my_rank==0) cout << " Number of grid points: " << normi << endl;
  
  u = u_exakt | boundary_points;
  u = 0.0     | interior_points;
  
  r = Cos(X()*M_PI)*Cos(Y()*M_PI)*Cos(Z()*M_PI)     // right hand side
    * -3.0 * M_PI * M_PI; 
  f = Helm_FE(r)    | interior_points;
  r = 0.0           | all_points;
  r = r             | interior_points;
  u_exakt = u_exakt | interior_points;
  
  f = f | interior_points;
  r = r | interior_points;
  e = e | interior_points;
  
  Gauss_Seidel(u,f,1,ResS);
  
  if(my_rank==0) cout << "\n MG: " << endl;
  if(my_rank==0) cout << "-------------------- " << endl;
  
  start_MG_time = MPI_Wtime();    
  
  for(i=0;i<iteration;++i) {
    if(my_rank==0) cout << i << ". ";
    MG(u,f,r,u.Give_grid()->Max_level(),iter_relax,ResS);
    
    e = u - u_exakt;
    sp = L_infty(e);
    if(my_rank==0) cout << "Iteration: Error: Max:" << sp;
    sp = product(e,e);
    if(my_rank==0) cout << ",  L2: " << sqrt(sp/normi) << endl;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  end_calc_time = MPI_Wtime();
  
  //  if(my_rank==0) DATEI.open("domain.dat",std::ios:: out);
  //  Print_UCD(u,&DATEI);
  //  if(my_rank==0) DATEI.close();
      
  if(my_rank==0) cout << " \n End: " << endl;
      
  endtime = MPI_Wtime();
      
      
  if(my_rank==0) {
    cout << "\n total time: " << endtime-starttime
	 << " sec. " << endl;
    cout << "\n calc time: " << end_calc_time-start_calc_time
	 << " sec. " << " \n end " << endl;
    cout << "\n elastime time: " << end_calc_time-start_MG_time
	 << " sec. " << " \n end " << endl;
  }
  MPI_Finalize();
}
