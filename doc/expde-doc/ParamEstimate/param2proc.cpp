/*
  
param2proc.cpp 

*/

#include <fstream>
#include <new>
#include <iostream>
#include <sstream>
#include "Ipl.h"

using namespace std;
  
#define T double
#define Dim 3


class ElemGeometry {
  
public:
  ElemGeometry() {};
  ElemGeometry(double ai, double bi, double li) : a(ai),b(bi),l(li) {};
  void print() { 
    Inform msg("ElemGeometry ");
    msg << "a= " << a << " b= " << b << " l= " << l << endl; 
  };
  double a;
  double b;
  double l;
};


int main(int argc, char ** argv)
{
  Inform msg("param2proc ");

  double a; double b; double l;   // geometry
  double maxErr;
  int iter_relax;
  double o_square, s_square; 
  int p_test;
  int r_parallel;
  long nPart;
  double normalConst;
  double hInit, r;
  int nMax;
  
  int maxNumProc;
  if (argc == 1) {
    cout << "Usage: <paramFile> <nProc> | <h> <rparallel> <osquare> <ssquare> <a> <b> <l> <nProc>";
    cout << endl;
    exit(1);  
  }
  else {
    if (argc == 3) {    
      ifstream fs;
      string filename(argv[1]);
      fs.open(filename.c_str(),std::ios::in);
      fs >> hInit  >> maxErr >> iter_relax >> r_parallel 
	 >> o_square >> s_square >> p_test >> a >> b 
	 >> l >> nPart >> r >> normalConst;
      fs.close();
      maxNumProc = atoi(argv[2]);
      cout << "read from " << filename;
      cout << " h = " << hInit << " rparallel= " << r_parallel << endl;
    }
    else {
      hInit = atof(argv[1]);
      r_parallel = atoi(argv[2]);
      o_square = atof(argv[3]);
      s_square = atof(argv[4]);
      a = atof(argv[5]);
      b = atof(argv[6]);
      l = atof(argv[7]);
      maxNumProc = atoi(argv[8]);
    }
  }

  ElemGeometry geom(a,b,l);
  geom.print();

  Grid_gen_parameters ggen_parameter;
  ggen_parameter.Set_r_parallel(r_parallel);    
  ggen_parameter.Set_offset_square(o_square);
  ggen_parameter.Set_stretch_square(s_square);          
  
  /*
    Domain is a cylinder stretched 
    by Vs and translated by Vm
  */ 

  D3vector Vs(geom.a,geom.b,geom.l);
  D3vector Vm(-geom.a/2.0,-geom.b/2.0,-geom.l/2.0);
  Cylinder cylinder;
  
  Domain geomDomain(&cylinder*Vs+Vm); 

  if (hInit>0.0) {
    Parallel_Info test_parallel_info1(hInit,&geomDomain,maxNumProc,ggen_parameter);
    test_parallel_info1.process_info();
  }else {
    nMax = (int) abs(hInit);
    Parallel_Info test_parallel_info2(nMax,&geomDomain,maxNumProc,ggen_parameter);
    test_parallel_info2.process_info();
  }
  return 0;
}
