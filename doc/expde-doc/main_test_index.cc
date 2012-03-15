// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <iostream.h>
#include <fstream.h>
    
#include <math.h>
#include <time.h>
      
#include "source/extpde.h"


                      

int main()
{
  cout.precision(12);
  cout.setf(ios::fixed,ios::floatfield);
  ifstream PARAMETER;
  int iteration, i, n, iter_relax;
  double normi;

  PARAMETER.open("para.dat",ios :: in);
  PARAMETER >> n >> iteration >> iter_relax;
  PARAMETER.close();


  Index1D I,J;

  cout << " ccord: " << I.coordinate() << endl;


  
  I = I.son(Ld);
  cout << " ccord: " << I.coordinate() << endl;
  I = I.son(Ld);
  cout << " ccord: " << I.coordinate() << endl;
  

  cout << " Test Neigh" << endl;
  J = I.neighbour(Ld);
  cout << " ccord: " << J.coordinate() << endl;

  J = I.neighbour(Rd);
  cout << " ccord: " << J.coordinate() << endl;

  cout << " Test Next 2" << endl;
  cout << " ccord: " << I.coordinate() << endl;
  J = I.next(LOrt,2);
  cout << " ccord: " << J.coordinate() << endl;

  J = I.next(ROrt,2);
  cout << " ccord: " << J.coordinate() << endl;


  cout << "\n end " << endl;
}
