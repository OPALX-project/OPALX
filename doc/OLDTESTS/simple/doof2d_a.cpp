// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include <iostream.h>
#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);

  cout << "input the number of vnodes" << endl;
  int vnodes;
  cin >> vnodes; 
  int sizeX;
  cout << "input sizeX" << endl;
  cin >> sizeX;
  int sizeY;
  cout << "input sizeY" << endl;
  cin >> sizeY;
  const unsigned Dim=2;
  Index I(sizeX);
  Index J(sizeY);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL, vnodes);
  Field<double,Dim> A(layout);
  FieldBlock<double,Dim> io("atest",layout);
  FieldView<double,Dim> plotA(A);
  
  int iter;
  do {
    cout << "input iter" << endl;
    cin >> iter;
    if( iter < 0 ) break;
    if( iter < io.get_NumRecords() ) {
      io.read(A,0,iter);
      plotA.view();
    }
  } while(iter >= 0 );

  cout << "enter to continue" << endl;
  int in;
  cin >> in;
  return 0;
}
/***************************************************************************
 * $RCSfile: doof2d_a.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:47 $
 * IPPL_VERSION_ID: $Id: doof2d_a.cpp,v 1.1.1.1 2000/11/30 21:05:47 adelmann Exp $ 
 ***************************************************************************/
