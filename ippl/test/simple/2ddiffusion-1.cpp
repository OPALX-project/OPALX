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
 * Note: need directory data on the same level than the program is executed
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>

#include "Ippl.h"

#define Dim 2



void saveField(Field<double,Dim> A, int iteration, double dx, double dy) { 

    /**
       Works only serial!
    */

    const NDIndex<Dim> lDom = A.getLayout().getLocalNDIndex();

    const int nx = lDom[0].length();
    const int ny = lDom[1].length();
    
    std::ofstream vtkout;
    vtkout.precision(10);
    vtkout.setf(ios::scientific, ios::floatfield);

    std::stringstream fname;
    fname << "data/u_";
    fname << setw(4) << setfill('0') << iteration;
    fname << ".vtk";

    vtkout.open(fname.str().c_str(), ios::out);

    vtkout << "# vtk DataFile Version 2.0" << endl;
    vtkout << "#u12" << endl;
    vtkout << "#ASCII" << endl;
    vtkout << "#DATASET STRUCTURED_POINTS" << endl;
    vtkout << "#DIMENSIONS " << nx << " " << ny << "   " << 1 << endl;
    vtkout << "#ORIGIN 0 0 0" << endl;
    vtkout << "#SPACING " << 1 << " " << 1 << "   " << 1 << endl;
    vtkout << "#POINT_DATA " << nx*ny << endl;
    vtkout << "#SCALARS TEMPERATURE float" << endl;
    vtkout << "#LOOKUP_TABLE default" << endl;

    for(int y=lDom[1].first(); y<lDom[1].last(); y++) {
	for(int x=lDom[0].first(); x<lDom[0].last(); x++) {
	    const double val = A[x][y].get();
	    vtkout << x << "  " << y << "  " << val << endl;
	}
    }
    vtkout.close();
} 


int main(int argc, char *argv[])
{
    Ippl ippl(argc,argv);
    Inform msg ("u12 ");

    const int sizeX      = 100;
    const int sizeY      = 100;
    const int iterations = 100;

    Index I(sizeX);  Index J(sizeY);
  
    FieldLayout<Dim> layout(I,J);

    Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));

    A = 0.0;
    A[sizeX/2][sizeY/2] = 1.0;

    double fact = 1.0/9.0;

    for(int iter = 0 ; iter < iterations ; iter++ ) {
	/*
	A[I][J] = fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] + 
			A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] + 
			A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1]);
	*/
	
	A[I][J] = fact*(A[I+1][J  ] + A[I-1][J  ] - 4*A[I][J] + 
			A[I  ][J+1] + A[I  ][J-1]); 

	saveField(A, iter, sizeX, sizeY);  
    
    }
    
    msg << "Done ; A[sizeX/2][sizeX/2] = " << A[sizeX/2][sizeX/2] << endl;
	
    return 0;
}
/***************************************************************************
 * $RCSfile: u10.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: u10.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
