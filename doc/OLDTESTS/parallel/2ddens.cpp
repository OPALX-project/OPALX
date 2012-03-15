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

// test program to demonstrate use of parallel reduction
#include "Ippl.h"

int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);
  Inform msg("2ddens ");

  int myN = Ippl::myNode();
  //
  // 2d array reduction 
  //
  const unsigned int BINS = 3;
  
  double **localDen;
  double **globalDen;

  localDen = (double **)malloc(BINS * sizeof(double *));
  globalDen = (double **)malloc(BINS * sizeof(double *));
  for (unsigned int i=0; i<BINS; i++) {
     localDen[i] = (double *)malloc(BINS * sizeof(double)); 
     globalDen[i] = (double *)malloc(BINS * sizeof(double));
  }

  for (unsigned int i=0; i<BINS; i++) {
    for (unsigned int j=0; j<BINS; j++) {
      localDen[i][j] =  myN+1.0;
      globalDen[i][j] = 0.0;
    }
  }
    
  reduce(&(localDen[0][0]), &(localDen[0][0]) + BINS*BINS,
	 &(globalDen[0][0]), OpAddAssign());
  
  for (unsigned int i=0; i<BINS; i++) {
    for (unsigned int j=0; j<BINS; j+=3) {
      msg << globalDen[i][j]   << " " 
	  << globalDen[i][j+1] << " " 
	  << globalDen[i][j+2];
    }
    msg << endl;
  }

  free(localDen);
  free(globalDen);
  return 0;
}

/***************************************************************************
 * $RCSfile: reduce.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:02:26 $
 * IPPL_VERSION_ID: $Id: reduce.cpp,v 1.1.1.1 2000/11/30 21:02:26 adelmann Exp $ 
 ***************************************************************************/
