/***************************************************************************
 *
 * 
 * 
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
 * $RCSfile: t5.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t5.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $ 
 ***************************************************************************/
