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

int main(int nnn,char *aaa[]) { 
  int **a,i1,i2,size1,size2; /* Declare variables. */ 
  size1=50; size2=40; i1=20; i2=30; /* Set variables. */ 
  a=malloc((size_t)(size1*size2*sizeof(int))); /* Allocate memory. */ 
  if(!a) Error(); /* See if we're out of memory. */ 
  a[ size2*i1 + i2 ] = 10; /* Set "a[i1,i2]" to 10. */ 
  free(a); /* Free memory. */ 
} 


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

