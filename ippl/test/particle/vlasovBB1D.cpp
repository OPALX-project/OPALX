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

 vlasovBB1D.cpp

 Based more or less on  Maria-Paz Zorzano code.
 
 This code integrates the coupled system of two Vlasov equations         
 describing the time evolution of the distribution function              
 of two colliding beams.                                          
 An alternating finite difference scheme is used.                    
 The machine is assumed to be linear (betatron rotation).             
 No damping is included in this case.                                 
 So far the only nonlinearity comes from the beam-beam force, acting periodically
 at the IP.                                                              
 The beam-beam force is evaluated exactly, for a flat beam.        


 (y=p_x=Q x')                                                             
 System of two coupled Vlasov equations ( if i=1-> j=2; if i=2 -> j=1):   
 d rho_i(x,y,n)/dn=                                                       
  - Q_i y d rho_i/dx +                                                    
  + Q_i x d rho_i/dy +                                                    
  + delta(t,IP) 4*pi*xi (p.v. int rho_j (x')/(x-x') dx')  d rho_i/dy      

 d rho_i(x,y,n)/dn=-A(x,t)d rho_i/dy - B(y)d rho_i/dx                     
 A(x,t)= -Q_i x out of the IP
 A(x,t)= -4*pi*xi (p.v. int rho_j (x')/(x-x') dx')  d rho_i/dy at the IP 
 B(y)= Q_i y  

 usage: vlasovBB1D 
 
***************************************************************************/

#include "Ippl.h"
#include <iostream>


// dimension of f
const unsigned Dim = 2;

// presision
typedef double T;

const T pi = acos(-1.0);

#define numx 81
#define numy 81

//grid spacing in x and y
#define deltax 0.15
#define deltay 0.15

#define numturns 16400    // number of turns
#define N 108             // every N steps equals one turn
#define numt (numturns*N) // number of times steps 

//the first 'a' step(s) are used for the beam-beam kick at the IP 
//(N-a) steps to integrate the linear part (rotation only)
#define a 3
#define dt (2.0*pi/(N-a)) //In these units T=2*pi for one turn, 

//BeamBeam parameters

//Natural frequencies of the system

#define Q1 0.145
#define Q2 0.109
#define Nb1 (2*1.05e11)
#define Nb2 Nb1
#define Eb 7.0e6 //MeV  particle energy
#define Eo 928.0 //MeV  rest energy
#define gamma (Eb/Eo)
#define rp 1.551e-18
#define betah 0.5 // m
#define betav 0.5 // m
#define rms_x 16.0e-6//m
//4*pi*xi= a*rp*Nb1/gamma/rms_x/rms_x*betah=4*pi*0.02 for these parameters 

//Initial condition
#define sigmax2 1.0 //Gaussian distribution variance in x
#define sigmay2 1.0 //Gaussian distribution variance in y
#define cix1 0.03//<x> for beam 1 Gaussian distrib.
#define ciy1  0.0//<y> for beam 1 Gaussian distrib.
#define cix2  0.0//<x> for beam 2 Gaussian distrib.
#define ciy2  0.0//<y> for beam 2 Gaussian distrib.

#define ti 0.0 //starting time




int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  
  n = atoi(argv[1]);

  msg << "Dim=" << Dim << " n= "<< n << " M= " << pow(Dim,n) << endl;

  // create layout objects for max DIM==6
  const Index X1(n), X2(n);          // spacial
  const Index P1(n), P2(n);          // momenta

  // Initialize domain objects
  
  NDIndex<Dim> domain1, domain2;

  domain[0] = n;
  domain[1] = n;
  domain[2] = n;
  domain[3] = n;
  
  domain1[0] = n+1;
  domain1[1] = n+1;
  domain1[2] = n+1;
  domain1[3] = n+1;

  if (Dim==6) {
    domain[4] = n;
    domain[5] = n;
    domain1[4] = n+1;
    domain1[5] = n+1;
  }

  UniformCartesian<Dim> mymesh(domain1);
  FieldLayout<Dim> FL(domain);

  int testNum = 1;
  double eps = 1.0e-07;

  // Flags to keep track of pass/fail results:
  bool passed = true;           // True at end means passed all tests
  bool passedSingleTest = true; // For individual tests

  Timer timer1;

  timer1.clear();
  timer1.start();
  Field<T,Dim> A(mymesh,FL);
  Field<T,Dim> B(mymesh,FL);
  Field<T,Dim> C(mymesh,FL);
  timer1.stop();
  msg << "Static field objects created t= " << timer1.clock_time() << endl;

  Field<double,Dim>::iterator fIt;
  T a,b,c;
  

  // TEST 1: testing accumulation: Field += Field
  a = 0.0;   b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;   B = 1.0;
  A += B;
  timer1.stop();

  a += b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 2: testing substraction: Field -= Field  
  a = 0.0;   b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;   B = 1.0;
  A -= B;
  timer1.stop();

  a -= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 3: testing multipplication: Field *= Field  
  a = 1.0;   b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;   B = 2.0;
  A *= B;
  timer1.stop();

  a *= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 4: testing multipplication: Field /= Field    
  a = 1.0;   b = 5.0;
  timer1.clear();
  timer1.start();
  A = 1.0;   B = 5.0;
  A /= B;
  timer1.stop();

  a /= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 5: testing scalar addition: Field += T
  a = 1.0; 
  timer1.clear();
  timer1.start();
  A = 0.0; 
  A += a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 6: testing scalar substraction: Field -= T
  a = 1.0; 
  timer1.clear();
  timer1.start();
  A = 2.0; 
  A -= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 7: testing scalar multipplication: Field *= T
  a = 2.0; 
  timer1.clear();
  timer1.start();
  A = 1.0; 
  A *= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 8: testing scalar division: Field /= T
  a = 2.0; 
  timer1.clear();
  timer1.start();
  A = 4.0; 
  A /= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 9: testing field expression and accumulation: Field += ExprElem  
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  C = 2.0;
  A += B + C;
  timer1.stop();

  a += b + c;
  
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 10: testing field expression and substraction: Field -= ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  C = 2.0;
  A -= B + C;
  timer1.stop();

  a -= b + c;
  
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 11: testing field expression and multipplication: Field *= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  C = 3.0;
  A *= B + C;
  timer1.stop();

  a *= b + c;
  
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 12: testing field expression and division: Field /= ExprElem
 
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  C = 3.0;
  A /= B + C;
  timer1.stop();

  a /= b + c;
  
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 13: testing accumulation with indexed field: IndexingField += IndexingField
  a = 0.0;
  b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  A[I][J][K][L][M][N] += B[I][J][K][L][M][N];
  timer1.stop();
  a += b;
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------
  

  // TEST 14: testing substraction with indexed field: IndexingField -= IndexingField
  a = 0.0;
  b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  A[I][J][K][L][M][N] -= B[I][J][K][L][M][N];
  timer1.stop();
  a -= b;
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 15: testing multipplication with indexed field: IndexingField *= IndexingField
  a = 1.0;
  b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  A[I][J][K][L][M][N] *= B[I][J][K][L][M][N];
  timer1.stop();
  a *= b;
  passedSingleTest = true;
  
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  
  // TEST 16: testing division with indexed field: IndexingField /= IndexingField
  a = 1.0;
  b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  A[I][J][K][L][M][N] /= B[I][J][K][L][M][N];
  timer1.stop();
  a /= b;
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 17: testing accumulation with indexed field: IndexingField += T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  A[I][J][K][L][M][N] += a;
  timer1.stop();
  
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 18: testing substraction with indexed field: IndexingField -= T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 2.0;
  A[I][J][K][L][M][N] -= a;
  timer1.stop();
  
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 19: testing multipplication with indexed field: IndexingField *= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  A[I][J][K][L][M][N] *= a;
  timer1.stop();
  
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 20: testing division with indexed field: IndexingField /= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 4.0;
  A[I][J][K][L][M][N] /= a;
  timer1.stop();
  
  passedSingleTest = true;
 
  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 21: testing accumulation index: IndexingField += Index
  a = 50.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  A[I][J][0][0][0][0] += I;
  timer1.stop();
  T sumval = sum(A);
  passedSingleTest = true;
  
  if (fabs(sumval - a) > eps) {
    passedSingleTest = false; passed = false;
    fIt = A.end();
  }
  
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 22: testing substraction index: IndexingField -= Index
  // TEST 23: testing multipplication index: IndexingField *= Index
  // TEST 24: testing division index: IndexingField /= Index

  testNum = 25;
  
  /*
    Note: have to use reduce in order to deduce the result
          of the "local" compares of the field values
  */

  NDIndex<Dim> idx = FL.getLocalNDIndex();
  NDIndex<Dim> elem;

  // TEST 25: testing accumulation with IndexingFields on rhs: IndexingField += ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] += B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();
  
  a += b + c;
  
  passedSingleTest = true;
   
  for (unsigned int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (unsigned int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
	passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  
  // TEST 26: testing substraction with IndexingFields on rhs: IndexingField -= ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] -= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();
  
  a -= b + c;
  
  passedSingleTest = true;
    
  for (unsigned int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (unsigned int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
	passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 27: testing multipplication with IndexingFields on rhs: IndexingField *= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] *= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();
  
  a *= b + c;
  
  passedSingleTest = true;
   
  for (unsigned int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (unsigned int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
	passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 28: testing divisiond with IndexingFields on rhs: IndexingField /= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] /= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();
  
  a /= b + c;
  
  passedSingleTest = true;
   
  for (unsigned int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (unsigned int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
	passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  Ippl::Comm->barrier();
  msg << "done ...." << endl;
  return 0;
}























/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

