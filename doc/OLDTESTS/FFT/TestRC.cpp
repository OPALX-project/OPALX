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
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// TestFFT.cpp , Tim Williams 1/27/1997
// Updated by Julian Cummings, 3/31/98

// Tests the use of the (parallel) FFT class.

#include "Ippl.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <complex>
using namespace std;
#else
#include <complex.h>
#endif

int main(int argc, char *argv[])
{
  TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);
  Ippl ippl(argc,argv);
  Inform testmsg(NULL,0);

  // Hardwire dimensionality by setting cpp definition ONED, TWOD, or THREED:
  // The preceding cpp definition causes compile-time setting of D:
  const unsigned D=3U;
  testmsg << "%%%%%%% Dimensionality: D = " << D << " %%%%%%%" << endl;
  
  unsigned ngrid[D];   // grid sizes

  // Used in evaluating correctness of results:
  double realDiff;
   
  // Various counters, constants, etc:
  int d;
  int Parent = 0;
  int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0);
  double pi = acos(-1.0);
  double twopi = 2.0*pi;
  // Timer:
  Timer timer;
  bool constInput;  // preserve input field in two-field transform

  // Layout information:
  unsigned vnodes;             // number of vnodes; input at cmd line
  e_dim_tag allParallel[D];    // Specifies SERIAL, PARALLEL dims
  for (d=0; d<D; d++) 
    allParallel[d] = PARALLEL;
  e_dim_tag serialParallel[D]; // Specifies SERIAL, PARALLEL dims
  serialParallel[0] = SERIAL;
  for (d=1; d<D; d++) 
    serialParallel[d] = PARALLEL;

  // Compression of temporaries:
  bool compressTemps;

  // ================BEGIN INTERACTION LOOP====================================
  while (true) {

    // read in vnodes etc. off of node 0
    if( Ippl::Comm->myNode() == Parent ) {
      bool vnodesOK = false;
      while (!vnodesOK) {
        cout << "Enter vnodes (0 to exit): ";
        cin >> vnodes;
        if (vnodes==0 || vnodes>=Ippl::getNodes())
          vnodesOK = true;
        else
          cout << "Number of vnodes should not be less than number of pnodes!"
               << endl;
      }
      if (vnodes > 0) {
	cout << "Compress temps? (Enter 1 for true or 0 false): ";
	int tempWorkaroundBoolIO;
	cin >> tempWorkaroundBoolIO;
	if (tempWorkaroundBoolIO == 0)
	  compressTemps = false;
        else
	  compressTemps = true;
	cout << "Constant input fields? (Enter 1 for true or 0 false): ";
	cin >> tempWorkaroundBoolIO;
	if (tempWorkaroundBoolIO == 0)
	  constInput = false;
        else
	  constInput = true;
	cout << "input nx, ny, and nz; space-delimited: ";
	for (d=0; d<D; d++) 
	  cin >> ngrid[d];
      }
      // now broadcast data to other nodes
      Message *mess = new Message();
      putMessage( *mess, vnodes );
      if (vnodes > 0) {
	putMessage(*mess, compressTemps);
	putMessage(*mess, constInput);
	for (d=0; d<D; d++) putMessage( *mess, ngrid[d] );
      }
      Ippl::Comm->broadcast_all(mess, tag);
    }
    // now each node recieves the data
    unsigned pe = Ippl::Comm->myNode();
    Message *mess = Ippl::Comm->receive_block(Parent, tag);
    PAssert(mess);
    getMessage( *mess, vnodes );
    if (vnodes <= 0) 
      break;
    getMessage(*mess, compressTemps);
    getMessage(*mess, constInput);
    for (d=0; d<D; d++) 
      getMessage( *mess, ngrid[d] );
    delete mess;
    
    // create standard domain
    NDIndex<D> ndiStandard;
    for (d=0; d<D; d++) 
      ndiStandard[d] = Index(ngrid[d]);
    // create new domain with axes permuted to match FFT output
    
    // create half-size domain for RC transform along zeroth axis
    NDIndex<D> ndiStandard0h = ndiStandard;
    ndiStandard0h[0] = Index(ngrid[0]/2+1);
  
     // all parallel layout, standard domain, normal axis order
    FieldLayout<D> layoutPPStan(ndiStandard,allParallel,vnodes);
    // zeroth axis serial, standard domain, normal axis order
    FieldLayout<D> layoutSPStan(ndiStandard,serialParallel,vnodes);
 
    // all parallel layout, zeroth axis half-size domain, normal axis order
    FieldLayout<D> layoutPPStan0h(ndiStandard0h,allParallel,vnodes);
    // zeroth axis serial, zeroth axis half-size domain, normal axis order
    FieldLayout<D> layoutSPStan0h(ndiStandard0h,serialParallel,vnodes);
       
    // create test Fields for complex-to-complex FFT
    BareField<dcomplex,D> CFieldPPStan(layoutPPStan);
   
    BareField<dcomplex,D> CFieldSPStan(layoutSPStan);
   
    BareField<double,D> diffFieldSPStan(layoutSPStan);
      
    // create test Fields for real-to-complex FFT
    BareField<double,D> RFieldSPStan(layoutSPStan);
    BareField<double,D> RFieldSPStan_save(layoutSPStan);
    BareField<dcomplex,D> CFieldSPStan0h(layoutSPStan0h);
 
    testmsg << "RFieldSPStan   layout= " << layoutSPStan << endl;
    testmsg << "CFieldSPStan0h layout= " << layoutSPStan0h << endl;


    // For calling FieldDebug functions from debugger, set up output format:
    setFormat(4,3);

    // Rather more complete test functions (sine or cosine mode):
    dcomplex sfact(1.0,0.0);      // (1,0) for sine mode; (0,0) for cosine mode
    dcomplex cfact(0.0,0.0);      // (0,0) for sine mode; (1,0) for cosine mode

    double xfact, kx, yfact, ky, zfact, kz;
    xfact = pi/(ngrid[0] + 1.0);
    yfact = 2.0*twopi/(ngrid[1]);
    zfact = 2.0*twopi/(ngrid[2]);
    kx = 1.0; ky = 2.0; kz = 3.0; // wavenumbers

    CFieldPPStan[ndiStandard[0]][ndiStandard[1]][ndiStandard[2]] = 
      sfact * ( sin( (ndiStandard[0]+1) * kx * xfact +
                      ndiStandard[1]    * ky * yfact +
                      ndiStandard[2]    * kz * zfact ) +
		sin( (ndiStandard[0]+1) * kx * xfact -
                      ndiStandard[1]    * ky * yfact -
                      ndiStandard[2]    * kz * zfact ) ) + 
      cfact * (-cos( (ndiStandard[0]+1) * kx * xfact +
                      ndiStandard[1]    * ky * yfact +
                      ndiStandard[2]    * kz * zfact ) + 
		cos( (ndiStandard[0]+1) * kx * xfact -
                      ndiStandard[1]    * ky * yfact -
                      ndiStandard[2]    * kz * zfact ) );
      
    // RC FFT tests
  
    RFieldSPStan = real(CFieldPPStan);
    CFieldSPStan0h = dcomplex(0.0,0.0);
   
    RFieldSPStan_save = RFieldSPStan;

    // create RC FFT object
    FFT<RCTransform,D,double> rcfft(ndiStandard, ndiStandard0h, compressTemps);

    // set direction names
    rcfft.setDirectionName(+1, "forward");
    rcfft.setDirectionName(-1, "inverse");

    //-------------------------------------------------------------------------

    testmsg << "RC transform using layout with zeroth dim serial ..." << endl;
    timer.start();
    
    rcfft.transform("forward", RFieldSPStan,  CFieldSPStan0h, constInput);
    rcfft.transform("inverse", CFieldSPStan0h, RFieldSPStan, constInput);
    timer.stop();

    diffFieldSPStan = Abs(RFieldSPStan - RFieldSPStan_save);
    realDiff = max(diffFieldSPStan);
    testmsg << "fabs(realDiff) = " << fabs(realDiff) << endl;
    testmsg << "CPU time used = " << timer.cpu_time() << " secs." << endl;
    timer.clear();
    // =================END INTERACTION LOOP=====================================
    return 0;
  }
}
/***************************************************************************
 * $RCSfile: TestFFT.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:02:12 $
 ***************************************************************************/

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

