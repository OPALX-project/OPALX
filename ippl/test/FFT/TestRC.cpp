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
 ***************************************************************************/

#include "Ippl.h"
#include <typeinfo>

#ifdef IPPL_USE_STANDARD_HEADERS
#include <complex>
using namespace std;
#else
#include <complex.h>
#endif


bool Configure(int argc, char *argv[],
	       unsigned int *nx, unsigned int *ny, unsigned int *nz,
	       int *serialDim, unsigned int *processes,  unsigned int *nLoop) 
{

  Inform msg("Configure ");
  Inform errmsg("Error ");

  /*
    string bc_str;
    string dist_str;
    *nx = 8; 
    *ny = 8;
    *nz = 16;
    *nLoop = 1; 
    *serialDim = 0;

    if (*serialDim == 0)
    msg << "Serial dimension is x" << endl;
    else if (*serialDim == 1)
    msg << "Serial dimension is y" << endl;
    else if (*serialDim == 2)
    msg << "Serial dimension is z" << endl;
    else {
    *serialDim = 0;
    msg << "Serial dimension is x" << endl;
    }
    *processes = Ippl::getNodes();
    */
  for (int i=1; i < argc; ++i) {
    string s(argv[i]);
    if (s == "-grid") {
      *nx = atoi(argv[++i]);
      *ny = atoi(argv[++i]);
      *nz = atoi(argv[++i]);
    }   else if (s == "-Loop") {
      *nLoop = atoi(argv[++i]);
    } else if (s == "-Decomp") {
      *serialDim = atoi(argv[++i]);
    } 
    else {
      errmsg << "Illegal format for or unknown option '" << s.c_str() << "'.";
      errmsg << endl;
    }
  }
  if (*serialDim == 0)
    msg << "Serial dimension is x" << endl;
  else if (*serialDim == 1)
    msg << "Serial dimension is y" << endl;
  else if (*serialDim == 2)
    msg << "Serial dimension is z" << endl;
  else {
    msg << "All parallel" << endl;
    *serialDim = -1;
  }

  *processes = Ippl::getNodes();


  return true;
}

int main(int argc, char *argv[])
{
  
  Ippl ippl(argc,argv);
  Inform testmsg(NULL,0);

  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
  static IpplTimings::TimerRef fftTimer = IpplTimings::getTimer("fftTimer");

  IpplTimings::startTimer(mainTimer);
  
  const unsigned D=3U;
  bool compressTemps = false;
  bool constInput    = true;  // preserve input field in two-field transform


  unsigned int processes;
  int serialDim;
  unsigned int nx,ny,nz;
  unsigned int nLoop;

  Configure(argc, argv, &nx, &ny, &nz, &serialDim, &processes, &nLoop); 

  int vnodes = processes;
  unsigned ngrid[D];   // grid sizes
  ngrid[0] = nx;
  ngrid[1] = ny;
  ngrid[2] = nz;

  // Used in evaluating correctness of results:
  double realDiff;
   
  // Various counters, constants, etc:
  
  double pi = acos(-1.0);
  double twopi = 2.0*pi;

  // Timer:
  Timer timer;

  e_dim_tag allParallel[D];    // Specifies SERIAL, PARALLEL dims
  for (unsigned int d=0; d<D; d++) 
    allParallel[d] = PARALLEL;

  e_dim_tag serialParallel[D]; // Specifies SERIAL, PARALLEL dims
  for (unsigned int d=0; d<D; d++) 
    serialParallel[d] = PARALLEL;
  serialParallel[serialDim] = SERIAL;

  // create standard domain
  NDIndex<D> ndiStandard;
  for (unsigned int d=0; d<D; d++) 
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
  
  INFOMSG("RFieldSPStan   layout= " << layoutSPStan << endl;);
  INFOMSG("CFieldSPStan0h layout= " << layoutSPStan0h << endl;);
  
  // Rather more complete test functions (sine or cosine mode):
  dcomplex sfact(1.0,0.0);      // (1,0) for sine mode; (0,0) for cosine mode
  dcomplex cfact(0.0,0.0);      // (0,0) for sine mode; (1,0) for cosine mode
  
  double xfact, kx, yfact, ky, zfact, kz;
  xfact = pi/(ngrid[0] + 1.0);
  yfact = 2.0*twopi/(ngrid[1]);
  zfact = 2.0*twopi/(ngrid[2]);
  kx = 1.0; ky = 2.0; kz = 32.0; // wavenumbers

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
  
  cout << "TYPEINFO:" << endl;
  cout << typeid(RFieldSPStan[0][0][0]).name() << endl;
  cout << typeid(CFieldPPStan[0][0][0]).name() << endl;

  // RC FFT tests  

  //RFieldSPStan = real(CFieldPPStan);
  for(int x = ndiStandard[0].first(); x <= ndiStandard[0].last(); x++) {
    for(int y = ndiStandard[1].first(); y <= ndiStandard[1].last(); y++) {
      for(int z = ndiStandard[2].first(); z <= ndiStandard[2].last(); z++) {
	RFieldSPStan[x][y][z] = real(CFieldPPStan[x][y][z].get());
      }
    }
  }
 
  CFieldSPStan0h = dcomplex(0.0,0.0);

  /*
    Inform fo1(NULL,"realField.dat",Inform::OVERWRITE);
   
    for(int x = ndiStandard[0].first(); x <= ndiStandard[0].last(); x++) {
    for(int y = ndiStandard[1].first(); y <= ndiStandard[1].last(); y++) {
    for(int z = ndiStandard[2].first(); z <= ndiStandard[2].last(); z++) {
    fo1 << x << " " << y << " " << z << " " <<  RFieldSPStan[x][y][z].get() << endl;
    }
    }
    }
  */

  // create RC FFT object
  FFT<RCTransform,D,double> rcfft(ndiStandard, ndiStandard0h, compressTemps);

  // set direction names
  rcfft.setDirectionName(+1, "forward");
  rcfft.setDirectionName(-1, "inverse");

  Inform fo2(NULL,"FFTrealField.dat",Inform::OVERWRITE);
  
  testmsg << "RC transform using layout with zeroth dim serial ..." << endl;


  rcfft.transform("forward", RFieldSPStan,  CFieldSPStan0h, constInput);
  rcfft.transform("inverse", CFieldSPStan0h, RFieldSPStan, constInput);

  for (unsigned i=0; i<nLoop; i++) {
    RFieldSPStan_save = RFieldSPStan;

    IpplTimings::startTimer(fftTimer);
    rcfft.transform("forward", RFieldSPStan,  CFieldSPStan0h, constInput);
    /*
      for(int x = ndiStandard0h[0].first(); x <= ndiStandard0h[0].last(); x++) {
      for(int y = ndiStandard0h[1].first(); y <= ndiStandard0h[1].last(); y++) {
      for(int z = ndiStandard0h[2].first(); z <= ndiStandard0h[2].last(); z++) {
      fo2 << x << " " << y << " " << z << " " <<  real(CFieldSPStan0h[x][y][z].get()) << " " << imag(CFieldSPStan0h[x][y][z].get()) << endl;
      }
      }
      }   
    */

    rcfft.transform("inverse", CFieldSPStan0h, RFieldSPStan, constInput);
    IpplTimings::stopTimer(fftTimer);

    //double total_time = 0;
    //total_time+= timer.cpu_time();
    /*
      Inform fo2(NULL,"FFTrealResult.dat",Inform::OVERWRITE);
      for(int x = ndiStandard[0].first(); x <= ndiStandard[0].last(); x++) {
      for(int y = ndiStandard[1].first(); y <= ndiStandard[1].last(); y++) {
      for(int z = ndiStandard[2].first(); z <= ndiStandard[2].last(); z++) {
      fo2 << x << " " << y << " " << z << " " <<  RFieldSPStan[x][y][z].get() << endl;
      }
      }
      }
    */

  
    diffFieldSPStan = Abs(RFieldSPStan - RFieldSPStan_save);
    realDiff = max(diffFieldSPStan);
    testmsg << "fabs(realDiff) = " << fabs(realDiff) << endl;
  }

  IpplTimings::stopTimer(mainTimer);
  IpplTimings::print();
  IpplTimings::print(std::string("TestRC.timing"));

  return 0;
}
/***************************************************************************
 * $RCSfile: TestRC.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:36 $
 ***************************************************************************/

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

