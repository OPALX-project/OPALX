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
#include "Utilities/Timer.h"

#include <fstream>

#include <complex>
using namespace std;

#define THREED

enum FsT {FFTSolver,Pt2PtSolver,TreeSolver,None};
enum InterPolT {NGP,CIC};
enum BCT {OOO,OOP,PPP,DDD,DDO,DDP};   // OOO == all dim, open BC
enum TestCases {test1};



void writeMemoryHeader(std::ofstream &outputFile)
{

    std::string dateStr("no time");
    std::string timeStr("no time");
    std::string indent("        ");

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    outputFile << "SDDS1" << std::endl;
    outputFile << "&description\n"
               << indent << "text=\"Memory statistics '"
               << "TestFFT" << "' "
               << dateStr << "" << timeStr << "\",\n"
               << indent << "contents=\"stat parameters\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of TestFFT\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"n.a\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=memory,\n"
               << indent << "type=double,\n"
               << indent << "units=" + memory->getUnit() + ",\n"
               << indent << "description=\"2 Total Memory\"\n"
               << "&end\n";

    unsigned int columnStart = 3;

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        outputFile << "&column\n"
                   << indent << "name=processor-" << p << ",\n"
                   << indent << "type=double,\n"
                   << indent << "units=" + memory->getUnit() + ",\n"
                   << indent << "description=\"" << columnStart
                   << " Memory per processor " << p << "\"\n"
                   << "&end\n";
        ++columnStart;
    }

    outputFile << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    outputFile << Ippl::getNodes() << std::endl;
    outputFile << "IPPL test" << " " << "no version" << " git rev. #" << "no version " << std::endl;
    outputFile << "IPPL test" << std::endl;
}

void open_m(std::ofstream& os, const std::string& fileName, const std::ios_base::openmode mode) {
    os.open(fileName.c_str(), mode);
    os.precision(15);
    os.setf(std::ios::scientific, std::ios::floatfield);
}

void writeMemoryData(std::ofstream &os_memData, unsigned int pwi, unsigned int step)
{
    os_memData << step << "\t";     // 1

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    int nProcs = Ippl::getNodes();
    double total = 0.0;
    for (int p = 0; p < nProcs; ++p) {
        total += memory->getMemoryUsage(p);
    }

    os_memData << total << std::setw(pwi) << "\t";

    for (int p = 0; p < nProcs; p++) {
        os_memData << memory->getMemoryUsage(p)  << std::setw(pwi);

        if ( p < nProcs - 1 )
            os_memData << "\t";

    }
    os_memData << std::endl;
}

bool Configure(int argc, char *argv[], InterPolT *interPol, 
	       unsigned int *nx, unsigned int *ny, unsigned int *nz,
	       TestCases *test2Perform,
	       int *serialDim, unsigned int *processes,  unsigned int *nLoop) 
{

  Inform msg("Configure ");
  Inform errmsg("Error ");

  string bc_str;
  string interPol_str;
  string dist_str;

  for (int i=1; i < argc; ++i) {
    string s(argv[i]);
    if (s == "-grid") {
      *nx = atoi(argv[++i]);
      *ny = atoi(argv[++i]);
      *nz = atoi(argv[++i]);
    } else if (s == "-test1") {
      *test2Perform = test1;
    } else if (s == "-NGP") {
      *interPol = NGP;
      msg << "Interploation NGP" << endl;
    } else if (s == "-CIC") {
      *interPol = CIC;
      msg << "Interploation CIC" << endl;
    } else if (s == "-Loop") {
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

  unsigned int processes;
  int serialDim;
  TestCases test2do;
  unsigned int nx,ny,nz;
  unsigned int nLoop;
  InterPolT interPol;

  Configure(argc, argv, &interPol, &nx, &ny, &nz, 
            &test2do, &serialDim, &processes, &nLoop); 
  std::string baseFn = std::string(argv[0])  + 
    std::string("-mx=") + std::to_string(nx) +
    std::string("-my=") + std::to_string(ny) +
    std::string("-mz=") + std::to_string(nz) +
    std::string("-p=") + std::to_string(processes) +
    std::string("-ddec=") + std::to_string(serialDim) ;

  unsigned int pwi = 10;
  std::ios_base::openmode mode_m = std::ios::out;

  std::ofstream os_memData;
  open_m(os_memData, baseFn+std::string(".mem"), mode_m);

  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");           
  IpplTimings::startTimer(mainTimer);                                                    

  static IpplTimings::TimerRef fftTimer = IpplTimings::getTimer("FFT");           
  static IpplTimings::TimerRef ifftTimer = IpplTimings::getTimer("IFFT");           
  static IpplTimings::TimerRef fEvalTimer = IpplTimings::getTimer("fEval");           
  static IpplTimings::TimerRef fInitTimer = IpplTimings::getTimer("init");           

  writeMemoryHeader(os_memData);

  // The preceding cpp definition causes compile-time setting of D:
  const unsigned D=3U;

  testmsg << "Dimensionality: D= " << D << " P= " << processes;
  testmsg << " nx= " << nx << " ny= " << ny << " nz= " << nz << endl;
  
  unsigned ngrid[D];   // grid sizes

  // Used in evaluating correctness of results:
  double realDiff;
  
  // Various counters, constants, etc:
  unsigned int d;
  
  Ippl::Comm->next_tag(IPPL_APP_TAG0);
  double pi = acos(-1.0);
  double twopi = 2.0*pi;
  // Timer:
  Timer timer;
  
  // Layout information:
  unsigned vnodes;             // number of vnodes; input at cmd line
  e_dim_tag allParallel[D];    // Specifies SERIAL, PARALLEL dims
  for (d=0; d<D; d++) 
    allParallel[d] = PARALLEL;

  // Compression of temporaries:
  bool compressTemps;
  vnodes = processes;
  compressTemps = true;
  ngrid[0]=nx;
  ngrid[1]=ny;
  ngrid[2]=nz;

  testmsg << "In-place CC transform using all-parallel layout ..." << endl;

  // ================BEGIN INTERACTION LOOP====================================
 
    //------------------------complex<-->complex-------------------------------
    // Complex test Fields
    // create standard domain
  IpplTimings::startTimer(fInitTimer);
  NDIndex<D> ndiStandard;
  for (d=0; d<D; d++) 
    ndiStandard[d] = Index(ngrid[d]);
  
  // create new domain with axes permuted to match FFT output
  NDIndex<D> ndiPermuted;
  ndiPermuted[0] = ndiStandard[D-1];
  for (d=1; d<D; d++) 
    ndiPermuted[d] = ndiStandard[d-1];

  // create half-size domain for RC transform along zeroth axis
  NDIndex<D> ndiStandard0h = ndiStandard;
  ndiStandard0h[0] = Index(ngrid[0]/2+1);
  // create new domain with axes permuted to match FFT output
  NDIndex<D> ndiPermuted0h;
  ndiPermuted0h[0] = ndiStandard0h[D-1];
  for (d=1; d<D; d++) 
    ndiPermuted0h[d] = ndiStandard0h[d-1];

  // create half-size domain for sine transform along zeroth axis
  // and RC transform along first axis
  NDIndex<D> ndiStandard1h = ndiStandard;
  ndiStandard1h[1] = Index(ngrid[1]/2+1);
  // create new domain with axes permuted to match FFT output
  NDIndex<D> ndiPermuted1h;
  ndiPermuted1h[0] = ndiStandard1h[D-1];
  for (d=1; d<D; d++) 
    ndiPermuted1h[d] = ndiStandard1h[d-1];

  // all parallel layout, standard domain, normal axis order
  FieldLayout<D> layoutPPStan(ndiStandard,allParallel,vnodes);
    
  FieldLayout<D> layoutPPStan0h(ndiStandard0h,allParallel,vnodes);
    
  FieldLayout<D> layoutPPStan1h(ndiStandard1h,allParallel,vnodes);
    
  // create test Fields for complex-to-complex FFT
  BareField<dcomplex,D> CFieldPPStan(layoutPPStan);
  BareField<dcomplex,D> CFieldPPStan_save(layoutPPStan);
  BareField<double,D> diffFieldPPStan(layoutPPStan);
  

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

  CFieldPPStan_save = CFieldPPStan; // Save values for checking later
    
  // CC FFT tests
  // Instantiate complex<-->complex FFT object
  // Transform in all directions
  FFT<CCTransform,D,double> ccfft(ndiStandard, compressTemps);
  IpplTimings::stopTimer(fInitTimer);

  // set direction names
  ccfft.setDirectionName(+1, "forward");
  ccfft.setDirectionName(-1, "inverse");
  
  testmsg << " Start test " << endl;
                                                    
  IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
  memory->sample();
  writeMemoryData(os_memData, pwi, 0);

  for (unsigned int i=0; i<nLoop; i++)  {
    // Test complex<-->complex transform (simple test: forward then inverse transform, see if get back original values.

    IpplTimings::startTimer(ifftTimer);                                                    
    ccfft.transform( "inverse" , CFieldPPStan);
    IpplTimings::stopTimer(ifftTimer);                                                         
    IpplTimings::startTimer(fftTimer);                                                    
    ccfft.transform( "forward" , CFieldPPStan);
    IpplTimings::stopTimer(fftTimer);                                                    


    IpplTimings::startTimer(fEvalTimer);                                                    
    diffFieldPPStan = Abs(CFieldPPStan - CFieldPPStan_save);
    realDiff = max(diffFieldPPStan);

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
    CFieldPPStan_save = CFieldPPStan;
    IpplTimings::stopTimer(fEvalTimer);                
                                    
    testmsg << "FFT->IFFT: CC <-> CC: fabs(realDiff) = " << fabs(realDiff) << endl;

    memory->sample();
    writeMemoryData(os_memData, pwi, i+1);
  }
  testmsg << " Stop test " << endl;
  IpplTimings::stopTimer(mainTimer);
  IpplTimings::print();
  IpplTimings::print(std::string(baseFn+std::string(".timing")));
    
  return 0;
}

/***************************************************************************
 * $RCSfile: TestFFT.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:36 $
 ***************************************************************************/

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

    /*
    ofstream myfile;
    myfile.open("IPPL_FFT", ofstream::app);
    myfile << nx << "\t" << timer.clock_time() << "\t" << timer.clock_time() / nLoop << endl;
    myfile.close();
    */
