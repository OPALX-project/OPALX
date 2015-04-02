/**
 * @file   dsmc.cpp
 * @author Andreas Adelmann <andreas.adelmann@psi.ch> and Helge Krueger
 * @date   Tue Jul  21 16:41:03 2005
 * 
 * @brief  This is a program to simulate Touschek scattering
 * 
 * Usage:  $MYMPIRUN -np 4 ./dsmc lat-1.dat [-restart] --processes 4 --commlib mpi
 * 
 */

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include "Ippl.h"

#ifdef TTProf
  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("TTmainTimer");
  static IpplTimings::TimerRef pairListTimer = IpplTimings::getTimer("TTpairListTimer");
  static IpplTimings::TimerRef applyMapTimer = IpplTimings::getTimer("TTapplyMapTimer");
  static IpplTimings::TimerRef statWritingTimer = IpplTimings::getTimer("TTstatWritingTimer");
  static IpplTimings::TimerRef doTouschekTimer = IpplTimings::getTimer("TTdoTouschekTimer");
  static IpplTimings::TimerRef runTouschekTimer = IpplTimings::getTimer("TTrunTouschekTimer");
  static IpplTimings::TimerRef restartWrTimer = IpplTimings::getTimer("TTrestartWrTimer");
  static IpplTimings::TimerRef restartRdTimer = IpplTimings::getTimer("TTrestartRdTimer");
#endif

#include "GTPhysics.h"
#include "GTChargedParticles.hh"
// #include "GTDistribution.hh"
#include "GTConst.hh"
#include "GTConfigure.hh"
#include "GTElemData.hh"

#include "Algorithms/PartData.h"
#include "GTTimer.hh"

#include "Parser/FileStream.h"

#include <fstream>
#include <new>
#include <exception>

#include <vector>
#include <map>
#include <set> 
#include <iostream>
#include <sstream>

#include "GTDistribution.hh"
#include "GTLinearMapper.hh"

#include "FixedAlgebra/FStaticFP.h"
#include "FixedAlgebra/FNormalForm.h"
#include "FixedAlgebra/DragtFinnMap.h"

#include "GTElemData.hh"
#include "GTSbendMap.hh"
#include "GTDriftMap.hh"
#include "GTMultipoleMap.hh"
#include "GTRecquadMap.hh"
#include "GTRFCavityMap.hh"
#include "GTIBS.hh"
#include "GTRadiation.hh"

using namespace std;
using namespace SPECIS;

typedef ChargedParticles<playout_t>::Vector_t Vector3;

vector< ElMap >               ::iterator latticeIt;
vector< ElMap >               ::reverse_iterator latticeRevIt;
vector< FVpsT >               ::iterator taylorIt;
vector< FMatrix<double,6,6> > ::iterator smplMIt;

SimCfgData simCfg;

void applyMap(ChargedParticles<playout_t> *bunch,  FMatrix<double,6,6> M1);

bool   isSymplectic(FMatrix<double,6,6> totest);
double subDets(FMatrix<double,6,6> mat);

FMatrix<double,6,6> makeCompactCombinedFunctionMagnet(double radius, double fieldindex, double angle);
FMatrix<double,6,6> makeCavity(double lambda);
FMatrix<double,6,6> makeMatchMatrix(double dispersion);

void diagnostics(ChargedParticles<playout_t> *bunch);

/** 
 * Print information of the simulation in short and detailed format
 * 
 * @param full if true, a detailed status of the simulation is pronted
 * @param eBunch the container representing the particle bunch
 * @param simCfg all simulation related input parameters exapt lattice and distribution are stored here
 * @param distrData the definition of the distribution
 * @param pd the data on the particle (electrons in that case)
 * @param lattice the lattice to track 
 * @param fn1 where the simulation input is read from
 * @param title title
 * @param perL the lenght of one period of the simulations
 */
void printStatus(bool full, 
		 ChargedParticles<playout_t> *bunch, 
		 SimCfgData &simCfg,
		 DistrData &distrData,
		 PartData &pd,
		 vector< ElMap > &lattice,
		 string fn1, string title,
		 double perL)
{
  Inform m("dsmc ");  
  if (full) {
    m << "---------------------------------------------------------------------------------" << endl;
    m << " T O U S C H E K - T R A C K 6D " << title << endl;
    m << "---------------------------------------------------------------------------------" << endl;

    m << "Lattice input file " << fn1
      << " periodL= " << perL   << " [m] periods " << simCfg.periods 
      << " steps "  << 360/simCfg.Division << " [deg]" << endl;
    
    if(simCfg.doTouschek)
      m << " Touschek scattering enabled " << endl;
    else
      m << " Touschek scattering disabled " << endl;

    m << "Number of processors " << Ippl::getNodes() << " all dimensions are parallel " << endl;
    m << "Base filename for data save " << simCfg.diagfile << endl << endl; 
    
    m << "Particle Data information" << endl;
    m << "q= " << pd.getQ() << "e [Cb] mass= " << pd.getM() << " [eV/c^2]  beta= ";
    m << bunch->getBeta() << " gamma= " << bunch->getGamma() << endl;
    m << "Etot= " << pd.getE() << " [eV] P0= " << pd.getP() << " [eV/c] " << endl << endl;
  
  
    m << "Distribution informations " << endl;
    m << "------------------------- " << endl;      
    m << " Electrons Gauss:  sigx= " << distrData.sigx << " [m]\t sigpx= " <<  distrData.sigpx << " [m/s]\t xpx-corel= " << distrData.corxpx << " [1]" << endl; 
    m << "                 sigy= " << distrData.sigy << " [m]\t sigpy= " <<  distrData.sigpy << " [m/s]\t ypy-corel= " << distrData.corypy << " [1]" << endl; 
    m << "                 sigz= " << distrData.sigt << " [m]\t sigpz= " <<  distrData.sigpt << " [m/s]\t zpz-corel= " << distrData.cortpt << " [1]" << endl; 

    /*    m << "Lattice informations " << endl;
    m << "-------------------- " << endl;
        
    vector< ElMap >::iterator latticeIt;
    unsigned int latItListLen = 0;
    
    for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
      if (((*latticeIt).elemData.elType != STATMARKER) &&
	  ((*latticeIt).elemData.elType != PARTMARKER)) {
	latItListLen++;
	m << (*latticeIt).elemData.convStr();
	m << " entrance at (" << (*latticeIt).elemData.entrance;
	m << " exit at " << (*latticeIt).elemData.entrance+(*latticeIt).elemData.length << "(" 
	  << " a= "  << (*latticeIt).elemData.a 
	  << " b= "  << (*latticeIt).elemData.b 
	  << " L= "  << (*latticeIt).elemData.length 
	  << " kx= " << (*latticeIt).elemData.K0 
	  << " ky= " << (*latticeIt).elemData.K1
	  << " kz= " << (*latticeIt).elemData.K2 << endl << endl;	
      }
      }*/
     
    m.setf(ios::scientific);
    m << " Ne= " << bunch->getTotalNum() << " \t Qe= " << bunch->getQTot() << " [Cb] " << endl << endl;
    
    m << "Parameters from the electron bunches " << endl
      << "------------------------------------ " << endl
      << " Electrons maxX= " << bunch->get_maxExtend() << " [m] " << endl 
      << "           minX= " << bunch->get_origin() << " [m]" << endl;
  }
  else {
    
    
  }
}


/**
 * is Symplectic tests if a given 6 times 6 Matrix is symplectic
 * or not. In case of no, it outputs the determinant as given by
 * the function subDets. 
 * @param totest: a 6 x 6 matrix with double entries
 * @return true if the matrix is symplectic and false if not.
 */
bool isSymplectic(FMatrix<double,6,6> totest) {
  FMatrix<double,6,6> transpose;
  FMatrix<double,6,6> sym;
  FMatrix<double,6,6> test;
  int i,j;
  double sensi = 1e-20;

  for(i=0;i<6;i++)
    for(j=0;j<6;j++) {
      transpose[j][i] = totest[i][j];
      sym[i][j] = 0;
    }
  for(i=0;i<5;i+=2) {				
    sym[i+1][i] = -1;
    sym[i][i+1] = 1;
  }

  test = transpose * sym * totest - sym;

  if(-sensi < subDets(test) < sensi) {
    Inform msg("dsmc ");
    //	msg << test << endl;
    msg << subDets(test) << endl;
    
    return false;
  }

  return true;
}

/**
 * subDets takes a 6 x 6 matrix and computes the determinants
 * of the 2 x 2 submatrixs on the diagonal and adds them
 * used to compute the determinant of matrices with only
 * entries on those.
 * @param mat: a 6 x 6 matrix with double entries
 */
double subDets(FMatrix<double,6,6> mat) {
  int i;
  double result = 0.0;

  for(i=0;i<6;i+=2) {
    result += mat[i][i] * mat[i+1][i+1] - mat[i][i+1] * mat[i+1][i];
  }

  result /= 3;

  return result;
}


/**
 * applymap applies the matrix M1 to the beam of ChargedParticles
 * bunch.
 * @param bunch the bunch of charged particles to apply the matrix to
 * @param M1 the transformation matrix
 */

void applyMap(ChargedParticles<playout_t> *bunch,  FMatrix<double,6,6> M1) {
  FVector<double,6> z;
#ifdef TTProf
      IpplTimings::startTimer(applyMapTimer); 
#endif    
  for(unsigned i=0; i<bunch->getLocalNum(); i++) {
    z[0] = bunch->R[i](0);
    z[1] = bunch->R[i](3); 
    z[2] = bunch->R[i](1);  
    z[3] = bunch->R[i](4); 
    z[4] = bunch->R[i](2); 
    z[5] = bunch->R[i](5); 
    z = M1*z; 
    bunch->R[i](0) = z[0];
    bunch->R[i](3) = z[1];
    bunch->R[i](1) = z[2];
    bunch->R[i](4) = z[3];
    bunch->R[i](2) = z[4];
    bunch->R[i](5) = z[5];
  } 
  bunch->boundp();
#ifdef TTProf
      IpplTimings::stopTimer(applyMapTimer); 
#endif

}

/**
 * makeCompactCombinedFunctionMagnet creates the 6 x 6 matrix for a compact bending matrix....
 * @param radius Radius of the magnet
 * @param fieldindex the fieldindex
 * @param angle the angle of the part of the magnet to consider. 2 pi for one time around
 * @return a 6 x 6 matrix describing the magnet
 */
FMatrix<double,6,6> makeCompactCombinedFunctionMagnet(double radius, double fieldindex, double angle) {
  FMatrix<double,6,6> result;
  int i, j;

  for(i = 0; i < 6; i++)
    for(j = 0; j < 6; j++) {
      result[i][j] = 0;
    }

  double phix, phiy;

  /*
   * Layout taken from the solution to Exercise 3:
   *           "The AURORA compact light source"
   *
   * in "Lattice and Emittance" by Andreas Streun
   *
   * Adapted to 6D with (2.74) in "Computional Methods
   * in Accelerator Physics" by Robert D. Ryne
   */


  phix = angle * sqrt(1 - fieldindex);
  phiy = angle * sqrt(fieldindex);

  //  Inform msg("ttrack ");
  // msg << "angle: " << angle << " phix " << phix << endl;

  result[0][0] = cos(phix);
  result[0][1] = radius / sqrt(1 - fieldindex) * sin(phix);
  result[1][0] = - sqrt(1 - fieldindex) / radius * sin(phix);
  result[1][1] = cos(phix);

  // Not sure if in the right place, because of no delta p / p0
  result[0][5] = - radius / (1 - fieldindex) * (1 - cos(phix));
  result[1][5] = - sin(phix) / ( sqrt(1 - fieldindex));
  
  result[2][2] = cos(phiy);
  result[2][3] = radius / sqrt(fieldindex) * sin(phiy);
  result[3][2] = - sqrt(fieldindex) / radius * sin(phiy);
  result[3][3] = cos(phiy);

  // Still needing Data for the t, pt variables
  result[4][4] = 1;
  result[5][5] = 1;
  result[4][5] = - radius * angle / (1 - fieldindex) + radius/pow(1 - fieldindex,3.0 / 2) * sin(phix);

  result[4][0] = - 1 / sqrt( 1 - fieldindex) * sin(phix);
  result[4][1] = - radius / (1 - fieldindex) * (1 - cos(phix));

  return result;
}

/**
 * make Cavity creates an acceleration cavity
 * @param lambda the value (need better description)
 * @return the matrix for the cavity
 */

FMatrix<double,6,6> makeCavity(double lambda) {
  FMatrix<double,6,6> result;

  int i,j;

  for(i = 0; i < 6; i++)
    for(j = 0; j < 6; j++) {
      if(i == j) result[i][j] = 1;
      else result[i][j] = 0;
    }

  result[5][4] = - lambda;

  return result;
}


/**
 * make match matrix used to fit the beam to the layout
 * @return the match matrix
 */

FMatrix<double,6,6> makeMatchMatrix(double dispersion) {
  FMatrix<double,6,6> result;

  int i,j;

  for(i = 0; i < 6; i++)
    for(j = 0; j < 6; j++) {
      if(i == j) result[i][j] = 1;
      else result[i][j] = 0;
    }

  result[0][5] = dispersion;

  return result;
}

void diagnostics(ChargedParticles<playout_t> *bunch) {
  Inform msg("ttrack ");


}

int main(int argc, char *argv[]) {

  // initialize Ippl
  Ippl ippl(argc, argv);
  Inform msg2all(argv[0], INFORM_ALL_NODES);
  Inform msg(argv[0]);

  MTimer simtimer;
  string dateStr(simtimer.date());
  string timeStr(simtimer.time());

  string restartfile;
  bool restartMode = false;
  int turn = 0;

  ElemData   elemData; 
  DistrData  distrData;

  double gamma,beta,q,P0,mass,couplingConstant;

  vector< ElMap > lattice;

  string startup(argv[1]);
  simCfg.simInputFileName = string(argv[1]);
  FileStream::setEcho(false);
  FileStream *is;

#ifdef TTProf
  IpplTimings::startTimer(mainTimer); 
#endif

  try {
      is = new FileStream(startup);
  } catch (...) {
      is = 0;
      cerr << "No startup file \"" << startup << "\" found." << endl;
  }
  
  if (!FSConfigure(is, &simCfg, &distrData))
    ERRORMSG("Error in configuration section of file \"" << startup << "\" found." << endl);
  
  string simTitle = dateStr + " " + timeStr;
  
  // in classic 3.3 will this will be more difficult
  FTpsT::setGlobalTruncOrder(simCfg.order);
 
  if (restartMode = (argc==3 && (string(argv[2]) == string("-restart"))))
      restartfile = simCfg.datadir+"/"+string(argv[1]) + ".restart.h5";
  
  /*
    Physics::m_e, Physics::m_p rest mass GeV
    Ek: kinetic energy in Gev
  */
  
  gamma = 1.0 + (simCfg.Ek / Physics::m_e);
  beta = sqrt(1-(1/(gamma*gamma)));
  P0 = beta*gamma*Physics::m_e;
  q = -1.0;
  mass = Physics::m_e;
  couplingConstant = 1.0/(4.0*Physics::pi*Physics::epsilon_0);

  PartData pd(q,mass*1.0e9,P0*1.0e9);

  //  double   radius = atof(argv[2]);

  double radius = simCfg.interrad;

  /* ------------------------------------------------------------------*/
  /* Create Lattice                                                    */
  /* ------------------------------------------------------------------*/
  
  double perLength = buildLattice (is,lattice,simCfg,pd);

  ChargedParticles<playout_t> *eBunch = new ChargedParticles<playout_t>(pd,simCfg);

  /* ------------------------------------------------------------------------ */
  /* Create Distributions or load restart file                                */
  /* ------------------------------------------------------------------------ */

  Distribution *distr = new Distribution(&distrData,perLength);
  
  Vector_t xmean (0.0);
  Vector_t xstddev = Vector3(distrData.sigx, distrData.sigy, distrData.sigt,
			    distrData.sigpx,distrData.sigpy,distrData.sigpt);
  Vector_t angle (0.0);
    
  distr->generateGaussianTT(eBunch, simCfg.nEInit, xmean, xstddev, angle);

  bool printStatusFull = true;
  printStatus(printStatusFull,eBunch,simCfg,distrData,pd,lattice,string(argv[1]),simTitle, perLength);

  eBunch->boundp();

  msg2all << "Nloc= " << eBunch->getLocalNum() << endl;

  eBunch->Q = 1.0/eBunch->getTotalNum();
  eBunch->scatterQ();
  eBunch->checkScatteredQ();
  
  /* ------------------------------------------------------------------------- */
  /* Aurora simulation                                                         */
  /* ------------------------------------------------------------------------- */
  
  FMatrix<double,6,6> transfo, cav, tmp;

  unsigned i,j,r;

  transfo = makeCompactCombinedFunctionMagnet(simCfg.radius, simCfg.fieldindex, 2 * Physics::pi/simCfg.Division);
  cav     = makeCavity(simCfg.lambda);

  msg << "Transformation Matrix" << endl;
  msg << transfo << endl;
  msg << (isSymplectic(transfo) ? "symplectic" : "!!! NOT SYMPLECTIC !!!") << endl << endl;

  msg << "Cavity Matrix" << endl;
  msg << cav << endl;
  msg << (isSymplectic(cav) ? "symplectic" : "!!! NOT SYMPLECTIC !!!") << endl << endl;

  tmp = transfo;

  for(i=0; i <179; i++)
    tmp = tmp * transfo;
  tmp = tmp * cav;
  for(i=0; i <180; i++)
    tmp = tmp * transfo;
  
  msg << tmp << endl;

  Ippl::Comm->barrier();	// wait for all nodes to finish
  return 0;
}

/***************************************************************************
 * $RCSfile: pdbtest.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: pdbtest.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
