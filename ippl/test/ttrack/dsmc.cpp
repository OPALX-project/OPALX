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
#include "GTDistribution.hh"
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
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FLieGenerator.h"
#include "FixedAlgebra/FStaticFP.h"
#include "FixedAlgebra/FNormalForm.h"

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

void diagnostics(ChargedParticles<playout_t> *bunch);
double cross_checkLoss(ChargedParticles<playout_t> *bunch, int id1, int id2, double Phi, double Theta);
double cross_checkLoss(ChargedParticles<playout_t> *bunch, Vector3 p1, Vector3 p2, Vector3 r, double Theta, double Phi);
void cross_run(ChargedParticles<playout_t> *bunch, double length);

FMatrix<double,6,6> calcAccMatrix(vector< ElMap > &lattice);

/**
 *
 *
 *
 **/
FMatrix<double,6,6> calcAccMatrix(vector< ElMap > &lattice) {
  
  FMatrix<double,6,6>  macc;
  
  for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++) {
      if(i == j) 
	macc[i][j] = 1;
      else 
	macc[i][j] = 0;
    }

  //  Inform msg("ttt ");
  
  for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
    if (((*latticeIt).elemData.elType != STATMARKER) &&
	((*latticeIt).elemData.elType != PARTMARKER)) {
      for(int i = 0; i < (*latticeIt).elemData.slice; i++)
	macc = macc * (*latticeIt).M;
    }
    //    msg << macc << endl;
  }
  return macc;
}




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
    m << " T O U S C H E K - T R A C K " << title << endl;
    m << "---------------------------------------------------------------------------------" << endl;

    m << "Lattice input file " << fn1
      << " periodL= " << perL   << " [m] periods " << simCfg.periods 
      << " steps "  << 360/simCfg.Division << " [deg]" << endl;
    
    if(simCfg.doTouschek)
      m << " Touschek scattering enabled " << endl;
    else
      m << " Touschek scattering disabled " << endl;

    m << " Interaction radius " << simCfg.interrad << endl;
    m << " MeshMultFactor " << simCfg.meshmultfactor << endl;

    m << "Number of processors " << Ippl::getNodes() << " serial dimension is z " << endl;
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


    double length = 0;

    /*    
    for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
      if (((*latticeIt).elemData.elType != STATMARKER) &&
	  ((*latticeIt).elemData.elType != PARTMARKER)) {
	
	length += (*latticeIt).elemData.length;
	

        if ((*latticeIt).elemData.elType != RFCAVITY) {
	  m << (*latticeIt).elemData.convStr();
	  m << " entrance at (" << (*latticeIt).elemData.entrance;
	  m << " exit at " << (*latticeIt).elemData.entrance+(*latticeIt).elemData.length << "(" 
	    << " a= "  << (*latticeIt).elemData.a 
	    << " b= "  << (*latticeIt).elemData.b 
	    << " L= "  << (*latticeIt).elemData.length 
	    << " k0= " << (*latticeIt).elemData.K0 
	    << " k1= " << (*latticeIt).elemData.K1
	    << " k2= " << (*latticeIt).elemData.K2 << endl 
	    << "M= "  << endl  << (*latticeIt).M << endl << endl;	
	}
	else {
	  m << (*latticeIt).elemData.convStr();
	  m << " entrance at (" << (*latticeIt).elemData.entrance;
	  m << " exit at " << (*latticeIt).elemData.entrance+(*latticeIt).elemData.length << "(" 
	    << " a= "  << (*latticeIt).elemData.a 
	    << " b= "  << (*latticeIt).elemData.b 
	    << " L= "  << (*latticeIt).elemData.length 
	    << " V= " << (*latticeIt).elemData.v 
	    << " PHI= " << (*latticeIt).elemData.phi
	    << " FREQ= " << (*latticeIt).elemData.freq << endl 
	    << "M= " << endl   << (*latticeIt).M << endl << endl;	
	}
      }
      }*/

   FMatrix<double,6,6> full = calcAccMatrix(lattice);

    m << "macc= " << endl << full << endl;
    m << "betax= " << sqrt(- full[0][1] / full[1][0]) << "\tbetay =" << sqrt(- full[2][3] / full[3][2]) << "\tbetaz =" << sqrt(abs(- full[4][5] / full[5][4])) << endl;


    m << "lenght = " << length << endl;
         
    m.setf(ios::scientific);
    m << " Ne= " << bunch->getTotalNum() << " \t Qe= " << bunch->getQTot() << " [Cb] " << endl << endl;
    
    m << "Parameters from the electron bunches " << endl
      << "------------------------------------ " << endl
      << " Electrons maxX= " << bunch->get_maxExtend() << " [m]\t minX= " 
      <<  bunch->get_origin() << " [m]\t" << endl;
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

/*
 * Using the scattering angle, monte carlo integration approach, this function
 * returns gamma if particle 1 and 2 are lost under the scattering angles Phi und Theta
 * @param bunch the bunch of particles
 * @param p1 momentum of the first particle to collide
 * @param p2 momentum of the second particle to collide
 * @param r place of ths collision
 * @param Theta polar scattering angle
 * @param Phi azimutal scattering angle
 * @return the relative speed * crosssection transformed to bunch system
 */
double cross_checkLoss(ChargedParticles<playout_t> *bunch, Vector3 p1, Vector3 p2, Vector3 r, double Theta, double Phi) {
  // Momenta in the Bunch System
  p1(2) += 1;
  p1 *= bunch->getGamma() * Physics::EMASS * bunch->getBeta() * Physics::c;
  p2(2) += 1;
  p2 *= bunch->getGamma() * Physics::EMASS * bunch->getBeta() * Physics::c;

  // Lorentz Transform
  p1(2) = bunch->getGamma() * p1(2) - bunch->getBeta() * bunch->getGamma() * sqrt(dot(p1,p1) + pow(Physics::EMASS * Physics::c, 2));
  p2(2) = bunch->getGamma() * p2(2) - bunch->getBeta() * bunch->getGamma() * sqrt(dot(p2,p2) + pow(Physics::EMASS * Physics::c, 2));

  double v2 = Physics::c * sqrt(dot(p1,p1)) / sqrt(pow(Physics::EMASS * Physics::c,2) + dot(p1,p1));

  // Energies divided by c.
  double Ec1 = sqrt(dot(p1,p1) + pow(Physics::EMASS * Physics::c, 2));
  double Ec2 = sqrt(dot(p2,p2) + pow(Physics::EMASS * Physics::c, 2));

  Vector3 beta = (p1 + p2) / (Ec1 + Ec2);
  double gamma = 1 / sqrt(1 - dot(beta,beta));

  // Transformation to the center of mass system
  Vector3 q1 = p1 + beta * gamma * ( gamma / (gamma + 1) * dot(beta,p1) - Ec1);
  Vector3 q2 = p2 + beta * gamma * ( gamma / (gamma + 1) * dot(beta,p2) - Ec2);

  // Compute more stuff, angles and length ...
  double q = sqrt(dot(q1,q1));
  double theta = acos(q1(2) / q);
  double phi = Physics::pi / 2 + atan(q1(1) / q1(0));
  if(q1(0) < 0) phi += Physics::pi;

  // Do the scattering (here a rotation)
  Vector3 qafter, pafter;
  qafter(0) = q * (cos(phi) * sin(Theta) * cos(Phi) - cos(theta) * sin(phi) * sin(Theta) * sin(Phi) + sin(theta) * sin(phi) * cos(Theta));
  qafter(1) = q * (sin(phi) * sin(Theta) * cos(Phi) + cos(theta) * cos(phi) * sin(Theta) * sin(Phi) - sin(theta) * cos(phi) * cos(Theta));
  qafter(2) = q * (sin(theta) * sin(Theta) * sin(Phi) + cos(theta) * cos(Theta));




  // Transform back to the bunch system
  pafter= qafter - beta * gamma * ( gamma / (gamma + 1) * dot(beta * (-1),qafter) - sqrt(dot(qafter,qafter) + pow(Physics::EMASS * Physics::c, 2)));

  // Back to lab system
  pafter(2) = bunch->getGamma() * pafter(2) + bunch->getBeta() * bunch->getGamma() * sqrt(dot(pafter,pafter) + pow(Physics::EMASS * Physics::c, 2));
  pafter /= bunch->getGamma() * Physics::EMASS * bunch->getBeta() * Physics::c;
  pafter(2) -= 1;

  // Compute crosssection
  double v = Physics::c * q / sqrt(pow(Physics::EMASS * Physics::c,2) + q*q);
  double x = pow(Physics::c / v,2);
  double sigma = pow(Physics::r_e,2)/4 *( 1 - pow(v/Physics::c,2)) * (pow(x + 1,2)*(4/pow(sin(Theta),4) - 3/pow(sin(Theta),3)) + 1 + 4/pow(sin(Theta),2));

  // Momenta Acceptance
  if(abs(pafter(2)) > simCfg.accps) return v2 * sigma * gamma;

  return 0;
}

/**
 * loops to all particles in the bunch and scatters them with another particle at random
 * @param bunch the bunch of particles
 * @length length of the linear optics element
 */
void cross_run(ChargedParticles<playout_t> *bunch, double length) {
  unsigned second; double Theta, Phi;
  double sigma, v, x, rho1 = 0, rho2 = 0;
  double deltau = 0, deltauG = 0;
  double deltauana = 0, deltauanaG = 0;
  double deltausemiana = 0, deltausemianaG = 0;
  unsigned scat = 0;
  double rho1ana = 0, rho2ana = 0;
  double sigmamax = 0;
  unsigned i, iG;
  Vector3 r, p1, p2;
  double a = simCfg.sigmawidth;

  bunch->calcBeamParameters();

  for(i = 0; i < bunch->getLocalNum(); i ++) {
    // Chose scattering params
    second = (unsigned)floor(bunch->getLocalNum() * IpplRandom());
    Theta = (Physics::pi / 2 - Physics::pi/2 * simCfg.thetamin) * IpplRandom() + Physics::pi/2 * simCfg.thetamin;
    Phi = Physics::pi * IpplRandom();

    // 2 Random Particle Approach
    if(v = cross_checkLoss(bunch,bunch->P[i], bunch->P[second], bunch->R[i], Theta, Phi)) {
      deltau += 2 * v * sin(Theta) *  bunch->getGamma() * length;
    }

    // Semi Analytic approach
    p2 = Vector3(a * sqrt(bunch->get_emit()(0) / bunch->get_csbeta()(0)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(1) / bunch->get_csbeta()(1)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(2) / bunch->get_csbeta()(2)) * 2 * (IpplRandom() - 0.5));

    for(unsigned j=0;j<2;j++) {
      p1(j) -= bunch->R[i](j) * bunch->get_csalpha()(j)/bunch->get_csbeta()(j);
      p2(j) -= bunch->R[i](j) * bunch->get_csalpha()(j)/bunch->get_csbeta()(j);
    }
    
    if(v = cross_checkLoss(bunch,bunch->P[i], Vector3(p2(0),p2(1),p2(2) / bunch->getGamma()),bunch->R[i], Theta, Phi)) {
      rho1ana = bunch->getLocalDensityAnalytic(bunch->get_R_beam(i), p2);
      deltausemiana += 2 * v * sin(Theta) * rho1ana * length;
      scat ++;
    }

    // Analytic approach
    r  = Vector3(a * sqrt(bunch->get_emit()(0) * bunch->get_csbeta()(0)) * 2 * (IpplRandom() - 0.5), 
	 	 a * sqrt(bunch->get_emit()(1) * bunch->get_csbeta()(1)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(2) * bunch->get_csbeta()(2)) * 2 * (IpplRandom() - 0.5));
    p1 = Vector3(a * sqrt(bunch->get_emit()(0) / bunch->get_csbeta()(0)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(1) / bunch->get_csbeta()(1)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(2) / bunch->get_csbeta()(2)) * 2 * (IpplRandom() - 0.5));  
    p2 = Vector3(a * sqrt(bunch->get_emit()(0) / bunch->get_csbeta()(0)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(1) / bunch->get_csbeta()(1)) * 2 * (IpplRandom() - 0.5), 
		 a * sqrt(bunch->get_emit()(2) / bunch->get_csbeta()(2)) * 2 * (IpplRandom() - 0.5));

    for(unsigned j=0;j<2;j++) {
      p1(j) -= r(j) * bunch->get_csalpha()(j)/bunch->get_csbeta()(j);
      p2(j) -= r(j) * bunch->get_csalpha()(j)/bunch->get_csbeta()(j);
    }

    if(v = cross_checkLoss(bunch, Vector3(p1(0),p1(1),p1(2) / bunch->getGamma()), Vector3(p2(0),p2(1),p2(2) / bunch->getGamma()),r,Theta, Phi)) {
      // Get the densities:
      rho1ana = bunch->getLocalDensityAnalytic(r, p2);
      rho2ana = bunch->getLocalDensityAnalytic(r, p1);
      // Add up
      deltauana += 2 * v * sin(Theta) * rho1ana * rho2ana * length;
    }
  }

  reduce(deltau,deltauG,OpAddAssign());
  reduce(deltauana,deltauanaG,OpAddAssign());
  reduce(deltausemiana,deltausemianaG,OpAddAssign());
  reduce(i,iG,OpAddAssign());

  bunch->add_delta_u(deltauG);
  bunch->add_delta_u_analytic(deltauanaG);
  bunch->add_delta_u_semianalytic(deltausemianaG);
  bunch->add_events(iG);

  // Compute the part of the volume

  // Analytic Approach
  double dV = pow(2 * a, 9);
  for(unsigned i=0; i<3; i++) 
    dV *= sqrt(bunch->get_emit()(i) * bunch->get_csbeta()(i)) * pow(sqrt(bunch->get_emit()(i)/bunch->get_csbeta()(i)),2);
  dV *= pow(Physics::pi,2)/2;
  bunch->add_delta_v(dV);

  // Semi Analytic Approach
  double dV2 = pow(2 * a, 3);
  for(unsigned i=0; i<3; i++) 
    dV2 *= sqrt(bunch->get_emit()(i)/bunch->get_csbeta()(i));
  dV2 *= pow(Physics::pi,2)/2;
  bunch->add_delta_v2(dV2);
}


/**
 * doSynchrotron routine that takes care of simulating the effects of synchrotron
 * radiation by reducing the energy of the particles.
 * NOT IMPLEMENTED LIKE IT SHOULD
 * needs to be randomized to take care of photon and quantum stuff
 */
void doSynchrotron(ChargedParticles<playout_t> *bunch) {
  
  double loss = sqrt((pow(Physics::q_e,2) * pow(bunch->getGamma(),2) * 2.909e-11)/(6 * Physics::pi * Physics:: c * 0.25 * 9.109e-31));

  for(unsigned i = 0; i < bunch->getLocalNum(); i ++) {
    bunch->P[i](2) = bunch->P[i](2) - loss;
  }
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
    z[1] = bunch->P[i](0); 
    z[2] = bunch->R[i](1);  
    z[3] = bunch->P[i](1); 
    z[4] = bunch->R[i](2); 
    z[5] = bunch->P[i](2); 
    z = M1*z; 
    bunch->R[i](0) = z[0];
    bunch->P[i](0) = z[1];
    bunch->R[i](1) = z[2];
    bunch->P[i](1) = z[3];
    bunch->R[i](2) = z[4];
    bunch->P[i](2)  = z[5];
  } 
  bunch->boundp();
#ifdef TTProf
      IpplTimings::stopTimer(applyMapTimer); 
#endif

}

void skipLines(unsigned int nL, ifstream& infile)
{
    char ch = infile.get();
    for (int i=0;i<nL;i++) {
        while (ch != '\n')
            ch = infile.get();
        ch = infile.get();
    }
}

/**
 * Prints out the total Momentum and Total Angular Momentum
 * @param bunch the particle bunch
 */
void diagnostics(ChargedParticles<playout_t> *bunch) {
  Inform msg("ttrack ");


  //Calculate total momentum and angular momentum
   Vector3 TotalP_m = sum(bunch->P);
   Vector3  TotalL_m = sum(cross(bunch->R,bunch->P));

   msg << "TotalP " << TotalP_m << endl;
   msg << "TotalL " << TotalL_m << endl;
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

  /* ------------------------------------------------------------------------- */
  /* Aurora simulation                                                         */
  /* ------------------------------------------------------------------------- */

  FMatrix<double,6,6> transfo, cav, tmp;

  unsigned i,j,r;

  /* ------------------------------------------------------------------*/
  /* Create Lattice                                                    */
  /* ------------------------------------------------------------------*/
  
  double perLength = buildLattice (is,lattice,simCfg,pd);

  int parallelDim = 2;
  int meshMultFactor = simCfg.meshmultfactor;
  ChargedParticles<playout_t> *eBunch = new ChargedParticles<playout_t>(radius,pd,simCfg,parallelDim,meshMultFactor);

  /* ------------------------------------------------------------------------ */
  /* Create Distributions or load restart file                                */
  /* ------------------------------------------------------------------------ */

  if(restartMode) {
#ifdef TTProf
      IpplTimings::startTimer(restartRdTimer);
#endif
      turn = eBunch->readRestartInfo(simCfg.datadir+"/"+string(argv[1])+".restart");
      msg << "sum(R)= " << sum(eBunch->R) << " sum(P)= " << sum(eBunch->P) << endl;
#ifdef TTProf
      IpplTimings::stopTimer(restartRdTimer);
#endif
  }
  else {
    Distribution *distr = new Distribution(&distrData,perLength);
    
    Vector3 xmean (0.0);
    Vector3 pmean (0.0);
    Vector3 xstddev = Vector3(distrData.sigx, distrData.sigy, distrData.sigt);
    Vector3 pstddev = Vector3(distrData.sigpx,distrData.sigpy,distrData.sigpt);
    Vector3 angle (0.0);
    
    switch(simCfg.distributiontype) {
    case 1:
      distr->generateGaussianFlatTop(eBunch, simCfg.nEInit, xmean, pmean, xstddev, pstddev, angle);
      break;
    case 2:
      distr->generateGaussianFlatTopHam(eBunch, simCfg.nEInit, xmean, pmean, xstddev, pstddev, angle);
      break;
    case 3:
      distr->generateGaussianTT(eBunch, simCfg.nEInit, xmean, pmean, xstddev, pstddev, angle);
      break;
    }

    for(unsigned i = 0; i < eBunch->getLocalNum(); i++) 
	if (eBunch->ID[i] < 10) {
	    msg2all << "r= " << eBunch->R[i] << " p= " << eBunch->P[i] << endl; 
	}



    /*   for(unsigned i = 0; i < eBunch->getLocalNum(); i++) {
      eBunch->R[i](0) += simCfg.dispersion * cos(IpplRandom() *2 * Physics::pi * sqrt(1 - simCfg.fieldindex)) * eBunch->P[i](2);
      eBunch->P[i](0) -= simCfg.dispersion * sin(IpplRandom() *2 * Physics::pi * sqrt(1 - simCfg.fieldindex)) * sqrt(1 - simCfg.fieldindex) / simCfg.radius * eBunch->P[i](2);
      }*/
    
    bool printStatusFull = true;
    printStatus(printStatusFull,eBunch,simCfg,distrData,pd,lattice,string(argv[1]),simTitle, perLength);
  
  }
  
  /* ------------------------------------------------------------------------ */
  /* Save initial conditions                                                  */   
  /* ------------------------------------------------------------------------ */
  
  if(simCfg.doStatistics) {
    eBunch->openFiles((simCfg.datadir+"/"+simCfg.diagfile),simCfg.title,restartMode);
    eBunch->writeStatistics();
  }
  if(simCfg.doPhaseSpace)
    eBunch->writePhaseSpaceSDDS((simCfg.datadir+"/"+simCfg.diagfile));
  if(simCfg.doInbalance)
    eBunch->getInbalance(simCfg.datadir+"/"+simCfg.diagfile+".inbal", true);

  if(simCfg.doLostPlace)
    eBunch->init_lose_particle(simCfg.datadir+"/"+simCfg.diagfile+".loss.dat");

  /* ------------------------------------------------------------------------ */
  /* Initial Mesh Configuration                                               */   
  /* ------------------------------------------------------------------------ */

  eBunch->boundp();

  eBunch->Q = 1.0/eBunch->getTotalNum();
  
  eBunch->scatterQ();
  eBunch->checkScatteredQ();

  /* ------------------------------------------------------------------------ */
  /* Simulation loop                                                          */   
  /* ------------------------------------------------------------------------ */

  Timer timer;
  timer.clear();
  timer.start();

  double n = 0;
  
  for (unsigned l=turn; l < simCfg.periods ; ++l) {

    for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
      if (((*latticeIt).elemData.elType != STATMARKER) &&
	  ((*latticeIt).elemData.elType != PARTMARKER)) {
	for(int i = 0; i < (*latticeIt).elemData.slice; i++) {
	  applyMap(eBunch, (*latticeIt).M);  
	  
	  if(simCfg.doStatistics)
	    eBunch->writeStatistics();
	  eBunch->set_spos(eBunch->get_spos() +(*latticeIt).elemData.length /(*latticeIt).elemData.slice );
	  if(simCfg.doInbalance)
	    eBunch->getInbalance(simCfg.datadir+"/"+simCfg.diagfile+".inbal", false);
	  MTimer simtimer;
	  
	  if(simCfg.doCrosssection && (*latticeIt).elemData.length) { // runs the crosssection code
	    eBunch->scatterQ();
	    cross_run(eBunch, (*latticeIt).elemData.length/(*latticeIt).elemData.slice);
	    n += eBunch->getTotalNum();
	    double circ = eBunch->get_spos();
//	    msg << "n " << n << " circ " << circ << endl;
	    msg.setf(ios::scientific);
            msg << "alpha mesh: " << pow(Physics::pi,2) / 2 / n / eBunch->getGamma() * eBunch->get_delta_u(); 
	    msg << " analytic: " << eBunch->get_delta_v()/ n / circ / eBunch->getGamma()  * eBunch->get_delta_u_analytic();
	    msg << " semianalytic: " << eBunch->get_delta_v2()/ circ /n / eBunch->getGamma()  * eBunch->get_delta_u_semianalytic() << endl;
	  }
	}
      }
    }
    msg << "Turn " << l+1 << " @ " << r*360/simCfg.Division  << " [deg] " << simtimer.time() << endl;
	
    
    if(simCfg.doPhaseSpace)
      eBunch->writePhaseSpaceSDDS((simCfg.datadir+"/"+simCfg.diagfile));
    if(simCfg.doRestart){
#ifdef TTProf
      IpplTimings::startTimer(restartWrTimer);
#endif
      eBunch->writeRestartInfo((simCfg.datadir+"/"+string(argv[1])+".restart.h5"), l + 1);
      msg << "sum(R)= " << sum(eBunch->R) << " sum(P)= " << sum(eBunch->P) << endl;
#ifdef TTProf
      IpplTimings::stopTimer(restartWrTimer);
#endif
    }
    
  }
  timer.stop();
  

  // Final computation of alpha

  msg.setf(ios::scientific);
  msg << "deltav " << eBunch->get_delta_v() << endl;
  msg.setf(ios::scientific);
  msg << "deltau " << eBunch->get_delta_u() << endl;
  msg << "deltau analztic " << eBunch->get_delta_u_analytic() << endl;

  //  double n = eBunch->getTotalNum() * simCfg.Division * simCfg.periods;
  double circ = Physics::pi * simCfg.periods;
  
  msg << "n: " << n << " events " << eBunch->get_events() << endl;

  msg << "alpha " << eBunch->get_delta_v()/n / eBunch->getGamma() / circ * eBunch->get_delta_u() << endl;

  msg << "alpha analytic " << eBunch->get_delta_v()/ n / eBunch->getGamma() / circ * eBunch->get_delta_u_analytic() << endl;


  msg << eBunch->get_csalpha() << endl;
  msg << eBunch->get_csbeta() << endl;
  msg << eBunch->get_csgamma() << endl;
  msg << eBunch->get_emit() << endl;


  /*
    Some statistice: partilce push per uS and turn etc.

  */
  unsigned numOfMappsApplied = simCfg.periods * (1 + (unsigned)simCfg.Division*2);
  unsigned partPushes = numOfMappsApplied*eBunch->getTotalNum();
  double partPushesPerUs = partPushes/(timer.cpu_time()*1e6);

  msg << "Using " << Ippl::getNodes() << " nodes we make " << partPushesPerUs << " particle pushes per 1e-6 seconds (cpu time) "
      << " we need for one turn " << timer.cpu_time()/simCfg.periods << " seconds" << endl; 
  //  eBunch->writePhaseSpaceSDDS((simCfg.datadir+"/"+simCfg.diagfile));

  eBunch->closeFiles();

#ifdef TTProf
  IpplTimings::stopTimer(mainTimer); 
#endif

#ifdef TTProf
  string name = simCfg.datadir+"/"+string(argv[0])+".timing";
  msg << "Process " << Ippl::myNode() << " writeout timing to file: " << name << endl;
  IpplTimings::print(name); // use  IpplTimings::print() to print on stdout
#endif    
  Ippl::Comm->barrier();	// wait for all nodes to finish
  return 0;
}

/***************************************************************************
 * $RCSfile: pdbtest.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: pdbtest.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
