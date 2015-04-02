#ifndef CfgData_HH
#define CfgData_HH
// ------------------------------------------------------------------------
// $RCSfile: GTElemData.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTElemData.hh,v 1.5 2004/04/28 08:07:36 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.5 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
//
// $Date: 2004/04/28 08:07:36 $
// $Author: adelmann $
// $Log: GTElemData.hh,v $
// Revision 1.5  2004/04/28 08:07:36  adelmann
// Add MAPS #define
//
// Revision 1.4  2004/04/06 13:23:26  adelmann
// *** empty log message ***
//
// Revision 1.3  2003/05/02 13:58:23  adelmann
// First electron cloud simulation 3d without SC and probable still with the wrong
// drive beam integration scheme. Otherwise :=) the program runs
// sureprisingly stable and the first inspection of the data looks promissing
//
// Revision 1.2  2003/04/17 14:21:01  adelmann
// *** empty log message ***
//
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------


#ifdef MAPS
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FLieGenerator.h"
#include "FixedAlgebra/DragtFinnMap.h"
typedef FVps<double,2*DIMENSION> FVpsT;
typedef FTps<double,2*DIMENSION> FTpsT;
#endif

#include <string>
#include <iostream>
#include <vector>
using namespace std;

#include "GTConst.hh"

/*
  Need a better geometry description and a 
  space charge scaling factor
*/

class ElemData
{
public:
  ElemData():
    K0(0.0),
    K0S(0.0),
    K1(0.0),
    K1S(0.0),
    K2(0.0),
    K2S(0.0),
    K3(0.0),
    K3S(0.0),
    length(0.0),
    slice(1.0),
    alpha(0.0),
    beta(0.0),
    fieldindex(0.0),
    lambda(0.0),
    e1(0.0),e2(0.0),angle(0.0),fint(0.0),hgap(0.0),
    a(0.0),b(0.0),
    order(0),
    r1(0.0),  
    r2(0.0),
    v(0.0),
    phi(0.0),
    freq(0.0),
    entrance(0.0),
    isMisaligned(false),
    missX(0.0),
    missY(0.0),
    elemNo(-1)
  {
  }

  void print(ostream &os) const
  {
    os << "Type= ";
    switch (elType) {
    case DRIFT: os << "Drift "; break;
    case QUADRUPOLE: os << "Quadrupole "; break;
    case SEXTUPOLE: os << "Sextupole "; break;
    case SBEND: os << "Sbend "; break;
    case MAP: os << "Map "; break;
    case PARTMARKER: os << "Partmarker "; break;
    case STATMARKER: os << "Statmarker "; break;
    case SCMAP: os << "Space Charge Kick "; break;
    case RECQUADRUPOLE: os << "Rare Earth Cobalt Quadrupole  "; break;	
    case IBS: os << "Intrabeam scattering  "; break;
    case RADIATION: os << "Radiation damping "; break;
    case RFCAVITY: os << "RFCavity "; break;
    default:    os << "Not known ";
    }
    if (length != 0.0)
      os << " L= " << length;
    if (slice != 1.0)
      os << " S= " << slice;
    if (alpha != 0.0)
      os << " alpha= " << alpha;
    if (beta != 0.0)
      os << " beta= " << beta;
    if (lambda != 0.0)
      os << " lambda= " << lambda;
    if (r1 != 0.0)
      os << " r1= " << r1;
    if (r2 != 0.0)
      os << " r2= " << r2;
    if (K0 != 0.0)
      os << " K0= " << K0;
    if (K0S != 0.0)
      os << " K0S= " << K0S;
    if (K1 != 0.0)
      os << " K1= " << K1;
    if (K1S != 0.0)
      os << " K1S= " << K1S;
    if (K2 != 0.0)
      os << " K2= " << K2;
    if (K2S != 0.0)
      os << " K2S= " << K2S;
    if (K3 != 0.0)
      os << " K3= " << K3;
    if (K3S != 0.0)
      os << " K0S= " << K3S; 
    if (e1 != 0.0)
      os << " e1= " << e1; 
    if (e2 != 0.0)
      os << " e2= " << e2;
    if (angle != 0.0)
      os << " angle= " << angle; 
    if (fint != 0.0)
      os << " fint= " << fint; 
    if (hgap != 0.0)
      os << " hgap= " << hgap; 
    os << " a= " << a;
    os << " b= " << b;
    os << endl;
    
  }

  string convStr() {
    string s;
    switch (elType) {
    case DRIFT: s="Drift"; break;
    case QUADRUPOLE: s="Quadrupole"; break;
    case SEXTUPOLE: s="Sextupole"; break;
    case SBEND: s="Sbend"; break;
    case MAP: s="Map"; break;
    case PARTMARKER: s="Partmarker"; break;
    case STATMARKER: s="Statmarker"; break;
    case SCMAP: s="Space charge kick"; break;
    case RECQUADRUPOLE: s="RECQuadrupole"; break;
    case IBS: s="Intrabeam scattering"; break;  
    case RADIATION: s="Radiation damping"; break;
    case RFCAVITY: s="RFCavity "; break;
    case DIPOLE: s="Dipole "; break;
    case SOLENOID: s="Solenoid "; break;
    case CONSTFOC: s="Constant focussing channel "; break;
    case ROTATION: s="Rotation "; break;
    default: s= "Not known";
    }
    return s;
  }

  double K0;
  double K0S;
  double K1;
  double K1S;
  double K2;
  double K2S;
  double K3;
  double K3S;
  double length;
  double slice;
  double alpha;
  double beta;
  double fieldindex;
  double e1; double e2; double angle, fint, hgap;
  unsigned int order;
  double a,b;   
  ElementT elType;
  UnitsT uType;
  double r1;
  double r2;
  double entrance;
  double v;
  double phi;
  double freq;
  double lambda;
  bool isMisaligned;
  double missX;
  double missY;
  int elemNo;
};

// Output operator.
inline ostream &operator<<(ostream &os, const ElemData &data)
{
  data.print(os);
  return os;
}

class SimCfgData 
{
public:
  void print(ostream &os) const
  {
    os << "-------------- GenPaTra ----------------------- " << endl << endl;
    if (uType==MAD9)
      os << "Using MAD9  units " << endl;
    else if (uType==ELCL)
      os << "Using Electron-Cloud units " << endl;
    else
      os << "Using MARYLIE units " << endl;
    if (specis==PROTONS)
      os << "Particles PROTONS " << endl;
    else
      os << "Particles ELECTRONS " << endl;
    if (indepVar==TIME) {
      os << "Using time as the independent variable, integrator ";
      if (integrator==VERLET)
	os << "Verlet 3 order" << endl;
      else if (integrator==RK4)
	os << "RK4 " << endl;
      else 
	os << "LEAPFROG " << endl;
      if (timeinvers)
	os << "Timeinversion set " << endl;
      if (dt>0.0)
	os << "Using fixed dt= " << dt << " and " << ndt << " integrations steps per dt" << endl;
    }
    else 
      os << "Using arclenght as the independent variable, integrator SPLIT" << endl;
    if (trackSymplectic)
      os << "Track symplectic using f3 and f4" << endl;
    else
      os << "Track symplectic using linear map only" << endl;
    os << "Simulation input file " << simInputFileName << endl;
    os << "Separating particle 0 and 1 with " << deltaEPS << endl;
    if (fsType == MultigridFEM)
      os << "Multigrid FEM using Gauss Seidel " << endl;
    if (doIBS)
      os << "Intrabeamscattering using " << dipoleIBSKicks << " kicks per dipole" << endl; 
    if(doRadiation) {
      os << "Radiation damping using a ";
      switch (radiationProcess) {
        case QUANTUM: os << " quantum machanical approach" << endl; break;
        case CLASSIC: os << " classical approach" << endl; break;
        default: os << "process not known" << endl;
      }
    }
    os << "-----------End GenPaTra ----------------------- " << endl << endl;
  }
  double Ek;
  double frequ;
  double I0;
  unsigned long periods;
  unsigned int order;

  unsigned long nPInit;    // initial number of electrons
  unsigned long nEInit;    // initial number of protons
  unsigned int  subTSteps;  // number of subtimesteps for the electrons

  unsigned int particlesToSave;
  InterPolT interPol;
  BCT bc;
  DistT distrType;
  
  double centerx;
  double centery;
  double centert;

  unsigned int nx, ny, nz;
  
  double Qtot;
   
  FsT fsType;
    
  InteractionT interaction;
    
  string title;
  
  string diagfile;
  
  string datadir;
  
  string distrInputFile;
 
  string simInputFileName;
  
  string elclInputFile;             // electron cloud cfg parameter
  
  vector<string> vFns;              // all used (base) filenames

  bool Usequad;

  UnitsT uType;

  bool trackSymplectic;
  bool trackOneTurnMap;

  // defines the separation of the first read in particle
  // from partInput.dat with the second particle in terms of DBL_EPSILON
  double deltaEPS;  

  IndepVarT indepVar;
  
  IntegratorT integrator;

  SpecisT specis;

  double dt;             // defines the dt for the timeintegration
			 // overwrites the L on the beamline element 

  unsigned int ndt;       // defines the number of integrations steps per dt

  bool timeinvers;       // after n periods run the beam backwards
			 // only implemented with timeintegration

  bool latticeTest;      // integrate the drive beam only


  /*
     EXPDE variables
  */
  
  int    maxMgLevel;
  double hInit;

  double relErr;
  int    iter_relax;

  double grid_stretch;

  /*
    IBS parameters
  */

  bool doIBS;
  unsigned int dipoleIBSKicks;

  /*
    Radiation damping
  */
  
  bool doRadiation;
  RadiationT radiationProcess;


  /*
    Electron cloud related stuff

  */

  double lambdaP; // proton line density
  double lambdaE; // electron line density
  double neutralfac; // neutralisation factor

  /*
    TTrack stuff

  */

  bool doTouschek;
  double Division;
  double fieldindex;
  double radius;
  double accpx, accpy, accps;  // Momentum Acceptance
  double accx;
  double dispersion;
  double lambda;

  double interrad;
  unsigned meshmultfactor;

  bool doRestart;
  bool doStatistics;
  bool doInbalance;
  bool doPhaseSpace;
  bool diagnostics;
  bool doLostPlace;
  bool doAllInteractions;
  bool doCrosssection;

  unsigned scattering;
  unsigned distributiontype;

  double thetamin;
  double sigmawidth;
  double rfvoltage;
};

inline ostream &operator<<(ostream &os, const SimCfgData &data)
{
  data.print(os);
  return os;
}

class DistrData {
public:
    
    DistT distrType;
    unsigned long N;
    
    double sigx;
    double sigpx;
    double corxpx;
    double emitx;
    double mx;
  
    double sigy;
    double sigpy;
    double corypy;
    double emity;
    double my;
    
    double sigt;
    double sigpt;
    double cortpt;
    double emitt;
    double mt;
};


/* 
   Helper Class to store 
   1. Element descriotion/data
   2. Map
   3. s value 
*/

class ElMap {

public:

#ifdef MAPS
  vector< FVpsT > taylor;
  vector< FMatrix<double,2*DIMENSION,2*DIMENSION> > Msymp;
  vector< DragtFinnMap<DIMENSION> > dfMapSimple;
#endif

  FMatrix<double,2*DIMENSION,2*DIMENSION> M;

  double actualS;
  ElemData elemData;
  /**
     Index of last element in lattice
     lattice: 0 ..... lastElIdx
  */
  int lastElIdx;      
};
#endif	
