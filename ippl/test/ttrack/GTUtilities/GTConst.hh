#ifndef CONST_HH
#define CONST_HH
// ------------------------------------------------------------------------
// $RCSfile: GTConst.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTUtilities/GTConst.hh,v 1.2 2003/04/17 14:22:15 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
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
// $Date: 2003/04/17 14:22:15 $
// $Author: adelmann $
// $Log: GTConst.hh,v $
// Revision 1.2  2003/04/17 14:22:15  adelmann
// *** empty log message ***
//
// Revision 1.1.1.1  2003/01/23 09:13:58  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#define DIMENSION 3

namespace SPECIS {
  enum { PROTON, ELECTRON };
}

namespace SIXVect {
  enum { X, PX, Y, PY, TT, PT };
}

enum FsT {P2PSolver,FFTSolverOOO,TreeSolver,FFTSolverOOP,MGAMRSolver,MultigridFEM,None};
enum InterPolT {NGP,CIC};
enum IntegratorT {RK4,LEAPFROG,SPLIT,VERLET};
enum GreenT {Poor,Sin};
enum BCT {OOO,OOP,PPP,DDD,DDO,DDP,DIRICHLET,NEUMANN,INHOMDIRICHLET,MIXED};   

enum TestCaseT {SimpleDrift,Buncher};

enum DistT {RECTUNIFORM,ELLIPSOIDALUNIFORM,TWOELLIPSOIDALUNIFORM,READFROMFILE,
	    BHFORMAT,BHFORMATCOMPRESSED,BINOMINAL};

enum InteractionT {GRAVITY,COULOMB,ELECTRONCLOUD};


enum ElementT {PARTMARKER, STATMARKER, DRIFT, QUADRUPOLE, RECQUADRUPOLE, SEXTUPOLE, SBEND, VSBEND, 
	       RBEND, VRBEND, MAP, SCMAP, IBS, RADIATION, RFCAVITY, DIPOLE, SOLENOID, CONSTFOC, ROTATION, 
               CCFMAGNET, UNKNOWN};

enum UnitsT {MARYLIE, MAD9, ELCL};

enum IndepVarT {TIME,ARCLEN};

enum SpecisT {PROTONS,ELECTRONS,HMINUS};

enum PartStatT {ALIVE,DONE,ATTHEWALL,LOST};

enum RadiationT {QUANTUM,CLASSIC};


enum MatchState {    // State                  print for
  INTERNAL,          // Internal call          level > 2
  PROGRESS,          // New minimum found      level > 1
  START,             // First call             level > 0
  RESTART,           // Algorithm restarted    level > 0
  CHECK,             // Check required         level > 0
  CONVERGED,         // Match has converged    level > 0
  FAILED,            // Match failed           level > 0
  CALL_LIMIT,        // Call limit exceeded    level > 0
  ACCURACY_LIMIT,    // Tolerance too small    level > 0
  TERMINATED         // Terminated by user     level > 0
};

/*
  GenTrackE-ECL Constants
  
*/

#define   MAXSEC 10

#endif // CONST_HH
















