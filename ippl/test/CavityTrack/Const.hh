// ------------------------------------------------------------------------
// $RCSfile: Const.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/pcyclint/prog/Const.hh,v 1.4 2004/10/01 20:28:19 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// $Date: 2004/10/01 20:28:19 $
// $Author: adelmann $
// $Log: Const.hh,v $
// Revision 1.4  2004/10/01 20:28:19  adelmann
// removed GAUSS2
//
// Revision 1.3  2004/09/28 04:41:08  adelmann
// Add  structures for cloning
//
// Revision 1.2  2004/09/24 19:13:35  adelmann
// Add stuff for Gaussian Distribution
//
// Revision 1.1.1.1  2004/09/22 12:10:44  adelmann
// Imported sources pcyclubt
//
// Revision 1.2  2003/10/23 04:39:19  adelmann
// Add DEG2RAD and RAD2DEG macros
//
// Revision 1.1.1.1  2003/10/03 11:55:44  wittberger
// Marcus Wittbergers Summerwork 2003
//
// ------------------------------------------------------------------------

#ifndef CONST_HH
#define CONST_HH
const double EPS0  = 8.854187817e-12; // F m^(-1)
const double MUE0  = 4 * 3.14159265358979323846 * 10e-7; // N A^(-2) , magnetic field constant
const double Q0    = 1.602176462e-19; // C
const double MPMEV = 938.271998;        // MeV c^(-2)
const double MP    = 0.938271998;  // GeV c^(-1)
const double CLIGHT = 299792458.0;  // ms^-1
const double CHARGE = 1.6021773349e-19; // C
const double MASS = 1.6726231e-27;  // kg  
const double KBOLTZ = 1.380650324e-23; // JK^-1

enum FsT {FFTSolver,P2PSolver,TreeSolver,FFTPeriodicSolver,None};
enum InterPolT {NGP,CIC};
enum IntegratorT {LEAPFROG,RK45,RK78,VELVERLET,RK4,RK4A};
enum GreenT {Poor,Sin};
enum BCT {OOO,OOP,PPP,DDD,DDO,DDP};   // OOO == all dim, open BC
                         // OOP == open transverse periodic BC in longitudinal BC
                         // PPP == all dim. periodic BC

enum DistT {RECTUniform,ELLIPSOIDALUniform,TwoELLIPSOIDALUniform,READFromFile,
	    BHFormat,BHFormatCompressed,BeamLets,BINOMINAL,ELLIPSOIDALUNIFORM,READFROMFILE,GAUSS,NOTKNOWN};

enum InteractionT {GRAVITY,COULOMB};


#define DEG2RAD 2.0*pi/360.0
#define RAD2DEG 360.0*(2.0*pi)

#define MAXCLONE 5      // maximal number of clones

const int PWIC = 8;

#endif // CONST_HH



