// ------------------------------------------------------------------------
// $RCSfile: pcyclint.h,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/pcyclint/prog/pcyclint.h,v 1.1.1.1 2004/09/22 12:10:44 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// $Date: 2004/09/22 12:10:44 $
// $Author: adelmann $
// $Log: pcyclint.h,v $
// Revision 1.1.1.1  2004/09/22 12:10:44  adelmann
// Imported sources pcyclubt
//
// Revision 1.1.1.1  2003/10/03 11:55:44  wittberger
// Marcus Wittbergers Summerwork 2003
//
// ------------------------------------------------------------------------

#ifndef PCYCLINT_H
#define PCYCLINT_H

#define rk23 1
#define rk45 2
#define rk78 3

#define NEQ 6
#define COEFFDIM 19
#define TSTART 0.0
#define TEND 1.5e-5
//#define dT 1.0e-11

//#define TOL 1e-6
 
#define e 1.6021773349e-19
#define mproton 1.67262311e-27
#define E0 938.0


const double pi         = 3.14159265358979323846;
const double two_pi     = 2.0 * pi;
const double u_two_pi   = 1.0 / (2.0 * pi);
//const double e          = 2.7182818284590452354;
const double log10e     = 0.43429448190325182765;

// Resonator geometry
const double geometry   = 1e6 / (4.97664e-2);
const double Racos      = cos( (90.0 - 34.9239) * 2.0 * pi / 360.0  );
const double Rasin      = sin( (90.0 - 34.9239) * 2.0 * pi / 360.0  );
const double Rbcos      = cos( (90.0 - 55.0546) * 2.0 * pi / 360.0  );
const double Rbsin      = sin( (90.0 - 55.0546) * 2.0 * pi / 360.0  );
const double Rccos      = cos( (270.0 - 215.119) * 2.0 * pi / 360.0 );
const double Rcsin      = sin( (270.0 - 215.119) * 2.0 * pi / 360.0 );
const double Rdcos      = cos( (270.0 - 234.9047) * 2.0 * pi / 360.0);
const double Rdsin      = sin( (270.0 - 234.9047) * 2.0 * pi / 360.0);

// Collimator settings
const double rrelKIP1    = 0.0;
const double rrelKIP2    = 0.0;
const double zrelKIG3    = 0.0;
const double zrelKIV     = 0.0;
const double rrelKIP4    = 0.0;
const double apertureKIP4= 0.011;
    
 

// Universal physical constants
// ----------------------------------------------------------------------

// Global constants:
const double c          = 299792458.0;
const double mu_0       = 1.256637061e-06;
const double epsilon_0  = 8.854187817e-12;
const double h_bar      = 6.5821220e-25;
  
  // electromagnetic constants:
const double q_e        = 1.60217733e-19;
const double alpha      = 7.29735308e-03;
  
// constants for electrons:
const double m_e        = 0.51099906e-03;
const double r_e        = 2.81794092e-15;
const double lamda_e    = 3.86159323e-13;
const double a_e        = 1.159652193e-03;
  
// constants for protons:
const double m_p        = 0.93827231e+00;
const double r_p        = 1.53469857e-18;
const double lamda_p    = 2.10308937e-16;
const double a_p        = 1.792847386e+00;

#endif
