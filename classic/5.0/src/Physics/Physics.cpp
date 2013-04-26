// ------------------------------------------------------------------------
// $RCSfile: Physics.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Physics:
//   Contains the definitions for the mathematical and physical constants.
//
// ------------------------------------------------------------------------
// Class category: Physics
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Physics/Physics.h"

// Namespace Physics
// ------------------------------------------------------------------------

namespace Physics {

    // Universal mathematical constants
    // ----------------------------------------------------------------------
    const double pi         = 3.14159265358979323846;
    const double two_pi     = 2.0 * pi;
    const double u_two_pi   = 1.0 / (2.0 * pi);
    const double e          = 2.7182818284590452354;
    const double log10e     = 0.43429448190325182765;


    // Universal physical constants
    // ----------------------------------------------------------------------

    // Global constants:
    const double c          = 299792458.0;
    const double mu_0       = 1.256637061e-06;
    const double epsilon_0  = 8.854187817e-12;
    const double h_bar      = 6.5821220e-25;
    const double Avo        = 6.022e23;
    const double kB         = 8.6173324e-5;

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
    const double z_p        = 1;

    // constants for carbon:
    const double m_c        = 12 * 0.931494027e+00;

    // constants for uranium:
    const double m_u        = 238 * 0.931494027e+00;

    // constants for muon:
    const double m_mu       = 0.1057;

    // constants for minus hydrogen:
    const double m_hm       = 0.939277e+00;

    // constants for deuteron
    const double m_d        = 2*0.931494027e+00;

    // constants for xenon
    const double m_xe       = 124*0.931494027e+00;

    const double PMASS      = 1.6726231e-27;  // kg
    const double EMASS      = 9.1093897e-31; // kg
    const double PCHARGE    = 1.6021773349e-19; // C

    // BET constants
    const double Ia         = 17.045148e+03;
    const double e0mc       = 5.86679074042e+02;
    const double e0m        = 1.75881961664e+11;

};
