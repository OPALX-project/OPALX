#ifndef CLASSIC_Physics_HH
#define CLASSIC_Physics_HH

// ------------------------------------------------------------------------
// $RCSfile: Physics.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Physics
//
// ------------------------------------------------------------------------
// Class category: Physics
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class Physics
// ------------------------------------------------------------------------
/// A namespace defining various mathematical and physical constants.

namespace Physics {

  /// The value of $\pi$
  extern const double pi;

  /// The value of $2 \pi$
  extern const double two_pi;

  /// The value of $1 / (2 \pi)$
  extern const double u_two_pi;

  /// The value of $e$
  extern const double e;

  /// The logarithm of $e$ to the base 10
  extern const double log10e;

  /// The velocity of light in m/s
  extern const double c;

  /// The permeability of vacuum in Vs/Am
  extern const double mu_0;

  /// The permittivity of vacuum in As/Vm
  extern const double epsilon_0;

  /// The reduced Planck constant in GeVs
  extern const double h_bar;
 
  /// The elementary charge in As
  extern const double q_e;

  /// The fine structure constant, no dimension
  extern const double alpha;

  /// The electron rest mass in GeV
  extern const double m_e;

  /// The classical electron radius in m
  extern const double r_e;

  /// The reduced Compton wave length for electrons in m
  extern const double lamda_e;

  /// The magnetic momentum anomaly for electrons, no dimension
  extern const double a_e;

  /// The proton rest mass in GeV
  extern const double m_p;

  /// The carbon rest mass in GeV
  extern const double m_c;

  /// The uranium rest mass in GeV
  extern const double m_u;

  /// The carbon rest mass in GeV
  extern const double m_hm;

  /// The classical proton radius in m
  extern const double r_p;

  /// The reduced Compton wave length for protons in m
  extern const double lamda_p;

  /// The magnetic momentum anomaly for protons, no dimension
  extern const double a_p;

  extern const double EMASS;

  extern const double PMASS;
 
  extern const double PCHARGE;

};

#endif // CLASSIC_Physics_HH
