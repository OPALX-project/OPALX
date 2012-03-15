//## begin module%1.2%.codegen_version preserve=yes
//   Read the documentation to learn more about C++ code generator
//   versioning.
//## end module%1.2%.codegen_version

//## begin module%3AE7DDF3017C.cm preserve=no
/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\BeamDynamics\ParticleTracking\SynchRadParticleProcess.h
 * last modified 09/06/01 01:50:53 PM
 */
//## end module%3AE7DDF3017C.cm

//## begin module%3AE7DDF3017C.cp preserve=no
/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 2001 by The Merlin Collaboration.
 * - ALL RIGHTS RESERVED - 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */
//## end module%3AE7DDF3017C.cp

//## Module: SynchRadParticleProcess%3AE7DDF3017C; Package specification
//## Subsystem: Merlin::BeamDynamics::ParticleTracking%3729F894017C
//## Source file: C:\C++\Merlin\BeamDynamics\ParticleTracking\SynchRadParticleProcess.h

#ifndef SynchRadParticleProcess_h
#define SynchRadParticleProcess_h 1

//## begin module%3AE7DDF3017C.additionalIncludes preserve=no
#include "merlin_config.h"
//## end module%3AE7DDF3017C.additionalIncludes

//## begin module%3AE7DDF3017C.includes preserve=yes
//## end module%3AE7DDF3017C.includes

// ParticleBunchProcess
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"
//## begin module%3AE7DDF3017C.additionalDeclarations preserve=yes

// Global functions for photon spectrum generation

// Generator from Helmut Burkhardt (CERN)
// (default generator, see also CERN-LEP-Note 632)
double HBSpectrumGen(double u);

// Generator from Andy Wolski (LBL). 
// Faster than HBSpectrumGen but not as accurate.
double AWSpectrumGen(double u);

//## end module%3AE7DDF3017C.additionalDeclarations


//## Class: SynchRadParticleProcess%37CB8D2103C7
//	Models the effects of synchrotron radiation in dipoles
//	and (optionally) quadrupoles. The default behaviour is
//	to include only the classical (average) effects. Calling
//	GeneratPhotons(true) models the quantum effects due to
//	random photon generation. The photon spectrum used can
//	be changed by a call to SetPhotonGenerator(double
//	(*)(double u)). The default spectrum (dipole radiation)
//	has been provided by H. Burkhardt (CERN-LEP-Note 632).
//
//	The number of equally spaced steps to take through a
//	component can be specified (default = 1). The effect of
//	the energy loss on the particles can be specified in two
//	ways: first (default), only the dp/p are adjusted, and
//	the reference energy (momentum) of the bunch is left
//	unchanged; alternatively, the reference energy is adjusted to
//	reflect the mean energy loss, and the dp/p then reflect
//	the additional change with respect to the mean.
//## Category: Merlin::BeamDynamics::ParticleTracking%36E7A9AD0176
//## Subsystem: Merlin::BeamDynamics::ParticleTracking%3729F894017C
//## Persistence: Transient
//## Cardinality/Multipplicity: n

//## Uses: <unnamed>%37CB906500DA;SectorBend { -> }
//## Uses: <unnamed>%3B949CDF0074;RandomNG { -> }
//## Uses: <unnamed>%3B94AA2001F9;RectMultipole { -> }

class SynchRadParticleProcess : public ParticleBunchProcess  //## Inherits: <unnamed>%3AE7D8A10118
{
  public:
    //## Class: PhotonGenerator%3B94842B0201
    //## Category: Merlin::BeamDynamics::ParticleTracking%36E7A9AD0176
    //## Subsystem: Merlin::BeamDynamics::ParticleTracking%3729F894017C
    //## Persistence: Transient
    //## Cardinality/Multipplicity: n

    typedef double (*PhotonGenerator)(double);

    //## Constructors (specified)
      //## Operation: SynchRadParticleProcess%37CB8E790187
      //	Constructor taking the process priority (>=0), the
      //	number of equal steps to take through a component
      //	(nstep). If the flag q==true, then quantum excitation is
      //	modelled.
      SynchRadParticleProcess (int prio, int nstep = 1, bool q = false);


    //## Other Operations (specified)
      //## Operation: SetCurrentComponent%37CB8F370270
      //	Sets the current accelerator component. If component is
      //	a SectorBend, then the process becomes active.
      virtual void SetCurrentComponent (AcceleratorComponent& component);

      //## Operation: DoProcess%37CB8F3702DF
      //	Preform the process for the specified step ds.
      virtual void DoProcess (double ds);

      //## Operation: GetMaxAllowedStepSize%37CB9DBB03BE
      //	Returns the current maximum step length for this process.
      virtual double GetMaxAllowedStepSize () const;

      //## Operation: IncludeQuadRadiation%3B9485C401EB
      //	Include radiation effects in Quadrupoles and Skew
      //	Quadrupoles.
      void IncludeQuadRadiation (bool quadsr);

      //## Operation: SetPhotonGenerator%3B9485D8013F
      //	Sets the function to be called to generate the photon
      //	spectrum.
      static void SetPhotonGenerator (PhotonGenerator pg)
      {
        //## begin SynchRadParticleProcess::SetPhotonGenerator%3B9485D8013F.body preserve=yes
		  pgen = pg;
        //## end SynchRadParticleProcess::SetPhotonGenerator%3B9485D8013F.body
      }

      //## Operation: SetNumComponentSteps%3B9486210271
      //	Sets the minimum number of equal steps to take through
      //	the component.
      void SetNumComponentSteps (int n);

      //## Operation: GeneratePhotons%3B94888E004D
      //	If gp==true, include quantum excitation effects via
      //	photon generation.
      void GeneratePhotons (bool gp);

      //## Operation: AdjustBunchReferenceEnergy%3B9761CC0218
      //	If flg==true, the reference energy (momentum) of the
      //	ParticleBunch is adjusted to the mean of the particle
      //	energies. If false, then only the dp/p are adjusted.
      void AdjustBunchReferenceEnergy (bool flg)
      {
        //## begin SynchRadParticleProcess::AdjustBunchReferenceEnergy%3B9761CC0218.body preserve=yes
		  adjustEref = flg;
        //## end SynchRadParticleProcess::AdjustBunchReferenceEnergy%3B9761CC0218.body
      }

    // Data Members for Associations

      //## Association: Merlin::BeamDynamics::ParticleTracking::class details::<unnamed>%3B948179006E
      //## Role: SynchRadParticleProcess::pgen%3B948179024F
      //	Pointer to the function which generates the random
      //	photon spectrum. The function should have the form
      //
      //	double (*)(double u)
      //
      //	where u is the critical photon energy. Returns the
      //	energy of the photon in GeV.
      //## begin SynchRadParticleProcess::pgen%3B948179024F.role preserve=no  public: static SynchRadParticleProcess::PhotonGenerator { -> 1VHgAN}
      static PhotonGenerator pgen;
      //## end SynchRadParticleProcess::pgen%3B948179024F.role

  public:
  protected:
  private:
    // Data Members for Class Attributes

      //## Attribute: ns%37CB8F090260
      //## begin SynchRadParticleProcess::ns%37CB8F090260.attr preserve=no  private: int {UA} 
      int ns;
      //## end SynchRadParticleProcess::ns%37CB8F090260.attr

      //## Attribute: quantum%3B948E66008A
      //## begin SynchRadParticleProcess::quantum%3B948E66008A.attr preserve=no  private: PhotonGenerator {UA} 
      PhotonGenerator quantum;
      //## end SynchRadParticleProcess::quantum%3B948E66008A.attr

      //## Attribute: incQ%3B948F590238
      //## begin SynchRadParticleProcess::incQ%3B948F590238.attr preserve=no  private: bool {UA} 
      bool incQ;
      //## end SynchRadParticleProcess::incQ%3B948F590238.attr

      //## Attribute: adjustEref%3B9761B400BF
      //## begin SynchRadParticleProcess::adjustEref%3B9761B400BF.attr preserve=no  private: bool {UA} 
      bool adjustEref;
      //## end SynchRadParticleProcess::adjustEref%3B9761B400BF.attr

  private: //## implementation
    // Data Members for Class Attributes

      //## Attribute: ds2%37CB90890050
      //	The rms momentum addition per unit length.
      //## begin SynchRadParticleProcess::ds2%37CB90890050.attr preserve=no  implementation: double {UA} 
      double ds2;
      //## end SynchRadParticleProcess::ds2%37CB90890050.attr

      //## Attribute: ds1%37CB909E038F
      //	The mean energy loss per unit length.
      //## begin SynchRadParticleProcess::ds1%37CB909E038F.attr preserve=no  implementation: double {UA} 
      double ds1;
      //## end SynchRadParticleProcess::ds1%37CB909E038F.attr

      //## Attribute: dL%37CB95180238
      //## begin SynchRadParticleProcess::dL%37CB95180238.attr preserve=no  implementation: double {UA} 
      double dL;
      //## end SynchRadParticleProcess::dL%37CB95180238.attr

      //## Attribute: nk1%37CB99B80263
      //## begin SynchRadParticleProcess::nk1%37CB99B80263.attr preserve=no  implementation: int {UA} 
      int nk1;
      //## end SynchRadParticleProcess::nk1%37CB99B80263.attr

      //## Attribute: intS%37CB99C10162
      //## begin SynchRadParticleProcess::intS%37CB99C10162.attr preserve=no  implementation: double {UA} 
      double intS;
      //## end SynchRadParticleProcess::intS%37CB99C10162.attr

    // Data Members for Associations

      //## Association: Merlin::BeamDynamics::ParticleTracking::<unnamed>%3B94A81D03C5
      //## Role: SynchRadParticleProcess::currentField%3B94A81D03C7
      //## begin SynchRadParticleProcess::currentField%3B94A81D03C7.role preserve=no  implementation: MultipoleField { -> RHAN}
      MultipoleField* currentField;
      //## end SynchRadParticleProcess::currentField%3B94A81D03C7.role

};

//## begin module%3AE7DDF3017C.epilog preserve=yes
//## end module%3AE7DDF3017C.epilog


#endif
