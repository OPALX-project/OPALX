/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file ParticleProcesses\SynchRadParticleProcess.h
 * last modified 05/10/00  03:19:12
 */

/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED. 
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

#ifndef IBSProcess_h
#define IBSProcess_h 1

#include "merlin_config.h"

// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"
// ParticleBunchProcess
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"

class IBSProcess : public ParticleBunchProcess
{
  public:

      IBSProcess (int prio, int nkick = 1);

      //	Initialise this process with bunch. bunch must be a
      //	ParticleBunch object, otherwise the process becomes
      //	inactive.
      virtual void InitialiseProcess (Bunch& bunch);

      //	Sets the current accelerator component. If component is
      //	a SectorBend, then the process becomes active.
      virtual void SetCurrentComponent (AcceleratorComponent& component);

      //	Preform the process for the specified step ds.
      virtual void DoProcess (double ds);

      //	Returns the current maximum step length for this process.
      virtual double GetMaxAllowedStepSize () const;

	  double DIntegralX(double sxp, double syp);

	  double DIntegralY(double sxp, double syp);

	  double DIntegralZ(double sxp, double syp);

	  double DIntegral1(double sxp, double syp);

	  double DIntegral2(double sxp, double syp);

  protected:
  private:
	  int nk, nk1;
	  double dL, intS;
      ParticleBunch* currentBunch;
};

#endif
