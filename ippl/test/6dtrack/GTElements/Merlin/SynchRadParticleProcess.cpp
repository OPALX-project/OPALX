//## begin module%1.2%.codegen_version preserve=yes
//   Read the documentation to learn more about C++ code generator
//   versioning.
//## end module%1.2%.codegen_version

//## begin module%3AE7DE070366.cm preserve=no
/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\BeamDynamics\ParticleTracking\SynchRadParticleProcess.cpp
 * last modified 09/06/01 01:50:54 PM
 */
//## end module%3AE7DE070366.cm

//## begin module%3AE7DE070366.cp preserve=no
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
//## end module%3AE7DE070366.cp

//## Module: SynchRadParticleProcess%3AE7DE070366; Package body
//## Subsystem: Merlin::BeamDynamics::ParticleTracking%3729F894017C
//## Source file: C:\C++\Merlin\BeamDynamics\ParticleTracking\SynchRadParticleProcess.cpp

//## begin module%3AE7DE070366.includes preserve=yes
#include <cmath>
#include <algorithm>
#include "NumericalUtils/utils.h"
//## end module%3AE7DE070366.includes

// SynchRadParticleProcess
#include "BeamDynamics/ParticleTracking/SynchRadParticleProcess.h"
// SectorBend
#include "AcceleratorModel/StdComponent/SectorBend.h"
// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"
// RandomNG
#include "Random/RandomNG.h"
//## begin module%3AE7DE070366.additionalDeclarations preserve=yes
namespace {

	struct ApplySR {

		typedef double (*spec_gen)(double);

		const MultipoleField& Bf;
		double dL;
		double P0;
		spec_gen photgen;

		double meanU;
		size_t n;

		ApplySR(const MultipoleField& field, double dl, double p0, spec_gen sg =0)
			: Bf(field),dL(dl),P0(p0),photgen(sg),meanU(0),n(0) {};

		double MeanEnergyLoss() const { return meanU/n; }

		void operator()(PSvector& v)
		{
			double B = abs(Bf.GetField2D(v.x(),v.y()));
			double g = P0*(1 + v.dp())/510.999e-6;
			double uc = 1.73651479874e-13*B*g*g;
			double u = 0;

			if(photgen) {
				int nphot = static_cast<int>(RandomNG::poisson(6.17938614443*B*dL));
				for(int n=0; n<nphot; n++) 
					u += photgen(uc);
			}
			else
				u = 1.90275746875 * B * dL * uc;

			meanU += u;
			v.dp() -= u / P0;
			n++;
		}

	};
}; // end of anonymous namespace
//## end module%3AE7DE070366.additionalDeclarations


// Class SynchRadParticleProcess 











//## begin SynchRadParticleProcess::pgen%3B948179024F.role preserve=no  public: static SynchRadParticleProcess::PhotonGenerator { -> 1VHgAN}
SynchRadParticleProcess::PhotonGenerator SynchRadParticleProcess::pgen = HBSpectrumGen;
//## end SynchRadParticleProcess::pgen%3B948179024F.role

//## Operation: SynchRadParticleProcess%37CB8E790187
SynchRadParticleProcess::SynchRadParticleProcess (int prio, int nstep, bool q)
  //## begin SynchRadParticleProcess::SynchRadParticleProcess%37CB8E790187.initialization preserve=yes
  : ParticleBunchProcess("SYNCHROTRON RADIATION",prio),ns(nstep),incQ(false),adjustEref(true)
  //## end SynchRadParticleProcess::SynchRadParticleProcess%37CB8E790187.initialization
{
  //## begin SynchRadParticleProcess::SynchRadParticleProcess%37CB8E790187.body preserve=yes
	GeneratePhotons(q);
  //## end SynchRadParticleProcess::SynchRadParticleProcess%37CB8E790187.body
}



//## Other Operations (implementation)
//## Operation: SetCurrentComponent%37CB8F370270
void SynchRadParticleProcess::SetCurrentComponent (AcceleratorComponent& component)
{
  //## begin SynchRadParticleProcess::SetCurrentComponent%37CB8F370270.body preserve=yes
	SectorBend* bend =0;
	RectMultipole* rmult =0;

	if(bend = dynamic_cast<SectorBend*>(&component)) 
		currentField =  &(bend->GetField());
	else if(incQ && (rmult = dynamic_cast<RectMultipole*>(&component)))
		currentField = &(rmult->GetField());
	else
		currentField = 0;
	
	dL = component.GetLength()/ns;
	nk1=0;
	intS=0;

	active = currentField && currentBunch;
  //## end SynchRadParticleProcess::SetCurrentComponent%37CB8F370270.body
}

//## Operation: DoProcess%37CB8F3702DF
void SynchRadParticleProcess::DoProcess (double ds)
{
  //## begin SynchRadParticleProcess::DoProcess%37CB8F3702DF.body preserve=yes
	if(fequal(intS+=ds,(nk1+1)*dL))
	{
		double E0 = currentBunch->GetReferenceMomentum();
		double meanU = for_each(
			currentBunch->begin(),
			currentBunch->end(),
			ApplySR(*currentField,dL,E0,quantum)).MeanEnergyLoss();

		// Finally we adjust the reference of the
		// bunch to reflect the mean energy loss
		// of all the particles
		if(adjustEref)
			currentBunch->AdjustRefMomentum(-meanU/E0);
		nk1++;
	}
	active = nk1!=ns;
  //## end SynchRadParticleProcess::DoProcess%37CB8F3702DF.body
}

//## Operation: GetMaxAllowedStepSize%37CB9DBB03BE
double SynchRadParticleProcess::GetMaxAllowedStepSize () const
{
  //## begin SynchRadParticleProcess::GetMaxAllowedStepSize%37CB9DBB03BE.body preserve=yes
	return (nk1+1)*dL-intS;
  //## end SynchRadParticleProcess::GetMaxAllowedStepSize%37CB9DBB03BE.body
}

//## Operation: IncludeQuadRadiation%3B9485C401EB
void SynchRadParticleProcess::IncludeQuadRadiation (bool quadsr)
{
  //## begin SynchRadParticleProcess::IncludeQuadRadiation%3B9485C401EB.body preserve=yes
	incQ = quadsr;
  //## end SynchRadParticleProcess::IncludeQuadRadiation%3B9485C401EB.body
}

//## Operation: SetNumComponentSteps%3B9486210271
void SynchRadParticleProcess::SetNumComponentSteps (int n)
{
  //## begin SynchRadParticleProcess::SetNumComponentSteps%3B9486210271.body preserve=yes
	ns = n;
  //## end SynchRadParticleProcess::SetNumComponentSteps%3B9486210271.body
}

//## Operation: GeneratePhotons%3B94888E004D
void SynchRadParticleProcess::GeneratePhotons (bool gp)
{
  //## begin SynchRadParticleProcess::GeneratePhotons%3B94888E004D.body preserve=yes
	quantum = gp ? pgen : 0;
  //## end SynchRadParticleProcess::GeneratePhotons%3B94888E004D.body
}

//## begin module%3AE7DE070366.epilog preserve=yes
//## end module%3AE7DE070366.epilog
