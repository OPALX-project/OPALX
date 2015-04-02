/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file ParticleProcesses\SynchRadParticleProcess.cpp
 * last modified 05/10/00  03:19:13
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

#include "NumericalUtils/utils.h"
#include "Random/RandomNG.h"

// IBSProcess
#include "IBSProcess/IBSProcess.h"

#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/PhysicalUnits.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

namespace
{
	struct IBSKick
	{
		double xpScale, ypScale, zpScale;
		IBSKick(double xpscale, double ypscale, double zpscale) : 
			xpScale(xpscale), ypScale(ypscale), zpScale(zpscale) {};

		void operator()(PSvector& v)
		{
			double du = RandomNG::normal(0,zpScale);

			v.dp() += du;
//			v.xp() /= (1. + du);
//			v.yp() /= (1. + du);
//			cout<<std::setw(14)<<du/v.dp()<<endl;
		}

	};
};

IBSProcess::IBSProcess (int prio, int nkick)
  : ParticleBunchProcess("INTRABEAM SCATTERING",prio),dL(0),nk(nkick) {}

void IBSProcess::InitialiseProcess (Bunch& bunch)
{
	currentBunch = dynamic_cast<ParticleBunch*>(&bunch);
}

void IBSProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	active = currentBunch ? true : false;
	dL = component.GetLength()/nk;
	nk1 = 0;
	intS = 0;
}

void IBSProcess::DoProcess (double ds)
{
	static const double cIbs = 5.0419e-31; //3.565e-31;
	static const double nPart = 1.5e10;
	static const double gamma = 3874;
	static const double cLog = 21;

	if(fequal(intS+=ds,(nk1+1)*dL) && (ds!=0.0))
	{
		PSmoments sigma2;
		currentBunch->GetMoments(sigma2);

		double vol = sigma2(0,0);
		for(int i=1; i<5; i++) vol *= sigma2(i,i);
		vol = sqrt(vol);
		double gamma3 = gamma * gamma * gamma;

		double xpScale = cIbs * nPart * cLog * ds / vol / gamma3;
		double ypScale = xpScale;
		double zpScale = xpScale;

	//	cout<<std::setw(14)<<sqrt(sigma2(5,5));
	//	cout<<std::setw(14)<<sqrt(zpScale)<<endl;
	//	xpScale *= DIntegralX( sqrt(sigma2(1,1)), sqrt(sigma2(3,3)) ) / sigma2(1,1);
	//	ypScale *= DIntegralY( sqrt(sigma2(1,1)), sqrt(sigma2(3,3)) ) / sigma2(3,3);
		zpScale *= DIntegral1( sqrt(sigma2(1,1)), sqrt(sigma2(3,3)) );

//		cout<<zpScale*SpeedOfLight/2/sigma2(5,5)/ds<<endl;

		for_each(currentBunch->begin(),currentBunch->end(),IBSKick(xpScale,ypScale,zpScale));
		nk1++;
	}
	active = nk1!=nk;

}

double IBSProcess::GetMaxAllowedStepSize () const
{
	return (nk1+1)*dL-intS;
}

double IBSProcess::DIntegral1(double sxp, double syp)
{
	double s = 2 * syp * syp / sxp / sxp;

	double c0 = -0.03613409093665682;
	double c1 = 0.00819261310489227;
	double c2 = 0.4491579033390567;
	double c3 = -3.579833752175229;
	double c4 = 10.88825269207325;
	double c5 = -13.53400438330434;
	double c6 = 5.880436109436087;

	double r = c0 + 
		c1*pow(s,1./2) + 
		c2*pow(s,1./4) + 
		c3*pow(s,1./6) + 
		c4*pow(s,1./8) +
		c5*pow(s,1./10) +
		c6*pow(s,1./12);

	return syp / r / 2.;
}

double IBSProcess::DIntegral2(double sxp, double syp)
{
	double s = 2 * syp * syp / sxp / sxp;

	double c0 = -1.965825326396189;
	double c1 = 0.01583009512725006;
	double c2 = -2.446379409023215;
	double c3 = 65.21655565099193;
	double c4 = -432.8300959567538;
	double c5 = 1127.770210422389;
	double c6 = -1255.995360621286;
	double c7 = 500.4016992120437;

	double r = c0 + 
		c1*pow(s,1./2) + 
		c2*pow(s,1./4) + 
		c3*pow(s,1./6) + 
		c4*pow(s,1./8) +
		c5*pow(s,1./10) +
		c6*pow(s,1./12) +
		c7*pow(s,1./14);

	return syp / r / 2.;
}

double IBSProcess::DIntegralX(double sxp, double syp)
{
	return DIntegral1(sxp,syp) - 3 * DIntegral2(syp,sxp);
}

double IBSProcess::DIntegralY(double sxp, double syp)
{
	return DIntegral1(sxp,syp) - 3 * DIntegral2(sxp,syp);
}

double IBSProcess::DIntegralZ(double sxp, double syp)
{
	return DIntegral1(sxp,syp);
}
