#ifndef RADIATIONMAP_HH
#define RADIATIONMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTRadiation.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTRadiation.hh,v 1.1.1.1 2003/01/23 09:13:53 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
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
// $Date: 2003/01/23 09:13:53 $
// $Author: adelmann $
// $Log: GTRadiation.hh,v $
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

/*
 * Adapted from Merlin C++ Class Library for Charged Particle Accelerator Simulations 
 * Class library version 2.0 (1999)
 * Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED. 
*/


#include "GTChargedParticles.hh"
#include "GTConst.hh"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"

using namespace Physics;

template <class K, unsigned int Dimension>
class GTRadiationProcess {
private:
  
  ChargedParticles<playout_t> *beam_m;
  gsl_rng * r_m;

  unsigned n;
  K mean_m;
  
  double AWSpectrumGen(double uc)
    {
      double u1 = IpplRandom();
      double u2 = 0.57 * uc * pow(u1,-1./3) * pow((1-u1),pi);
      return u2;
    }

public:
  GTRadiationProcess(ChargedParticles<playout_t> *bunch) :
    beam_m(bunch)
  { 
    
    /* create a generator chosen by the 
       environment variable GSL_RNG_TYPE */

    gsl_rng_env_setup();
    r_m = gsl_rng_alloc (gsl_rng_default);
      
    };
  
  ~GTRadiationProcess() { };
  
  void applyRadiation (ElemData cfgData) 
    {
      /*
	Some words on units:

	m_e = 0.510998902e-03 /GeV*c^{-2}
	Bv                    /T
	PO  =                 /eV*c^{-1};
	Q   = -1.0            
      */
      Inform msg("applyRadiation ");
      
      BMultipoleField field;
      double multipoleScalFactor = beam_m->pdata_m.getP() / Physics::c;
      double h = cfgData.angle/cfgData.length;
      field.setNormalComponent(1, multipoleScalFactor * h);          // ok 
      field.setSkewComponent  (1, multipoleScalFactor * cfgData.K0S);
      field.setNormalComponent(2, multipoleScalFactor * cfgData.K1);
      field.setSkewComponent  (2, multipoleScalFactor * cfgData.K1S);
      field.setNormalComponent(3, multipoleScalFactor * cfgData.K2  / 2.0);
      field.setSkewComponent  (3, multipoleScalFactor * cfgData.K2S / 2.0);
      field.setNormalComponent(4, multipoleScalFactor * cfgData.K3  / 6.0);
      field.setSkewComponent  (4, multipoleScalFactor * cfgData.K3S / 6.0);
      double P0 = beam_m->pdata_m.getP();
      double scale = Physics::c*beam_m->pdata_m.getQ()/P0;
      double dL = cfgData.length;
      double meanU = 0.0;
      unsigned n = 0;
      
      double globalmeanU = 0.0;
      msg << "P0= " << P0 << " scale= " << scale << " Q= " <<   beam_m->pdata_m.getQ() << endl;    
      for(unsigned i=0; i<beam_m->getLocalNum(); i++) {
	// Bv in Tesla
	BVector Bv = field.Bfield(Point3D(beam_m->R[i](0),beam_m->R[i](1),beam_m->R[i](2)));
	double B = sqrt(pow(Bv.getBx(),2)+pow(Bv.getBy(),2)); 
	double gamma = P0*(1.0 + beam_m->P[i](2))/Physics::m_e;
	double uc = 1.73651479874e-13*B*gamma*gamma;
	double u = 0;
	
	if(beam_m->simCfg_m.radiationProcess == QUANTUM) {
	  int nphot = static_cast<int>(gsl_ran_poisson(r_m, (6.17938614443*B*dL)));
	  for(int n=0; n<nphot; n++) 
	    u += AWSpectrumGen(uc);
	}
	else
	  u = 1.90275746875 * B * dL * uc;
	beam_m->P[i](2) -= u / P0;
	meanU += u;
	n++;
      }

      // mean energy loss
      meanU /=n;
      reduce(meanU, globalmeanU, OpAddAssign());
      
      // ajust reference energy
      
    }
};
#endif


