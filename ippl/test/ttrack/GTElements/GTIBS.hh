#ifndef IBSMAP_HH
#define IBSMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTIBS.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTIBS.hh,v 1.3 2004/04/28 08:07:37 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
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
// $Date: 2004/04/28 08:07:37 $
// $Author: adelmann $
// $Log: GTIBS.hh,v $
// Revision 1.3  2004/04/28 08:07:37  adelmann
// Add MAPS #define
//
// Revision 1.2  2004/04/05 13:08:26  adelmann
// Many changes to use new IPPL instead of POOMA
//
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

/*
 * Adapted from Merlin C++ Class Library for Charged Particle Accelerator Simulations 
 * Class library version 2.0 (1999)
 * Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED. 
*/

#include "Ippl.h"
//#include "GTChargedParticles.hh"

//typedef ParticleInteractLayout<double, 3> playout_t;


#include "GTConst.hh"
#include "Physics/Physics.h"

using namespace Physics;

class GTIBSProcess {
private:

  //  ChargedParticles< playout_t > *beam_m;
  
  double DIntegral1(double sxp, double syp)
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
  
  double DIntegral2(double sxp, double syp)
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

  double DIntegralX(double sxp, double syp)
  {
    return DIntegral1(sxp,syp) - 3 * DIntegral2(syp,sxp);
  }
  
  double DIntegralY(double sxp, double syp)
  {
    return DIntegral1(sxp,syp) - 3 * DIntegral2(sxp,syp);
  }
  
  double DIntegralZ(double sxp, double syp)
  {
    return DIntegral1(sxp,syp);
  }
  
public:
  GTIBSProcess()  //ChargedParticles<playout_t>  *bunch) //:
    //    beam_m(bunch)
  { 
  };
  
  ~GTIBSProcess() { };

  void applyIBSKick (double ds) 
  {
    static const double cIbs = 5.0419e-31; //3.565e-31;
    static const double nPart = 1.5e10;
    static const double gamma = 3874;
    static const double cLog = 21;
        
    //    beam_m->calcMoments();
    	
    double vol; // = beam_m->getMoment(0,0);
    //    for(int i=1; i<5; i++) vol *= beam_m->getMoment(i,i);
    vol = sqrt(vol);
    double gamma3 = gamma * gamma * gamma;
    
    double xpScale = cIbs * nPart * cLog * ds / vol / gamma3;
    double ypScale = xpScale;
    double zpScale = xpScale;
    
    //xpScale *= DIntegralX( sqrt(sigma2(1,1)), sqrt(sigma2(3,3)) ) / sigma2(1,1);
    //ypScale *= DIntegralY( sqrt(sigma2(1,1)), sqrt(sigma2(3,3)) ) / sigma2(3,3);

    //    zpScale *= DIntegral1(sqrt(beam_m->getMoment(1,1)), 
    //			  sqrt(beam_m->getMoment(3,3)));
    
    /*
      Note on the sign  beam_m->P[i](2) -=IpplRandom()*zpScale;
      
      For Merlin, the 5 variable (ct) means: ct > 0, the particle is ahead
      of the reference particle. In MAD9 convention, this is exact the oposite.
    */

    IpplRandom();
    //for(unsigned i=0; i<beam_m->getLocalNum(); i++) {
    //  beam_m->P[i](2) -=IpplRandom()*zpScale;
    //}
  }
};
#endif
