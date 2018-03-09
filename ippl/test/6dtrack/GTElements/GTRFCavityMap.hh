#ifndef RFCAVITYMAP_HH
#define RFCAVITYMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTRFCavityMap.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTRFCavityMap.hh,v 1.1.1.1 2003/01/23 09:13:53 adelmann Exp $
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
// $Log: GTRFCavityMap.hh,v $
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include "GTConst.hh"
#include "GTElemData.hh"

#include "Physics/Physics.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/DragtFinnMap.h"

using namespace Physics;

class RFCavityMap {
private:
  ElemData cfgData_m;
  PartData pdata_m;
  
public:
  RFCavityMap(ElemData cfgData, PartData pdata)
  { 
    cfgData_m = cfgData;
    pdata_m = pdata; 
    FTpsT::setGlobalTruncOrder(cfgData.order);
  };

  ~RFCavityMap() { };
    
  DragtFinnMap<DIMENSION> getDragtFinnMap(DragtFinnMap<DIMENSION> map);
  FVpsT getTaylorSeries(FVpsT map);

private:
  FTpsT getH();
};
 
FTpsT RFCavityMap::getH() {
  
  double freq = cfgData_m.freq*1.0e6;
  double P0 = pdata_m.getP();
  double peak = cfgData_m.v*1.0e6 / P0;
  double kin = pdata_m.getM() / P0;
  
  // Step 1: Define variables.
  FTpsT H;
  FTpsT t = FTpsT::makeVariable(SIXVect::TT);
  
  //Step 2: Kinematic terms.
  if (cfgData_m.uType==MAD9) {
    FTpsT pt = FTpsT::makeVariable(SIXVect::PT) + 1.0;
    FTpsT speed = (Physics::c * pt) / sqrt(pt*pt + kin*kin);
    FTpsT phase = cfgData_m.phi + (freq * t) / speed;
    H =  peak * sin(phase) / pt;
  } else if (cfgData_m.uType==MARYLIE) {
    FTpsT pt = FTpsT::makeVariable(SIXVect::PT);
    FTpsT speed = (Physics::c * pt) / sqrt(pt*pt + kin*kin);
    FTpsT phase = cfgData_m.phi + (freq * t) / speed;
    H =  peak * sin(phase) / pt;
  }
  return H;
}

DragtFinnMap<DIMENSION> RFCavityMap::getDragtFinnMap(DragtFinnMap<DIMENSION> map)
{
  FTpsT H;
  H = getH();
  double parlen = cfgData_m.length;
  DragtFinnMap<DIMENSION> theMap = DragtFinnMap<DIMENSION>::factorSimple(H * parlen);
  return  theMap;
}

FVpsT RFCavityMap::getTaylorSeries(FVpsT map)
{
  double parlen = cfgData_m.length;
  return ExpMap(- getH() * parlen).substitute(map);
}
#endif
