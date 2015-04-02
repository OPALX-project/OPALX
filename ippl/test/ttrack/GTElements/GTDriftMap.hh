#ifndef DRIFTMAP_HH
#define DRIFTMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTDriftMap.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTDriftMap.hh,v 1.1.1.1 2003/01/23 09:13:53 adelmann Exp $
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
// $Log: GTDriftMap.hh,v $
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

class DriftMap {
private:
  ElemData cfgData_m;
  PartData pdata_m;
  
public:
  DriftMap(ElemData cfgData, PartData pdata)
  { 
    cfgData_m = cfgData;
    pdata_m = pdata; 
    FTpsT::setGlobalTruncOrder(cfgData.order);
  };

  ~DriftMap() { };
     
  DragtFinnMap<DIMENSION> getDragtFinnMap(DragtFinnMap<DIMENSION> map);
  FVpsT getTaylorSeries(FVpsT map);

private:
  FTpsT getH();
};


FTpsT DriftMap::getH() {
  //  FVps<double,6> itsMap; //(cfgData_m.order,cfgData_m.order);
  double beta =  pdata_m.getBeta();
   
    
  // Step 1: Define variables.
  FTpsT px = FTpsT::makeVariable(SIXVect::PX);
  FTpsT py = FTpsT::makeVariable(SIXVect::PY);
  FTpsT H;
  //Step 2: Kinematic terms.
  if (cfgData_m.uType==MAD9) {
    FTpsT pt = FTpsT::makeVariable(SIXVect::PT) + 1.0;
    FTpsT pz = sqrt(pt*pt - px*px - py*py);
    double kin = pdata_m.getM() / pdata_m.getP();
    FTpsT E = sqrt(pt*pt + kin*kin) / beta;
    H = E - pz;
    // cout << "brr: kin= " << kin << " refT= " << cfgData_m.length/beta << " bet= " << beta;
    // cout << " m= " << pdata_m.getM() << " p= " << pdata_m.getP() << endl;
    
  } else if (cfgData_m.uType==MARYLIE) {
    FTpsT pt = FTpsT::makeVariable(SIXVect::PT);
    FTpsT Pt = 2.0*pt/beta;
    FTpsT Ptt = pt/beta;
    FTpsT pz = sqrt(1.0 + Pt + pt*pt - px*px - py*py);
    H = - pz;
  }
  return H;
}


DragtFinnMap<DIMENSION> DriftMap::getDragtFinnMap(DragtFinnMap<DIMENSION> map)
{
  FTpsT H;
  H = getH();
  double parlen = cfgData_m.length;
  DragtFinnMap<DIMENSION> theMap = DragtFinnMap<DIMENSION>::factorSimple(H * parlen);
  return  theMap;
  //return map.catenate(theMap);
}

FVpsT DriftMap::getTaylorSeries(FVpsT map)
{
  double parlen = cfgData_m.length;
  return ExpMap(- getH() * parlen).substitute(map);
    
}

#endif
