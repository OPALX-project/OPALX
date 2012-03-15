#ifndef MULTIPOLEMAP_HH
#define MULTIPOLEMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTMultipoleMap.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTMultipoleMap.hh,v 1.2 2004/04/05 13:08:27 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
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
// $Date: 2004/04/05 13:08:27 $
// $Author: adelmann $
// $Log: GTMultipoleMap.hh,v $
// Revision 1.2  2004/04/05 13:08:27  adelmann
// Many changes to use new IPPL instead of POOMA
//
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"

#include "GTElemData.hh"
#include "Algorithms/PartData.h"

using namespace Physics;

typedef FTps<double,6> Series;

class MultiPoleMap {
private:
  ElemData cfgData_m;
  PartData pdata_m;
  BMultipoleField field_m;
  double scale_m;
  double h_m;
public:
  MultiPoleMap(ElemData cfgData, PartData pdata)
  { 
    cfgData_m = cfgData;
    pdata_m = pdata; 
    FTps<double,6>::setGlobalTruncOrder(cfgData.order);
  
    double multipoleScalFactor = pdata_m.getP() / Physics::c;
    h_m = cfgData_m.angle/cfgData_m.length;
    field_m.setNormalComponent(1, multipoleScalFactor * h_m);          // ok 
    field_m.setSkewComponent  (1, multipoleScalFactor * cfgData_m.K0S);
    field_m.setNormalComponent(2, multipoleScalFactor * cfgData_m.K1);
    field_m.setSkewComponent  (2, multipoleScalFactor * cfgData_m.K1S);
    field_m.setNormalComponent(3, multipoleScalFactor * cfgData_m.K2  / 2.0);
    field_m.setSkewComponent  (3, multipoleScalFactor * cfgData_m.K2S / 2.0);
    field_m.setNormalComponent(4, multipoleScalFactor * cfgData_m.K3  / 6.0);
    field_m.setSkewComponent  (4, multipoleScalFactor * cfgData_m.K3S / 6.0);
    scale_m = Physics::c*pdata_m.getQ()/pdata_m.getP();
  };
  
  ~MultiPoleMap() { };
  
  DragtFinnMap<3> getDragtFinnMap(DragtFinnMap<3> map);
 
  FVps<double,6> getTaylorSeries(FVps<double,6> map);
  
  FVps<double,6> getEntranceFringeMap(double angle, double curve,
				     const BMultipoleField &field,
				     double scale, FVps<double,6> map);

  FVps<double,6> getExitFringeMap(double angle, double curve,
				 const BMultipoleField &field,
				 double scale, FVps<double,6> map);
  
  FVps<double,6> applyThinMultipole(const BMultipoleField &field, double scale, FVps<double,6> map);

private:
  
  Series getH();
  
  Series buildMultipoleVectorPotential(const BMultipoleField &);
  

};


DragtFinnMap<3> MultiPoleMap::getDragtFinnMap(DragtFinnMap<3> map)
{
  Series H;
  H = getH();
  double parlen = cfgData_m.length;
  DragtFinnMap<3> theMap = DragtFinnMap<3>::factorBerzForestIrwin(H * parlen);
  return  theMap;
  //return map.catenate(theMap);
}

FVps<double,6> MultiPoleMap::getTaylorSeries(FVps<double,6> map)
{

  Series  H;
  H = getH();
  double parlen = cfgData_m.length;
  return ExpMap(- getH() * parlen).substitute(map);
    
}


Series MultiPoleMap::getH()
{
  double beta =  pdata_m.getBeta();

 // std::cout << "Particle Data information" << std::endl;
 // std::cout << "Q= " << pdata_m.getQ() << " m= " << pdata_m.getM() << " beta= ";
 // std::cout << pdata_m.getBeta() << " gamma= " << pdata_m.getGamma() << std::endl;
 // std::cout << "E= " << pdata_m.getE() << " P0= " << pdata_m.getP() << std::endl << std::endl;
  	 
  // Step 1: Define variables.
  Series px = Series::makeVariable(SIXVect::PX);
  Series py = Series::makeVariable(SIXVect::PY);
  Series H;
  
  if (cfgData_m.uType==MAD9) {
    
    Series pt = Series::makeVariable(SIXVect::PT) + 1.0;
    
    // Step 2: Kinematic terms.
    Series pz = sqrt(pt*pt - px*px - py*py);
    double kin = pdata_m.getM() / pdata_m.getP();
    Series E = sqrt(pt*pt + kin*kin) / beta;
    
    // Step 3: Vector potential in straight reference.
    Series As = buildMultipoleVectorPotential(field_m) * scale_m;
    
    // Step 4: Finish Hamiltonian, build exponential map,
    //         and substitute previous map.
    H = As + E - pz;
    
  } else if (cfgData_m.uType==MARYLIE) {
    
    Series pt = Series::makeVariable(SIXVect::PT);
    Series Pt = 2.0*pt/beta;
    Series Ptt = pt/beta;
    Series pz = sqrt(1.0 + Pt + pt*pt - px*px - py*py);
    // Step 3: Vector potential in straight reference.
    Series As = buildMultipoleVectorPotential(field_m) * scale_m;
    // Step 4: Finish Hamiltonian, build exponential map,
    //         and substitute previous map.
    H =  As - pz + Ptt;
  }
  return H;
}


FVps<double,6> MultiPoleMap::getEntranceFringeMap(double angle, double curve,
					    const BMultipoleField &field,
					    double scale, FVps<double,6> map)
{
  // *** MISSING *** Higher order terms for entrance fringe.
  double hx = scale * field.normal(1);
  double ex = hx * tan(angle);
  double ey = hx * tan(angle + map[SIXVect::PX][0]);
  map[SIXVect::PX] += ex * map[SIXVect::X];
  map[SIXVect::PY] -= ey * map[SIXVect::Y];
  return map;
}

FVps<double,6> MultiPoleMap::getExitFringeMap(double angle, double curve,
					 const BMultipoleField &field,
					 double scale, FVps<double,6> map)
{
  // *** MISSING *** Higher order terms for exit fringe.
  double hx = scale * field.normal(1);
  double ex = hx * tan(angle);
  double ey = hx * tan(angle - map[SIXVect::PX][0]);
  map[SIXVect::PX] += ex * map[SIXVect::X];
  map[SIXVect::PY] -= ey * map[SIXVect::Y];
  return map;
}

FVps<double,6> MultiPoleMap::applyThinMultipole(const BMultipoleField &field, double scale, FVps<double,6> map) 
{
  cout << " not implemented " << endl;
  return map;
}

Series MultiPoleMap::buildMultipoleVectorPotential(const BMultipoleField &field)
{
  int order = field.order();
  
  if (order > 0) {
    static const Series x = Series::makeVariable(SIXVect::X);
    static const Series y = Series::makeVariable(SIXVect::Y);
    Series kx = + field.normal(order) / double(order);
    Series ky = - field.skew(order)   / double(order);
    
    while (order > 1) {
      Series kxt = x * kx - y * ky;
      Series kyt = x * ky + y * kx;
      order--;
      kx = kxt + field.normal(order) / double(order);
      ky = kyt - field.skew(order)   / double(order);
    }
    
    return (x * kx - y * ky);
  } else {
    return Series(0.0);
  }
}
#endif
