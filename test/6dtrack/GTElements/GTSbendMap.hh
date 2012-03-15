#ifndef SBENDMAP_HH
#define SBENDMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTSbendMap.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTSbendMap.hh,v 1.1.1.1 2003/01/23 09:13:53 adelmann Exp $
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
// $Log: GTSbendMap.hh,v $
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

class SBendMap {
private:
  ElemData cfgData_m;
  PartData pdata_m;
  BMultipoleField field_m;
  double scale_m;
  double h_m;
public:
  SBendMap(ElemData cfgData, PartData pdata)
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

  ~SBendMap() { };

  DragtFinnMap<3> getDragtFinnMap(DragtFinnMap<3> map);
 
  FVps<double,6> getTaylorSeries(FVps<double,6> map);
    
  FVps<double,6> getEntranceFringeMap(FVps<double,6> map);

  FVps<double,6> getExitFringeMap(FVps<double,6> map);

private:
 
  FTps<double,6> buildSBendVectorPotential(double h);
 
  Series getH();
};


DragtFinnMap<3> SBendMap::getDragtFinnMap(DragtFinnMap<3> map) {
  Series H;
  H = getH();
  double parlen = cfgData_m.length;
  DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H * parlen);
  return  theMap;
  //return map.catenate(theMap);
}

FVps<double,6> SBendMap::getTaylorSeries(FVps<double,6> map) {
  Series  H;
  H = getH();
  double parlen = cfgData_m.length;
  return ExpMap(- getH() * parlen).substitute(map);
}

FVps<double,6> SBendMap::getEntranceFringeMap(FVps<double,6> map)
{
  // *** MISSING *** Higher order terms for entrance fringe.
  double alpha = 0.0;
  if (cfgData_m.e1!=0.0) {
    double hx = scale_m * field_m.normal(1);
    double ex = hx * tan(cfgData_m.e1);
    double ey = hx * tan(cfgData_m.e1-alpha); // + map[PX][0]);
    map[SIXVect::PX] += ex * map[SIXVect::X];
    map[SIXVect::PY] -= ey * map[SIXVect::Y];
    map.substitute(map);
  }
  return map;
}

FVps<double,6> SBendMap::getExitFringeMap(FVps<double,6> map)
{
  // *** MISSING *** Higher order terms for exit fringe.
  double alpha = 0.0;
  if (cfgData_m.e2!=0.0) {
    double hx = scale_m * field_m.normal(1);
    double ex = hx * tan(cfgData_m.e2);
    double ey = hx * tan(cfgData_m.e2-alpha); // map[PX][0]);
    map[SIXVect::PX] += ex * map[SIXVect::X];
    map[SIXVect::PY] -= ey * map[SIXVect::Y];
    map.substitute(map);
  }
  return map;
}

Series  SBendMap::getH() {

  // Build Hamiltonian in local coordinates; substitution done later.
  // Step 1: Define variables.
  double beta =  pdata_m.getBeta();
  Series H;
  Series px = Series::makeVariable(SIXVect::PX);
  Series py = Series::makeVariable(SIXVect::PY);

  // Step 3: Vector potential  (1 + h*x) * As  in curved reference.
  Series As = buildSBendVectorPotential(h_m) * scale_m; 
  
  if (cfgData_m.uType==MAD9) {
    Series pt = Series::makeVariable(SIXVect::PT) + 1.0;
  
    // Step 2: Kinematic terms.
    Series x  = Series::makeVariable(SIXVect::X);
    Series pz = (1.0 + h_m*x) * sqrt(pt*pt - px*px - py*py);
    double kin = pdata_m.getM() / pdata_m.getP();
    Series E  = sqrt(pt*pt + kin*kin) / beta;
  
    // Step 4: Finish Hamiltonian
    H = As + E - pz;
  } else if (cfgData_m.uType==MARYLIE) {
    
    Series x  = Series::makeVariable(SIXVect::X);
    
    Series pt = -Series::makeVariable(SIXVect::PT);
    Series Pt = 2.0*pt/beta;
    Series Ptt = pt/beta;
    Series pz =  (1.0 + h_m*x) * sqrt(1.0 + Pt + pt*pt - px*px - py*py);
    H = As - pz + Ptt;
  }
  return H;
}

Series SBendMap::buildSBendVectorPotential(double h)
{
  Series As; 
  int order = field_m.order();

  if (order > 0) {
    static const Series x = Series::makeVariable(SIXVect::X);
    static const Series y = Series::makeVariable(SIXVect::Y);
    
    // Construct terms constant and linear in y.
    Series Ae = + field_m.normal(order); // Term even in y.
    Series Ao = - field_m.skew(order);   // Term odd  in y.
    
    for (int i = order; --i >= 1; ) {
      Ae = Ae * x + field_m.normal(i);
      Ao = Ao * x - field_m.skew(i);
    }
    
    Ae = + (Ae * (1.0 + h * x)).integral(SIXVect::X);
    Ao = - (Ao * (1.0 + h * x));
    
    // Add terms up to maximum order.
    As = Ae + y * Ao;
    int k = 2;
    
    if (k <= order) {
      Series yp = y * y / 2.0;
      Series factor = h / (1.0 + h*x);
      
      while (true) {
	// Terms even in y.
	Ae = Ae.derivative(SIXVect::X);
	Ae = (factor*Ae - Ae.derivative(SIXVect::X)) * yp;
	As += Ae;
	if (++k > order) break;
	yp *= y / double(k);
	
	// Terms odd in y.
	Ao = Ao.derivative(SIXVect::X);
	Ao = (factor*Ao - Ao.derivative(SIXVect::X)) * yp;
	As += Ao;
	if (++k > order) break;
	yp *= y / double(k);
      }
    }
  }
  return As;
}
#endif
