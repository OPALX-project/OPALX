#ifndef RECQUADMAP_HH
#define RECQUADMAP_HH
// ------------------------------------------------------------------------
// $RCSfile: GTRecquadMap.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTElements/GTRecquadMap.hh,v 1.1.1.1 2003/01/23 09:13:53 adelmann Exp $
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
// $Log: GTRecquadMap.hh,v $
// Revision 1.1.1.1  2003/01/23 09:13:53  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include "Physics/Physics.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"

#include "GTElemData.hh"
#include "Algorithms/PartData.h"

using namespace Physics;
//using namespace SIXVect;

typedef FTps<double,6> Series;

class RecQuadMap {
private:
  ElemData cfgData_m;
  PartData pdata_m;
  double r1_m;
  double r2_m;
  double scale_m;
  double h_m;
public:
  RecQuadMap(ElemData cfgData, PartData pdata)
  { 
    cfgData_m = cfgData;
    pdata_m = pdata;
    r1_m = cfgData.r1;
    r2_m = cfgData.r2;
    FTps<double,6>::setGlobalTruncOrder(cfgData.order);
        
    double multipoleScalFactor = pdata_m.getP() / Physics::c;
    scale_m = Physics::c*pdata_m.getQ()/pdata_m.getP();
  };

  ~RecQuadMap() { };

  DragtFinnMap<3> getDragtFinnMap(DragtFinnMap<3> map);
 
  FVps<double,6> getTaylorSeries(FVps<double,6> map);
  
 private:
 
  FTps<double,6> buildVectorPotential(unsigned int order);
 
  Series getH();
};


DragtFinnMap<3> RecQuadMap::getDragtFinnMap(DragtFinnMap<3> map) {
  Series H;
  H = getH();
  double parlen = cfgData_m.length;
  DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H * parlen);
  return  theMap;
  //return map.catenate(theMap);
}

FVps<double,6> RecQuadMap::getTaylorSeries(FVps<double,6> map) {
  Series  H;
  H = getH();
  double parlen = cfgData_m.length;
  return ExpMap(- getH() * parlen).substitute(map);
}

Series  RecQuadMap::getH() {
    
    cout << "r1= " << r1_m << " r2= " << r2_m << " order= " << cfgData_m.order << endl;
    
    // Build Hamiltonian in local coordinates; substitution done later.
    // Step 1: Define variables.
    double beta =  pdata_m.getBeta();
    Series H;
    Series px = Series::makeVariable(SIXVect::PX);
    Series py = Series::makeVariable(SIXVect::PY);
    
    // Step 3: Vector potential  
    Series As = buildVectorPotential(cfgData_m.order) * scale_m; 
  
    if (cfgData_m.uType==MAD9) {
    	Series pt = Series::makeVariable(SIXVect::PT) + 1.0;
    
	// Step 2: Kinematic terms.
	Series pz = sqrt(pt*pt - px*px - py*py);
	double kin = pdata_m.getM() / pdata_m.getP();
	Series E = sqrt(pt*pt + kin*kin) / beta;
        
	// Step 4: Finish Hamiltonian, build exponential map,
	//         and substitute previous map.
	H = As + E - pz;
    
    } else if (cfgData_m.uType==MARYLIE) {
    	Series pt = Series::makeVariable(SIXVect::PT);
	Series Pt = 2.0*pt/beta;
	Series Ptt = pt/beta;
	Series pz = sqrt(1.0 + Pt + pt*pt - px*px - py*py);
	// Step 4: Finish Hamiltonian, build exponential map,
	//         and substitute previous map.
	H =  As - pz + Ptt;
    }
    return H;
}

Series RecQuadMap::buildVectorPotential(unsigned int order)
{
    Series A; 
    if (order > 0) {
	static const Series Ax = Series::makeVariable(SIXVect::X);
	static const Series Ay = Series::makeVariable(SIXVect::Y);
	static const Series At = Series::makeVariable(SIXVect::TT);
    }  
    return A;
}
#endif
    /* // Construct terms constant and linear in y. */
/*     Series Ae = + field_m.normal(order); // Term even in y. */
/*     Series Ao = - field_m.skew(order);   // Term odd  in y. */
    
/*     for (int i = order; --i >= 1; ) { */
/*       Ae = Ae * x + field_m.normal(i); */
/*       Ao = Ao * x - field_m.skew(i); */
/*     } */
    
/*     Ae = + (Ae * (1.0 + h * x)).integral(SIXVect::X); */
/*     Ao = - (Ao * (1.0 + h * x)); */
    
/*     // Add terms up to maximum order. */
/*     As = Ae + y * Ao; */
/*     int k = 2; */
    
/*     if (k <= order) { */
/*       Series yp = y * y / 2.0; */
/*       Series factor = h / (1.0 + h*x); */
      
/*       while (true) { */
/* 	// Terms even in y. */
/* 	Ae = Ae.derivative(SIXVect::X); */
/* 	Ae = (factor*Ae - Ae.derivative(SIXVect::X)) * yp; */
/* 	As += Ae; */
/* 	if (++k > order) break; */
/* 	yp *= y / double(k); */
	
/* 	// Terms odd in y. */
/* 	Ao = Ao.derivative(SIXVect::X); */
/* 	Ao = (factor*Ao - Ao.derivative(SIXVect::X)) * yp; */
/* 	As += Ao; */
/* 	if (++k > order) break; */
/* 	yp *= y / double(k); */
/*       } */
/*     } */
/*   } */
