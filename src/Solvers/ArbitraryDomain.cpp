// ------------------------------------------------------------------------
// $Version: 1.2.1 $
// ------------------------------------------------------------------------
// Copyright & License: See Copyright.readme in src directory
// ------------------------------------------------------------------------
// Class ArbitraryDomain
//   Interface to iterative solver and boundary geometry
//   for space charge calculation
//
// ------------------------------------------------------------------------
// $Author: kaman $
// $Date: 2014 $
// ------------------------------------------------------------------------
//#define DEBUG_INTERSECT_RAY_BOUNDARY

#ifdef HAVE_SAAMG_SOLVER
#include "Solvers/ArbitraryDomain.h"

#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>

#include <math.h>
#define MIN2(a,b) (((a) < (b)) ? (a) : (b))
#define MAX2(a,b) (((a) > (b)) ? (a) : (b))

ArbitraryDomain::ArbitraryDomain(
	BoundaryGeometry * bgeom,
	Vector_t nr,
	Vector_t hr,
	std::string interpl) {

    	bgeom_m  = bgeom;
    	intersectMinCoords_m = bgeom->getmincoords();
    	intersectMaxCoords_m = bgeom->getmaxcoords();

    	setNr(nr);
    	setHr(hr);

   	startId = 0;

	if(interpl == "CONSTANT")
        	interpolationMethod = CONSTANT;
	else if(interpl == "LINEAR")
        	interpolationMethod = LINEAR;
	else if(interpl == "QUADRATIC")
        	interpolationMethod = QUADRATIC;
}

ArbitraryDomain::~ArbitraryDomain() {
    //nothing so far
}

void ArbitraryDomain::Compute(Vector_t hr) {

    setHr(hr);
}


void ArbitraryDomain::Compute(Vector_t hr, NDIndex<3> localId) {

    setHr(hr);
    setNr(nr);
    int zGhostOffsetLeft  = (localId[2].first()== 0) ? 0 : 1;
    int zGhostOffsetRight = (localId[2].last() == nr[2] - 1) ? 0 : 1;
    int yGhostOffsetLeft  = (localId[1].first()== 0) ? 0 : 1;
    int yGhostOffsetRight = (localId[1].last() == nr[1] - 1) ? 0 : 1;
    int xGhostOffsetLeft  = (localId[0].first()== 0) ? 0 : 1;
    int xGhostOffsetRight = (localId[0].last() == nr[0] - 1) ? 0 : 1;

    hasGeometryChanged_m = true;

    IntersectLoX.clear();
    IntersectHiX.clear();
    IntersectLoY.clear();
    IntersectHiY.clear();
    IntersectLoZ.clear();
    IntersectHiZ.clear();


    //calculate intersection
    Vector_t P, dir, I;
    for (int idz = localId[2].first()-zGhostOffsetLeft; idz <= localId[2].last()+zGhostOffsetRight; idz++) {
	 P[2] = (idz - (nr[2]-1)/2.0)*hr[2];

	 for (int idy = localId[1].first()-yGhostOffsetLeft; idy <= localId[1].last()+yGhostOffsetRight; idy++) {
	     P[1] = (idy - (nr[1]-1)/2.0)*hr[1];

    	     for (int idx = localId[0].first()-xGhostOffsetLeft; idx <= localId[0].last()+xGhostOffsetRight; idx++) {
	       	  P[0] = (idx - (nr[0]-1)/2.0)*hr[0];

       			std::tuple<int, int, int> pos(idx, idy, idz);

	        	dir = Vector_t(0,0,1);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectHiZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));

	        	dir = Vector_t(0,0,-1);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectLoZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));

	        	dir = Vector_t(0,1,0);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectHiY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));

	        	dir = Vector_t(0,-1,0);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectLoY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));

	        	dir = Vector_t(1,0,0);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectHiX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));

	        	dir = Vector_t(-1,0,0);
		        if (bgeom_m->intersectRayBoundary(P, dir, I))
       	      		 IntersectLoX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
		}
	 }
     }

    //number of ghost nodes to the right
    int znumGhostNodesRight = 0;
    if(zGhostOffsetRight != 0) {
        for(int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
            for(int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
                if(isInside(idx, idy, localId[2].last() + zGhostOffsetRight))
                    znumGhostNodesRight++;
            }
        }
    }

    //number of ghost nodes to the left
    int znumGhostNodesLeft = 0;
    if(zGhostOffsetLeft != 0) {
        for(int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
            for(int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
                if(isInside(idx, idy, localId[2].first() - zGhostOffsetLeft))
                    znumGhostNodesLeft++;
            }
        }
    }

    //number of ghost nodes to the right
    int ynumGhostNodesRight = 0;
    if(yGhostOffsetRight != 0) {
        for(int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
            for(int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
                if(isInside(idx, localId[1].last() + yGhostOffsetRight, idz))
                    ynumGhostNodesRight++;
            }
        }
    }

    //number of ghost nodes to the left
    int ynumGhostNodesLeft = 0;
    if(yGhostOffsetLeft != 0) {
        for(int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
            for(int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
                if(isInside(idx, localId[1].first() - yGhostOffsetLeft, idz))
                    ynumGhostNodesLeft++;
            }
        }
    }


    //number of ghost nodes to the right
    int xnumGhostNodesRight = 0;
    if(xGhostOffsetRight != 0) {
	for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
                if(isInside(localId[0].last() + xGhostOffsetRight, idy, idz))
                    xnumGhostNodesRight++;
            }
        }
    }

    //number of ghost nodes to the left
    int xnumGhostNodesLeft = 0;
    if(xGhostOffsetLeft != 0) {
       	for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for(int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
                if(isInside(localId[0].first() - xGhostOffsetLeft, idy, idz))
                    xnumGhostNodesLeft++;
            }
        }
    }
    //xy points in z plane
    int numxy;
    int numtotalxy = 0;

    numXY.clear();

    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
	numxy =0;
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx =localId[0].first(); idx <= localId[0].last(); idx++) {
                if (isInside(idx, idy, idz))
                   numxy++;
            }
        }
        numtotalxy += numxy;
    }

    startId = 0;
    MPI_Scan(&numtotalxy, &startId, 1, MPI_INTEGER, MPI_SUM, Ippl::getComm());

    startId -= numtotalxy;

    //build up index and coord map
    IdxMap.clear();
    CoordMap.clear();

    register int id = startId - xnumGhostNodesLeft - ynumGhostNodesLeft - znumGhostNodesLeft;
     for (int idz = localId[2].first()-zGhostOffsetLeft; idz <= localId[2].last()+zGhostOffsetRight; idz++) {
    	 for (int idy = localId[1].first()-yGhostOffsetLeft; idy <= localId[1].last()+yGhostOffsetRight; idy++) {
    	     for (int idx = localId[0].first()-xGhostOffsetLeft; idx <= localId[0].last()+xGhostOffsetRight; idx++) {
	            if (isInside(idx, idy, idz)) {
                    IdxMap[toCoordIdx(idx, idy, idz)] = id;
                    CoordMap[id] = toCoordIdx(idx, idy, idz);
                    id++;
                 }
             }
         }
     }
}

void ArbitraryDomain::Compute(Vector_t hr, NDIndex<3> localId, Vector_t globalMeanR, Vektor<double, 4> globalToLocalQuaternion){

    setHr(hr);
    globalMeanR_m = globalMeanR;
    globalToLocalQuaternion_m = globalToLocalQuaternion;
    localToGlobalQuaternion_m[0] = globalToLocalQuaternion[0];
    for (int i=1; i<4; i++)
		localToGlobalQuaternion_m[i] = -globalToLocalQuaternion[i];

    int zGhostOffsetLeft  = (localId[2].first()== 0) ? 0 : 1;
    int zGhostOffsetRight = (localId[2].last() == nr[2] - 1) ? 0 : 1;
    int yGhostOffsetLeft  = (localId[1].first()== 0) ? 0 : 1;
    int yGhostOffsetRight = (localId[1].last() == nr[1] - 1) ? 0 : 1;
    int xGhostOffsetLeft  = (localId[0].first()== 0) ? 0 : 1;
    int xGhostOffsetRight = (localId[0].last() == nr[0] - 1) ? 0 : 1;

    hasGeometryChanged_m = true;

    IntersectLoX.clear();
    IntersectHiX.clear();
    IntersectLoY.clear();
    IntersectHiY.clear();
    IntersectLoZ.clear();
    IntersectHiZ.clear();

    //calculate intersection
    Vector_t P, saveP, dir, I;
    // TODO: Find and set the reference point for any time of geometry
    //Reference Point inside the boundary for IsoDar Geo
    Vector_t P0 = Vector_t(0,0,bgeom_m->getmincoords()[2]+hr[2]);
    //Reference Point where the boundary geometry is cylinder
    //P0 = Vector_t(0,0,0);  // Uncomment for cylinder Benchmarking -DW
    for (int idz = localId[2].first()-zGhostOffsetLeft; idz <= localId[2].last()+zGhostOffsetRight; idz++) {
	 saveP[2] = (idz - (nr[2]-1)/2.0)*hr[2];
	 for (int idy = localId[1].first()-yGhostOffsetLeft; idy <= localId[1].last()+yGhostOffsetRight; idy++) {
	     saveP[1] = (idy - (nr[1]-1)/2.0)*hr[1];
    	     for (int idx = localId[0].first()-xGhostOffsetLeft; idx <= localId[0].last()+xGhostOffsetRight; idx++) {
	       	  saveP[0] = (idx - (nr[0]-1)/2.0)*hr[0];
		  P = saveP;
		  rotateWithQuaternion(P, localToGlobalQuaternion_m);
		  P += globalMeanR_m;

		  if (bgeom_m->fastIsInside(P0, P) % 2 == 0) {
		     P0 = P;

                     std::tuple<int, int, int> pos(idx, idy, idz);

		     rotateZAxisWithQuaternion(dir, localToGlobalQuaternion_m);
		     if (bgeom_m->intersectRayBoundary(P, dir, I)) {
//			setYRangeMax(MIN2(intersectMaxCoords_m(2),I[2]));
			I -= globalMeanR_m;
			rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		IntersectHiZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
		     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "zdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		     }

		     if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
//			setYRangeMin(MAX2(intersectMinCoords_m(2),I[2]));
		        I -= globalMeanR_m;
			rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		IntersectLoZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
		     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "zdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
	  	     }

	             rotateYAxisWithQuaternion(dir, localToGlobalQuaternion_m);
		     if (bgeom_m->intersectRayBoundary(P, dir, I)) {
//			 setZRangeMax(MIN2(intersectMaxCoords_m(1),I[1]));
			 I -= globalMeanR_m;
			 rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		 IntersectHiY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
	   	     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "ydir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		     }

		     if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
//	    	        setZRangeMin(MAX2(intersectMinCoords_m(1),I[1]));
		   	I -= globalMeanR_m;
			rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		IntersectLoY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
		     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "ydir=-1" << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		     }

	             rotateXAxisWithQuaternion(dir, localToGlobalQuaternion_m);
		     if (bgeom_m->intersectRayBoundary(P, dir, I)) {
//			setXRangeMax(MIN2(intersectMaxCoords_m(0),I[0]));
			I -= globalMeanR_m;
			rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		IntersectHiX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
		     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "xdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		     }

		     if (bgeom_m->intersectRayBoundary(P, -dir, I)){
//			setXRangeMin(MAX2(intersectMinCoords_m(0),I[0]));
			I -= globalMeanR_m;
			rotateWithQuaternion(I, globalToLocalQuaternion_m);
       	      		IntersectLoX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
		     } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "xdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		     }
		  } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
			   *gmsg << "OUTSIDE" << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
		  }
	     }
	 }
     }

    IdxMap.clear();
    CoordMap.clear();

    register int id=0;
    for (int idz = 0; idz < nr[2]; idz++) {
        for (int idy = 0; idy < nr[1]; idy++) {
            for (int idx = 0; idx < nr[0]; idx++) {
                IdxMap[toCoordIdx(idx, idy, idz)] = id;
                CoordMap[id++] = toCoordIdx(idx, idy, idz);
            }
         }
     }
}

// Conversion from (x,y,z) to index in xyz plane
inline int ArbitraryDomain::toCoordIdx(int idx, int idy, int idz) {
    return (idz * nr[1] + idy) * nr[0]  + idx;
}

// Conversion from (x,y,z) to index on the 3D grid
int ArbitraryDomain::getIdx(int idx, int idy, int idz) {
    return IdxMap[toCoordIdx(idx, idy, idz)];
}

// Conversion from a 3D index to (x,y,z)
inline void ArbitraryDomain::getCoord(int idxyz, int &idx, int &idy, int &idz) {

    int id = CoordMap[idxyz];
    idx = id % (int)nr[0];
    id /= nr[0];
    idy = id % (int)nr[1];
    id /= nr[1];
    idz = id;
}

inline bool ArbitraryDomain::isInside(int idx, int idy, int idz) {
    Vector_t P;

    P[0] = (idx - (nr[0]-1)/2.0)*hr[0];
    P[1] = (idy - (nr[1]-1)/2.0)*hr[1];
    P[2] = (idz - (nr[2]-1)/2.0)*hr[2];

    bool ret = false;
    int  countH, countL;
    std::multimap < std::tuple<int, int, int>, double >::iterator itrH, itrL;
    std::tuple<int, int, int> coordxyz(idx, idy, idz);

    //check if z is inside with x,y coords
    itrH = IntersectHiZ.find(coordxyz);
    itrL = IntersectLoZ.find(coordxyz);

    countH = IntersectHiZ.count(coordxyz);
    countL = IntersectLoZ.count(coordxyz);
    if(countH == 1 && countL == 1)
        ret = (P[2] <= itrH->second) && (P[2] >= itrL->second);

     //check if y is inside with x,z coords
    itrH = IntersectHiY.find(coordxyz);
    itrL = IntersectLoY.find(coordxyz);

    countH = IntersectHiY.count(coordxyz);
    countL = IntersectLoY.count(coordxyz);
    if(countH == 1 && countL == 1)
        ret = ret && (P[1] <= itrH->second) && (P[1] >= itrL->second);

    //check if x is inside with y,z coords
    itrH = IntersectHiX.find(coordxyz);
    itrL = IntersectLoX.find(coordxyz);

    countH = IntersectHiX.count(coordxyz);
    countL = IntersectLoX.count(coordxyz);
    if(countH == 1 && countL == 1)
        ret = ret && (P[0] <= itrH->second) && (P[0] >= itrL->second);

    return ret;
}

int ArbitraryDomain::getNumXY(int z) {

	return numXY[z];
}


void ArbitraryDomain::getBoundaryStencil(int idxyz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {
    int idx = 0, idy = 0, idz = 0;

    getCoord(idxyz, idx, idy, idz);
    getBoundaryStencil(idx, idy, idz, W, E, S, N, F, B, C, scaleFactor);
}

void ArbitraryDomain::getBoundaryStencil(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;
   // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod){
    	case CONSTANT:
        	ConstantInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        	break;
    	case LINEAR:
            //  LinearInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        	break;
    	case QUADRATIC:
            //  QuadraticInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        	break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void ArbitraryDomain::ConstantInterpolation(int idx, int idy, int idz, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    W = -1/(hr[0]*hr[0]);
    E = -1/(hr[0]*hr[0]);
    N = -1/(hr[1]*hr[1]);
    S = -1/(hr[1]*hr[1]);
    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);
    C = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    if(!isInside(idx-1,idy,idz))
        W = 0.0;
    if(!isInside(idx+1,idy,idz))
        E = 0.0;

    if(!isInside(idx,idy+1,idz))
        N = 0.0;
    if(!isInside(idx,idy-1,idz))
        S = 0.0;

    if(!isInside(idx,idy,idz-1))
	F = 0.0;
    if(!isInside(idx,idy,idz+1))
	B = 0.0;
}

/*
void ArbitraryDomain::LinearInterpolation(int idx, int idy, int idz, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor)
{

    scaleFactor = 1.0;

    double dx=-1, dy=-1, dz=-1;

    double cx = (idx - (nr[0]-1)/2.0)*hr[0];
    double cy = (idy - (nr[1]-1)/2.0)*hr[1];
    double cz = (idz - (nr[2]-1)/2.0)*hr[2];

    int    countH, countL;
    std::multimap < std::tuple<int, int, int>, double >::iterator itrH, itrL;

    std::tuple<int, int, int> coordxyz(idx, idy, idz);

    //check if z is inside with x,y coords
    itrH = IntersectHiZ.find(coordxyz);
    itrL = IntersectLoZ.find(coordxyz);

    countH = IntersectHiZ.count(coordxyz);
    countL = IntersectLoZ.count(coordxyz);

    if(countH == 1 && countL == 1)
        ret = (cz <= itrH->second) && (cz >= itrL->second);



    std::multimap< std::pair<int, int>, double >::iterator it;
    std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;

    std::pair<int, int> coordyz(idy, idz);
    ret = IntersectXDir.equal_range(coordyz);
    for(it = ret.first; it != ret.second; ++it) {
        if(fabs(it->second - cx) < hr[0]) {
            dx = it->second;
            break;
        }
    }

    std::pair<int, int> coordxz(idx, idz);
    ret = IntersectYDir.equal_range(coordxz);
    for(it = ret.first; it != ret.second; ++it) {
        if(fabs(it->second - cy) < hr[1]) {
            dy = it->second;
            break;
        }
    }


    double dw=hr[0];
    double de=hr[0];
    double dn=hr[1];
    double ds=hr[1];
    C = 0.0;

    // we are a right boundary point
    if(dx >= 0 && dx > cx)
    {
        C += 1/((dx-cx)*de);
        E = 0.0;
    } else {
        C += 1/(de*de);
        E = -1/(de*de);
    }

    // we are a left boundary point
    if(dx <= 0 && dx < cx)
    {
        C += 1/((cx-dx)*dw);
        W = 0.0;
    } else {
        C += 1/(dw*dw);
        W = -1/(dw*dw);
    }

    // we are a upper boundary point
    if(dy >= 0 && dy > cy)
    {
        C += 1/((dy-cy)*dn);
        N = 0.0;
    } else {
        C += 1/(dn*dn);
        N = -1/(dn*dn);
    }

    // we are a lower boundary point
    if(dy <= 0 && dy < cy)
    {
        C += 1/((cy-dy)*ds);
        S = 0.0;
    } else {
        C += 1/(ds*ds);
        S = -1/(ds*ds);
    }

    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);

    //XXX: In stand-alone only Dirichlet for validation purposes
    if(z == 0 || z == nr[2]-1) {

        // Dirichlet
        C += 2/hr[2]*1/hr[2];

        //C += 1/hr[2]*1/hr[2];

        // case where we are on the Neumann BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        //for test no neumann
        //C += 2/((hr[2]*nr[2]/2.0) * hr[2]);
        //
        //   double d = hr[2]*(nr[2])/2;
        //   C += 2/(d * hr[2]);


        ////neumann stuff
        //W /= 2.0;
        //E /= 2.0;
        //N /= 2.0;
        //S /= 2.0;
        //C /= 2.0;

        scaleFactor *= 0.5;

    } else
        C += 2*1/hr[2]*1/hr[2];
}
*/

void ArbitraryDomain::getNeighbours(int id, int &W, int &E, int &S, int &N, int &F, int &B) {

    int idx = 0, idy = 0, idz = 0;

    getCoord(id, idx, idy, idz);
    getNeighbours(idx, idy, idz, W, E, S, N, F, B);
}

void ArbitraryDomain::getNeighbours(int idx, int idy, int idz, int &W, int &E, int &S, int &N, int &F, int &B) {

    W = getIdx(idx - 1, idy, idz);
    E = getIdx(idx + 1, idy, idz);
    N = getIdx(idx, idy + 1, idz);
    S = getIdx(idx, idy - 1, idz);
    F = getIdx(idx, idy, idz - 1);
    B = getIdx(idx, idy, idz + 1);

    if(!isInside(idx+1,idy,idz))
        E = -1;

    if(!isInside(idx-1,idy,idz))
        W = -1;

    if(!isInside(idx,idy+1,idz))
        N = -1;

    if(!isInside(idx,idy-1,idz))
        S = -1;

    if(!isInside(idx,idy,idz-1))
	F = -1;

    if(!isInside(idx,idy,idz+1))
	B = -1;

}


inline void ArbitraryDomain::crossProduct(double A[], double B[], double C[]) {
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
}

inline void ArbitraryDomain::rotateWithQuaternion(Vector_t & v, Vektor<double, 4> const quaternion) {
    // rotates a Vector_t (3 elements) using a quaternion.
    // Flip direction of rotation by quaternionVectorcomponent *= -1

    Vector_t const quaternionVectorComponent = Vector_t(quaternion(1), quaternion(2), quaternion(3));
    double const quaternionScalarComponent = quaternion(0);

    v = 2.0 * dot(quaternionVectorComponent, v) * quaternionVectorComponent
        + (quaternionScalarComponent * quaternionScalarComponent
        -  dot(quaternionVectorComponent, quaternionVectorComponent)) * v
        + 2.0 * quaternionScalarComponent * cross(quaternionVectorComponent, v);
}

inline void ArbitraryDomain::rotateXAxisWithQuaternion(Vector_t & v, Vektor<double, 4> const quaternion) {
    // rotates the positive xaxis using a quaternion.

    v(0) = quaternion(0) * quaternion(0)
         + quaternion(1) * quaternion(1)
         - quaternion(2) * quaternion(2)
         - quaternion(3) * quaternion(3);

    v(1) = 2.0 * (quaternion(1) * quaternion(2) + quaternion(0) * quaternion(3));
    v(2) = 2.0 * (quaternion(1) * quaternion(3) - quaternion(0) * quaternion(2));
}

inline void ArbitraryDomain::rotateYAxisWithQuaternion(Vector_t & v, Vektor<double, 4> const quaternion) {
    // rotates the positive yaxis using a quaternion.

    v(0) = 2.0 * (quaternion(1) * quaternion(2) - quaternion(0) * quaternion(3));

    v(1) = quaternion(0) * quaternion(0)
         - quaternion(1) * quaternion(1)
         + quaternion(2) * quaternion(2)
         - quaternion(3) * quaternion(3);

    v(2) = 2.0 * (quaternion(2) * quaternion(3) + quaternion(0) * quaternion(1));
}

inline void ArbitraryDomain::rotateZAxisWithQuaternion(Vector_t & v, Vektor<double, 4> const quaternion) {
    // rotates the positive zaxis using a quaternion.
    v(0) = 2.0 * (quaternion(1) * quaternion(3) + quaternion(0) * quaternion(2));
    v(1) = 2.0 * (quaternion(2) * quaternion(3) - quaternion(0) * quaternion(1));

    v(2) = quaternion(0) * quaternion(0)
         - quaternion(1) * quaternion(1)
         - quaternion(2) * quaternion(2)
         + quaternion(3) * quaternion(3);

}
#endif //#ifdef HAVE_SAAMG_SOLVER
