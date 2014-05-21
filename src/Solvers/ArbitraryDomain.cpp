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

#ifdef HAVE_SAAMG_SOLVER
#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>

#include "ArbitraryDomain.h"

ArbitraryDomain::ArbitraryDomain(
	BoundaryGeometry * bgeom, 
	Vector_t nr, 
	Vector_t hr, 
	std::string interpl) {

    	bgeom_m  = bgeom;
    	Geo_mincoords_m = bgeom->getmincoords();
    	Geo_maxcoords_m = bgeom->getmaxcoords();
   
    	setNr(nr);
    	setHr(hr);
    	setMinMaxZ(Geo_mincoords_m[2],Geo_maxcoords_m[2]);
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

void ArbitraryDomain::Compute(Vector_t hr, NDIndex<3> localId) {

    setHr(hr);

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

// Conversion from (x,y,z) to index in xyz plane
inline int ArbitraryDomain::toCoordIdx(int idx, int idy, int idz) {
	return (idz * nr[1] + idy) * nr[0]  + idx;
}

// Conversion from (x,y,z) to index on the 3D grid
int ArbitraryDomain::getIdx(int idx, int idy, int idz) {
	if(isInside(idx, idy, idz) && idx>=0 && idy >=0 && idz >=0 )
       	       	 	return IdxMap[toCoordIdx(idx, idy, idz)];
    	else
        	return -1;
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
 /* Expensive computation to check */
 /*
	Vector_t P0, P;
    	P0 = Vector_t(0,0,0);

	P[0] = (idx - (nr[0]-1)/2.0)*hr[0];
	P[1] = (idy - (nr[1]-1)/2.0)*hr[1];
	P[2] = (idz - (nr[2]-1)/2.0)*hr[z];
	return (bgeom_m->fastIsInside(P0, P) % 2 == 0); 
 */
    bool ret = false;
    double cx = (idx - (nr[0]-1)/2.0)*hr[0];
    double cy = (idy - (nr[1]-1)/2.0)*hr[1];
    double cz = (idz - (nr[2]-1)/2.0)*hr[2];

    int    countHz, countLz, countHy, countLy, countHx, countLx;
    std::multimap < std::tuple<int, int, int>, double >::iterator itrHz, itrLz;
    std::multimap < std::tuple<int, int, int>, double >::iterator itrHy, itrLy;
    std::multimap < std::tuple<int, int, int>, double >::iterator itrHx, itrLx;

    std::tuple<int, int, int> coordxyz(idx, idy, idz);
             
    //check if z is inside with x,y coords
    itrHz = IntersectHiZ.find(coordxyz);
    itrLz = IntersectLoZ.find(coordxyz);

    countHz = IntersectHiZ.count(coordxyz);
    countLz = IntersectLoZ.count(coordxyz);

     //check if y is inside with x,z coords
    itrHy = IntersectHiY.find(coordxyz);
    itrLy = IntersectLoY.find(coordxyz);

    countHy = IntersectHiY.count(coordxyz);
    countLy = IntersectLoY.count(coordxyz);

    //check if x is inside with y,z coords
    itrHx = IntersectHiX.find(coordxyz);
    itrLx = IntersectLoX.find(coordxyz);

    countHx = IntersectHiX.count(coordxyz);
    countLx = IntersectLoX.count(coordxyz);

    if(countHz == 1 && countLz == 1)
        ret = (cz <= itrHz->second) && (cz >= itrLz->second);
    else if(countHy == 1 && countLy == 1)
        ret = ret && (cy <= itrHy->second) && (cy >= itrLy->second);
    else if(countHx == 1 && countLx == 1)
        ret = ret && (cx <= itrHx->second) && (cx >= itrLx->second);
       
    return ret; 
}

int ArbitraryDomain::getNumXY(int z) {
    
	return numXY[z];
}


void ArbitraryDomain::getBoundaryStencil(int idxyz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {
    int x = 0, y = 0, z = 0;

    getCoord(idxyz, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);
}

void ArbitraryDomain::getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;
   // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod){
    	case CONSTANT:
        	ConstantInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        	break;
    	case LINEAR:
        //	LinearInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        	break;
    	case QUADRATIC:
	    //  QuadraticInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        	break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void ArbitraryDomain::ConstantInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    W = -1/(hr[0]*hr[0]);
    E = -1/(hr[0]*hr[0]);
    N = -1/(hr[1]*hr[1]);
    S = -1/(hr[1]*hr[1]);
    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);
    C = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    if(!isInside(x+1,y,z)) 
        E = 0.0;

    if(!isInside(x-1,y,z))
        W = 0.0;

    if(!isInside(x,y+1,z))
        N = 0.0;

    if(!isInside(x,y-1,z)) 
        S = 0.0;
    
    if(!isInside(x,y,z-1)) 
	F = 0.0;	

    if(!isInside(x,y,z+1)) 
	B = 0.0;	

}
/*
void ArbitraryDomain::LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) 
{

    scaleFactor = 1.0;

    double dx=-1, dy=-1;

    double cx = x * hr[0];
    double cy = y * hr[1];

    std::multimap< std::pair<int, int>, double >::iterator it;
    std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;

    std::pair<int, int> coordyz(y, z);
    ret = IntersectXDir.equal_range(coordyz);
    for(it = ret.first; it != ret.second; ++it) {
        if(fabs(it->second - cx) < hr[0]) {
            dx = it->second;
            break;
        }
    }

    std::pair<int, int> coordxz(x, z);
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

    int x = 0, y = 0, z = 0;

    getCoord(id, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);
}

void ArbitraryDomain::getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B) {

    if(x > 0)
        W = getIdx(x - 1, y, z);
    else
        W = -1;
    if(x < nr[0] - 1)
        E = getIdx(x + 1, y, z);
    else
        E = -1;

    if(y < nr[1] - 1)
        N = getIdx(x, y + 1, z);
    else
        N = -1;
    if(y > 0)
        S = getIdx(x, y - 1, z);
    else
        S = -1;

    if(z > 0)
        F = getIdx(x, y, z - 1);
    else
        F = -1;
    if(z < nr[2] - 1)
        B = getIdx(x, y, z + 1);
    else
        B = -1;

}


inline void ArbitraryDomain::crossProduct(double A[], double B[], double C[]) {
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
}

#endif //#ifdef HAVE_SAAMG_SOLVER
