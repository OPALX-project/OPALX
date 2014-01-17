#ifdef HAVE_SAAMG_SOLVER
#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>

#include "ArbitraryDomain.h"


//IFF: simplified case where we have intersect = 2
//IFF: every node counts own number of gridpoints and then sum to global procs before MyPID (including ghost layers in z direction).
//IFF: triangles intersecting in z dir? (cathode?) -> boundary condition there? neumann? then we should take a normal and stuff!?

ArbitraryDomain::ArbitraryDomain(BoundaryGeometry * bgeom, Vector_t Geo_mincoords, Vector_t Geo_maxcoords, Vector_t nr, Vector_t hr, NDIndex<3> locidx, std::string interpl) {
    bgeom_m=bgeom;
    Geo_mincoords_m = bgeom->getmincoords();
    Geo_maxcoords_m = bgeom->getmaxcoords();
    setNr(nr);
    setHr(hr);
    localidx = locidx;
    startIdx = 0;

    std::cout << "Arbitrary domain TK localidx:" << localidx << std::endl; 
    std::cout << Geo_mincoords_m << std::endl;
    std::cout << Geo_maxcoords_m << std::endl;
    std::cout << "bgeom-maxcoords" << bgeom->getmaxcoords() << std::endl;	
    std::cout << "hr_m" << hr << std::endl;	

    if(interpl == "constant")
        interpolationMethod = CONSTANT;
    else if(interpl == "linear")
        interpolationMethod = LINEAR;
    else if(interpl == "quadratic")
        interpolationMethod = QUADRATIC;

}

ArbitraryDomain::~ArbitraryDomain() {
    //nothing so far
}

void ArbitraryDomain::Compute(Vector_t hr) {

    setHr(hr);

    hasGeometryChanged_m = true;

    //h5_float64_t edge1, edge2, t, q, p, P0;
    double edge1[3], edge2[3], t[3], q[3], p[3], P0[3];
    double det, invDet, u, v, tt;
    double dir[3], origin[3];
    int zGhostOffsetLeft = (localidx[2].first() == 0) ? 0 : 1;
    int zGhostOffsetRight = (localidx[2].last() == nr[2] - 1) ? 0 : 1;
    std::multimap< std::pair<int, int>, double >::iterator it;
    std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;
    bool hit;


    //forall triangles
    int numbfaces = bgeom_m->getNumBFaces();
    for (int i = 0; i < numbfaces ; i++) {

        Vector_t x1 = bgeom_m->getVertexCoord(i, 1); //geo3Dcoords_m[allbfaces_m[4 * face_id + vertex_id]]
        Vector_t x2 = bgeom_m->getVertexCoord(i, 2); //geo3Dcoords_m[allbfaces_m[4 * face_id + vertex_id]]
        Vector_t x3 = bgeom_m->getVertexCoord(i, 3); //geo3Dcoords_m[allbfaces_m[4 * face_id + vertex_id]]

	P0[0] = x2[0];
        P0[1] = x2[1];
        P0[2] = x2[2];

        edge1[0] = x1[0] - x2[0];
        edge1[1] = x1[1] - x2[1];
        edge1[2] = x1[2] - x2[2];

        edge2[0] = x3[0] - x2[0];
        edge2[1] = x3[1] - x2[1];
        edge2[2] = x3[2] - x2[2];

        //IFF: use normed dir -> t = intersection (?)
        //x dir rays = forall points in origin y,z with dir (1,0,0) (Y,Z -> COORD)
        dir[0] = 1;
        dir[1] = 0;
        dir[2] = 0;
        origin[0] = 0; //Geo_mincoords_m(0); //0;
        origin[1] = 0; //Geo_mincoords_m(1); //0;
	origin[2] = getMinZ(); //0;
        for(int y = 0; y < nr[1]; y++) {
            origin[1] = (y - (nr[1]-1)/2.0)*hr[1];
            //ONE LAYER MORE = GHOST LAYER
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].last() + zGhostOffsetRight; z++) {

                origin[2] = getMinZ()+z * hr[2];

                crossProduct(dir, edge2, p);
                det = dotProduct(edge1, p);

                if(det > -1e-5 && det < 1e-5)
                    continue;

                invDet = 1.0 / det;

                t[0] = origin[0] - P0[0];
                t[1] = origin[1] - P0[1];
                t[2] = origin[2] - P0[2];
                u = dotProduct(t, p) * invDet;
                if(u < 0.0 || u > 1.0)
                    continue;

                crossProduct(t, edge1, q);
                v = dotProduct(dir, q) * invDet;
                if(v < 0.0 || u + v > 1.0)
                    continue;

                //ray really intersects triangle
                tt = dotProduct(edge2, q) * invDet;

                //put intersection point in structure
                std::pair<int, int> tmp(y, z);
                //IntersectZDir[tmp] = tt;
                hit = false;
                ret = IntersectXDir.equal_range(tmp);
                for(it = ret.first; it != ret.second; ++it)
                    hit = hit || (fabs(it->second - tt) < 1e-15);
                if(!hit)
                    IntersectXDir.insert(std::pair< std::pair<int, int>, double >(tmp, tt));

            }
        }

        //y dir rays = forall points in origin x,z with dir (0,1,0) (X,Z -> COORD)
        dir[0] = 0;
        dir[1] = 1;
        dir[2] = 0;
	origin[0] = 0;//Geo_mincoords_m(0); //0;
        origin[1] = 0;//Geo_mincoords_m(1); //0;
        origin[2] = getMinZ(); //0;
        for(int x = 0; x < nr[0]; x++) {
            origin[0] = (x - (nr[0]-1)/2.0)*hr[0];
            //ONE LAYER MORE = GHOST LAYER
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].last() + zGhostOffsetRight; z++) {

                origin[2] = getMinZ() + z * hr[2];

                crossProduct(dir, edge2, p);
                det = dotProduct(edge1, p);

                if(det > -1e-5 && det < 1e-5)
                    continue;

                invDet = 1.0 / det;

                t[0] = origin[0] - P0[0];
                t[1] = origin[1] - P0[1];
                t[2] = origin[2] - P0[2];
                u = dotProduct(t, p) * invDet;
                if(u < 0.0 || u > 1.0)
                    continue;

                crossProduct(t, edge1, q);
                v = dotProduct(dir, q) * invDet;
                if(v < 0.0 || u + v > 1.0)
                    continue;

                //ray really intersects triangle
                tt = dotProduct(edge2, q) * invDet;

                //put intersection point in structure
                std::pair<int, int> tmp(x, z);
                hit = false;
                ret = IntersectYDir.equal_range(tmp);
                for(it = ret.first; it != ret.second; ++it)
                    hit = hit || (fabs(it->second - tt) < 1e-15);
                if(!hit)
                    IntersectYDir.insert(std::pair< std::pair<int, int>, double >(tmp, tt));

            }
        }

        //z dir rays = forall points in origin x,y with dir (0,0,1) (X,Y -> COORD)
        dir[0] = 0;
        dir[1] = 0;
        dir[2] = 1;
        origin[0] = 0; //Geo_mincoords_m(0); //0;
        origin[1] = 0; //Geo_mincoords_m(1); //0;
        origin[2] = getMinZ(); //0;
        for(int x = 0; x < nr[0]; x++) {
            origin[0] = (x - (nr[0]-1)/2.0)*hr[0];
            for(int y = 0; y < nr[1]; y++) {
            	origin[1] = (y - (nr[1]-1)/2.0)*hr[1];

                crossProduct(dir, edge2, p);
                det = dotProduct(edge1, p);

                if(det > -1e-5 && det < 1e-5)
                    continue;

                invDet = 1.0 / det;

                t[0] = origin[0] - P0[0];
                t[1] = origin[1] - P0[1];
                t[2] = origin[2] - P0[2];
                u = dotProduct(t, p) * invDet;
                if(u < 0.0 || u > 1.0)
                    continue;

                crossProduct(t, edge1, q);
                v = dotProduct(dir, q) * invDet;
                if(v < 0.0 || u + v > 1.0)
                    continue;

                //ray really intersects triangle
                tt = dotProduct(edge2, q) * invDet;

                //put intersection point in structure
                std::pair<int, int> tmp(x, y);
                hit = false;
                ret = IntersectZDir.equal_range(tmp);
                for(it = ret.first; it != ret.second; ++it)
                    hit = hit || (fabs(it->second - tt) < 1e-15);
                if(!hit)
                    IntersectZDir.insert(std::pair< std::pair<int, int>, double >(tmp, tt));
            }
        }

    }
    std::cout << "number of intersections in x-dir: " << IntersectXDir.size() << std::endl;
    std::cout << "number of intersections in y-dir: " << IntersectYDir.size() << std::endl;
    std::cout << "number of intersections in z-dir: " << IntersectZDir.size() << std::endl;

    //number of ghost nodes to the left
    int numGhostNodesLeft = 0;
    if(localidx[2].first() != 0) {
        for(int x = 0; x < nr[0]; x++) {
            for(int y = 0; y < nr[1]; y++) {

                if(isInside(x, y, localidx[2].first() - zGhostOffsetLeft))
                    numGhostNodesLeft++;

            }
        }
    }

    std::cout << "ghost nodes left: " << numGhostNodesLeft << std::endl;

    //xy points in z plane
    int numxy = 0;
    int numtotal = 0;
    numXY.clear();
    for(int idx = localidx[2].first(); idx <= localidx[2].last(); idx++) {
        numxy = 0;
        for(int x = 0; x < nr[0]; x++) {
            for(int y = 0; y < nr[1]; y++) {

                if(isInside(x, y, idx))
                    numxy++;

            }
        }

        numXY[idx-localidx[2].first()] = numxy;
        numtotal += numxy;

    }

    std::cout << "number of gridpoints: " << numtotal << std::endl;

    startIdx = 0;
    MPI_Scan(&numtotal, &startIdx, 1, MPI_INTEGER, MPI_SUM, Ippl::getComm());
    startIdx -= numtotal;

    std::cout << "start idx: " << startIdx << std::endl;

    //build up index and coord map
    IdxMap.clear();
    CoordMap.clear();
    register int idx = startIdx - numGhostNodesLeft;

	std::cout << idx << std::endl; 
    for(int x = 0; x < nr[0]; x++) {
        for(int y = 0; y < nr[1]; y++) {
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].last() + zGhostOffsetRight; z++) {

                if(isInside(x, y, z)) {
                    IdxMap[toCoordIdx(x, y, z)] = idx;
                    CoordMap[idx++] = toCoordIdx(x, y, z);
                }

            }
        }
    }

}

/// conversion from (x,y,z) to index in xyz plane
inline int ArbitraryDomain::toCoordIdx(int x, int y, int z) {
    return (z * nr[1] + y) * nr[0] + x;
}

/// conversion from (x,y,z) to index on the 3D grid
/*inline*/
int ArbitraryDomain::getIdx(int x, int y, int z) {
    if(isInside(x, y, z) && x >= 0 && y >= 0 && z >= 0)
        return IdxMap[toCoordIdx(x, y, z)];
    else
        return -1;
}


/// conversion from a 3D index to (x,y,z)
inline void ArbitraryDomain::getCoord(int idx, int &x, int &y, int &z) {

    int idxx = CoordMap[idx];

    x = idxx % (int)nr[0];
    idxx /= nr[0];
    y = idxx % (int)nr[1];
    idxx /= nr[1];
    z = idxx;

}


//IFF: at the moment this implementation only allows 2 intersections with geometry!!
inline bool ArbitraryDomain::isInside(int x, int y, int z) {

    bool ret = false;
    //double cx = x * hr[0];
    //double cy = y * hr[1];
    double cx = (x - (nr[0]-1)/2)*hr[0];
    double cy = (y - (nr[1]-1)/2)*hr[1];
 //   double cz = z * hr[2];
    double val1 = 0;
    double val2 = 0;
    std::multimap < std::pair<int, int>, double >::iterator itr;

    //check if x is inside with y,z coords
    std::pair<int, int> coordyz(y, z);
    itr = IntersectXDir.find(coordyz);
    int count = IntersectXDir.count(coordyz);
    if(count == 1)
        ret = (cx == itr->second);
    else {
        val1 = itr->second;
        val2 = (++itr)->second;
        ret = (cx <= std::max(val1, val2)) && (cx >= std::min(val1, val2));
    }

    //check if y is inside with x,z coords
    std::pair<int, int> coordxz(x, z);
    itr = IntersectYDir.find(coordxz);
    count = IntersectYDir.count(coordxz);
    if(count == 1)
        ret = ret && (cy == itr->second);
    else {
        val1 = itr->second;
        val2 = (++itr)->second;
        ret = ret && (cy <= std::max(val1, val2)) && (cy >= std::min(val1, val2));
    }

    //check if z is inside with x,y coords
/*
    std::pair<int,int> coordxy(x,y);
    itr = IntersectZDir.find(coordxy);
    count = IntersectZDir.count(coordxy);
    if(count == 1)
      ret = ret && (cz == itr->second);
    else {
      val1 = itr->second;
      val2 = (++itr)->second;
      ret = ret && (cz <= std::max(val1, val2)) && (cz >= std::min(val1, val2));
    }

    ret = ret && (cz >= 0 && cz <= nr[2]);
*/
    return ret;

}

int ArbitraryDomain::getNumXY(int z) {

    return numXY[z];

}

void ArbitraryDomain::getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

   // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod) {

    case CONSTANT:
        ConstantInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        break;

    case LINEAR:
        LinearInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        break;

    case QUADRATIC:
	std::cout << "TK implement this later" <<std::endl;
    //    QuadraticInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor);
        break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void ArbitraryDomain::ConstantInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    scaleFactor = 1.0;

    W = -1/(hr[0]*hr[0]);
    E = -1/(hr[0]*hr[0]);
    N = -1/(hr[1]*hr[1]);
    S = -1/(hr[1]*hr[1]);
    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);
    C = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    // we are a right boundary point
    if(!isInside(x+1,y,z)) {
        E = 0.0;
    }

    // we are a left boundary point
    if(!isInside(x-1,y,z)) {
        W = 0.0;
    }

    // we are a upper boundary point
    if(!isInside(x,y+1,z)) {
        N = 0.0;
    }

    // we are a lower boundary point
    if(!isInside(x,y-1,z)) {
        S = 0.0;
    }

    if(z == 0 || z == nr[2] - 1) {
        // case where we are on the Robin BC in Z-direction
        // where we distinguish two cases
        // IFF: this values should not matter because they
        // never make it into the discretization matrix
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        // add contribution of Robin discretization to center point
        // d the distance between the center of the bunch and the boundary
        //double cx = (x-(nr[0]-1)/2)*hr[0];
        //double cy = (y-(nr[1]-1)/2)*hr[1];
        //double cz = hr[2]*(nr[2]-1);
        //double d = sqrt(cx*cx+cy*cy+cz*cz);
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);
        //C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);

        // scale all stencil-points in z-plane with 0.5 (Robin discretization)
        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;
    }

}
/*
void ArbitraryDomain::LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) 
{

    scaleFactor = 1.0;

    double dx=-1, dy=-1;

    double cx = (x - (nr[0]-1)/2)*hr[0];
    double cy = (y - (nr[1]-1)/2)*hr[1];

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
        
        //  double d = hr[2]*(nr[2])/2;
        //   C += 2/(d * hr[2]);


        //neumann stuff
        // W /= 2.0;
        // E /= 2.0;
        // N /= 2.0;
        // S /= 2.0;
        //C /= 2.0;

       // scaleFactor *= 0.5;
       // 

    } else 
        C += 2*1/hr[2]*1/hr[2];
}
*/
void ArbitraryDomain::LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    double cx = x * hr[0];
    double cy = y * hr[1];
    //double cz = z * hr[2]; 

    std::multimap< std::pair<int, int>, double >::iterator it;
    std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;

    //since every vector here is only an array[2] we
    //can catch all cases manually
    double dx = -1.0, dy = -1.0;//, dz = -1.0; 
    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    double df = hr[2];
    double db = hr[2];
    C = 0.0;

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

/*
    std::pair<int, int> coordxy(x,y);
    ret = IntersectZDir.equal_range(coordxy);
    for(it=ret.first; it!=ret.second; ++it) {
      if(abs(it->second - cz) < hr[2]) {
        dz = it->second;
        break;
      }
    }
*/   

    //we are a right boundary point
    //if(!isInside(x+1,y,z))
    if(dx >= 0 && dx > cx)
        de = dx - cx;

    //we are a left boundary point
    //if(!isInside(x-1,y,z))
    if(dx >= 0 && dx < cx)
        dw = cx - dx;

    //we are a upper boundary point
    //if(!isInside(x,y+1,z))
    if(dy >= 0 && dy > cy)
        dn = dy - cy;

    //we are a lower boundary point
    //if(!isInside(x,y-1,z))
    if(dy >= 0 && dy < cy)
        ds = cy - dy;

    //we are a lower boundary point
/*
    //if(!isInside(x,y,z+1))
    if(dz >= 0 && dz > cz)
        df = dz - cz;

    //we are a lower boundary point
    //if(!isInside(x,y,z-1))
    if(dy >= 0 && dy <= cy)
        db = cz - dz;
 */

    //for regular gridpoints no problem with symmetry, just boundary

    if(dw != 0)
        W = -(df + db) * (dn + ds) / dw;
    else
        W = 0;
    if(de != 0)
        E = -(df + db) * (dn + ds) / de;
    else
        E = 0;
    if(dn != 0)
        N = -(df + db) * (dw + de) / dn;
    else
        N = 0;
    if(ds != 0)
        S = -(df + db) * (dw + de) / ds;
    else
        S = 0;
    if(df != 0)
        F = -(dw + de) * (dn + ds) / df;
    else
        F = 0;
    if(db != 0)
        B = -(dw + de) * (dn + ds) / db;
    else
        B = 0;

    //RHS scaleFactor for current 3D index
    //0.5* comes from discretiztaion
    //scaleFactor = 0.5*(dw+de)*(dn+ds)*(df+db);
    scaleFactor = 0.5;
    if((dw + de) != 0)
        scaleFactor *= (dw + de);
    if((dn + ds) != 0)
        scaleFactor *= (dn + ds);
    if((df + db) != 0)
        scaleFactor *= (df + db);

    //catch the case where a point lies on the boundary
    //IFF: do this more elegant!
    double m1 = dw * de;
    double m2 = dn * ds;
    if(de == 0)
        m1 = dw;
    if(dw == 0)
        m1 = de;
    if(dn == 0)
        m2 = ds;
    if(ds == 0)
        m2 = dn;
    //IFF: dn+ds || dw+de can be 0
    //C = 2*(dn+ds)*(dw+de)/hr[2];
    C = 2 / hr[2];
    if(dw != 0 || de != 0)
        C *= (dw + de);
    if(dn != 0 || ds != 0)
        C *= (dn + ds);
    if(dw != 0 || de != 0)
        C += (df + db) * (dn + ds) * (dw + de) / m1;
    if(dn != 0 || ds != 0)
        C += (df + db) * (dw + de) * (dn + ds) / m2;

    if(C <= 0)
      std::cout << "!!!!!!!!!! C is <= 0: " << C <<  std::endl;

        // handle boundary condition in z direction
    if(z == 0 || z == nr[2] - 1) {
	F = -1 / (hr[2] * hr[2]);
	B = -1 / (hr[2] * hr[2]);
	C += 2 / (hr[2] * hr[2]);

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;

        //if(C <= 0)
        //  cout << "!!!!!!!!!! C is <= 0" << endl;
    }
	
}

/*
 * OLD STENCIL O(h)
void EllipticDomain::getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C) {


  double cx = (x-floor(nr[0]/2.0))*hr[0];
  double cy = (y-floor(nr[1]/2.0))*hr[1];

  //since every vector here is only an array[2] we
  //can catch all cases manually
  double dx = 0.0;
  //std::vector <double> tmp = IntersectXDir.find(y)->second;
  multimap<int, double>::iterator it = IntersectXDir.find(y);
  if(cx < 0)
    it++;
  dx = it->second;
  //else
  //  dx = (++it)->second;

  double dy = 0.0;
  //std::vector <double> tmp = IntersectYDir.find(x)->second;
  it = IntersectYDir.find(x);
  if(cy < 0)
    it++;
  dy = it->second;
  //else
  //  dy = (++it)->second;

  double dw=hr[0];
  double de=hr[0];
  double dn=hr[1];
  double ds=hr[1];
  C = 0.0;

  bool hasW = true;
  bool hasE = true;
  bool hasN = true;
  bool hasS = true;

  //TODO: cx > 0 --> xc <= dx
  //TODO: remove isInside and replace with direct computation if we are inside

  if(true) { // z != 0 && z != nr[2]-1) {

    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    if(!isInside(x+1,y,z)) {
      //we are a right boundary point
      C += (dx-cx);
      de = 0.0;
    } else {
      C += de;
    }

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    if(!isInside(x-1,y,z)) {
      //we are a left boundary point
      C += (abs(dx)-abs(cx));
      dw = 0.0;
    } else {
      C += dw;
    }

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    if(!isInside(x,y+1,z)) {
      //we are a upper boundary point
      C += (dy-cy);
      dn = 0.0;
    } else {
      C += dn;
    }

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    if(!isInside(x,y-1,z)) {
      //we are a lower boundary point
      C += (abs(dy)-abs(cy));
      ds = 0.0;
    } else {
      C += ds;
    }

    W = -dw;
    E = -de;
    N = -dn;
    S = -ds;
    F = -hr[2];
    B = -hr[2];
    C += hr[2]+hr[2];

    if(C <= 0)
      cout << "!!!!!!!!!! C is <= 0" << endl;

    if(C <= 2*hr[2])
      cout << "!!!!!!!!!! C is smaller than surrounding nodes" << endl;

    //if(cx == dx && cy == dy) {
    //  C = 0.0;
    //}

  } else {

    //NEUMAN in Z!

    if(z == 0) {
      F = 0.0;
      B = -hr[2];
      C = 2*hr[2];
    } else {
      B = 0.0;
      F = -hr[2];
      C = 2*hr[2];
    }

    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    if(!isInside(x+1,y,z)) {
      //we are a right boundary point
      C += (dx-cx);
      de = 0.0;
    } else {
      C += de;
    }

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    if(!isInside(x-1,y,z)) {
      //we are a left boundary point
      C += (abs(dx)-abs(cx));
      dw = 0.0;
    } else {
      C += dw;
    }

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    if(!isInside(x,y+1,z)) {
      //we are a upper boundary point
      C += (dy-cy);
      dn = 0.0;
    } else {
      C += dn;
    }

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    if(!isInside(x,y-1,z)) {
      //we are a lower boundary point
      C += (abs(dy)-abs(cy));
      ds = 0.0;
    } else {
      C += ds;
    }

    //neumann stuff
    W = -dw/2.0;
    E = -de/2.0;
    N = -dn/2.0;
    S = -ds/2.0;
    C /= 2.0;

    if(C <= 0)
      cout << "!!!!!!!!!! C is <= 0" << endl;

    //if(cx == dx && cy == dy) {
    //  C = 0.0;
    //}
  }
}
*/

void ArbitraryDomain::getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    //IFF: reverse map seah5errh?
    //double mem or O(n) seah5errh in map to get x,y,z from idx

    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);

}

void ArbitraryDomain::getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B) {

    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

//IFF:should be ok
//getIdx returns -1 when not inside domain
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

    //we have ghost nodes that can handle nodes at the boundary for z-1 and z+1
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
