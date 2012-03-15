#ifdef HAVE_ML_SOLVER
#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>

//FIXME: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)

#include "EllipticDomain.h"

EllipticDomain::EllipticDomain(Vector_t nr, Vector_t hr) {
    setNr(nr);
    setHr(hr);
}

EllipticDomain::EllipticDomain(double semimajor, double semiminor, Vector_t nr, Vector_t hr, std::string interpl) {
    SemiMajor = semimajor;
    SemiMinor = semiminor;
    setNr(nr);
    setHr(hr);

    if(interpl == "constant")
        interpolationMethod = CONSTANT;
    else if(interpl == "linear")
        interpolationMethod = LINEAR;
    else if(interpl == "quadratic")
        interpolationMethod = QUADRATIC;
}

EllipticDomain::~EllipticDomain() {
    //nothing so far
}


// for this geometry we only have to calculate the intersection on
// one x-y-plane
// for the moment we center the ellipse around the center of the grid
// hr holds the grid-spacings (boundary ellipse embedded in hr-grid)
void EllipticDomain::Compute(Vector_t hr) {

    //there is nothing to be done if the mesh spacings have not changed
    if(hr[0] == getHr()[0] && hr[1] == getHr()[1] && hr[2] == getHr()[2]) {
        hasGeometryChanged_m = false;
        return;
    }

    setHr(hr);
    hasGeometryChanged_m = true;
    //reset number of points inside domain
    nxy_m = 0;


    // clear previous coordinate maps
    IdxMap.clear();
    CoordMap.clear();
    //clear previous intersection points
    IntersectYDir.clear();
    IntersectXDir.clear();

    // build a index and coordinate map
    register int idx = 0;
    register int x, y;
    for(x = 0; x < nr[0]; x++) {
        for(y = 0; y < nr[1]; y++) {

            if(isInside(x, y, 1)) {
                //IdxMap[toCoordIdx(x, y)] = idx++;
                IdxMap[toCoordIdx(x, y)] = idx;
                CoordMap[idx++] = toCoordIdx(x, y);
                nxy_m++;
            }

        }
    }

    switch(interpolationMethod) {

        case CONSTANT:
            break;
        case LINEAR:
        case QUADRATIC:

            double smajsq = SemiMajor * SemiMajor;
            double sminsq = SemiMinor * SemiMinor;
            double yd = 0.0;
            double xd = 0.0;
            double pos = 0.0;
            double mx = (nr[0] - 1) * hr[0] / 2.0;
            double my = (nr[1] - 1) * hr[1] / 2.0;

            //calculate intersection with the ellipse
            for(x = 0; x < nr[0]; x++) {
                pos = x * hr[0] - mx;
                yd = std::abs(sqrt(sminsq - sminsq * pos * pos / smajsq)); // + 0.5*nr[1]*hr[1]);
                IntersectYDir.insert(pair<int, double>(x, yd));
                IntersectYDir.insert(pair<int, double>(x, -yd));
            }

            for(y = 0; y < nr[1]; y++) {
                pos = y * hr[1] - my;
                xd = std::abs(sqrt(smajsq - smajsq * pos * pos / sminsq)); // + 0.5*nr[0]*hr[0]);
                IntersectXDir.insert(pair<int, double>(y, xd));
                IntersectXDir.insert(pair<int, double>(y, -xd));
            }
    }
}

void EllipticDomain::getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod) {
        case CONSTANT:
            ConstantInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
        case LINEAR:
            LinearInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
        case QUADRATIC:
            QuadraticInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void EllipticDomain::getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {


    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);
}


void EllipticDomain::getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B) {

    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

void EllipticDomain::getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B) {

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

void EllipticDomain::ConstantInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;

    W = -1 / hr[0] * 1 / hr[0];
    E = -1 / hr[0] * 1 / hr[0];
    N = -1 / hr[1] * 1 / hr[1];
    S = -1 / hr[1] * 1 / hr[1];
    F = -1 / hr[2] * 1 / hr[2];
    B = -1 / hr[2] * 1 / hr[2];
    C = 2 / hr[0] * 1 / hr[0] + 2 / hr[1] * 1 / hr[1] + 2 / hr[2] * 1 / hr[2];

    // we are a right boundary point
    if(!isInside(x + 1, y, z))
        E = 0.0;

    // we are a left boundary point
    if(!isInside(x - 1, y, z))
        W = 0.0;

    // we are a upper boundary point
    if(!isInside(x, y + 1, z))
        N = 0.0;

    // we are a lower boundary point
    if(!isInside(x, y - 1, z))
        S = 0.0;

    if(z == 1 || z == nr[2] - 2) {

        // case where we are on the Robin BC in Z-direction
        // where we distinguish two cases
        // IFF: this values should not matter because they
        // never make it into the discretization matrix
        if(z == 1)
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

void EllipticDomain::LinearInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;

    double cx = x * hr[0] - (nr[0] - 1) * hr[0] / 2.0;
    double cy = y * hr[1] - (nr[1] - 1) * hr[1] / 2.0;

    double dx = 0.0;
    std::multimap<int, double>::iterator it = IntersectXDir.find(y);
    if(cx < 0)
        it++;
    dx = it->second;

    double dy = 0.0;
    it = IntersectYDir.find(x);
    if(cy < 0)
        it++;
    dy = it->second;


    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    C = 0.0;

    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        C += 1 / ((dx - cx) * de);
        E = 0.0;
    } else {
        C += 1 / (de * de);
        E = -1 / (de * de);
    }

    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        C += 1 / ((std::abs(dx) - std::abs(cx)) * dw);
        W = 0.0;
    } else {
        C += 1 / (dw * dw);
        W = -1 / (dw * dw);
    }

    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        C += 1 / ((dy - cy) * dn);
        N = 0.0;
    } else {
        C += 1 / (dn * dn);
        N = -1 / (dn * dn);
    }

    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        C += 1 / ((std::abs(dy) - std::abs(cy)) * ds);
        S = 0.0;
    } else {
        C += 1 / (ds * ds);
        S = -1 / (ds * ds);
    }

    F = -1 / (hr[2] * hr[2]);
    B = -1 / (hr[2] * hr[2]);
    C += 2 / (hr[2] * hr[2]);

    // handle boundary condition in z direction
    if(z == 0 || z == nr[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        //hr[2]*(nr2[2]-1)/2 = radius
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);

        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;

    }

}

void EllipticDomain::QuadraticInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    double cx = (x - (nr[0] - 1) / 2.0) * hr[0];
    double cy = (y - (nr[1] - 1) / 2.0) * hr[1];

    // since every vector for elliptic domains has ALWAYS size 2 we
    // can catch all cases manually
    double dx = 0.0;
    std::multimap<int, double>::iterator it = IntersectXDir.find(y);
    if(cx < 0)
        it++;
    dx = it->second;

    double dy = 0.0;
    it = IntersectYDir.find(x);
    if(cy < 0)
        it++;
    dy = it->second;

    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    W = 1.0;
    E = 1.0;
    N = 1.0;
    S = 1.0;
    F = 1.0;
    B = 1.0;
    C = 0.0;

    //TODO: = cx+hr[0] > dx && cx > 0
    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    ////we are a right boundary point
    ////if(!isInside(x+1,y,z)) {
    //de = dx-cx;
    //}

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    ////we are a left boundary point
    ////if(!isInside(x-1,y,z)) {
    //dw = std::abs(dx)-std::abs(cx);
    //}

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    ////we are a upper boundary point
    ////if(!isInside(x,y+1,z)) {
    //dn = dy-cy;
    //}

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    ////we are a lower boundary point
    ////if(!isInside(x,y-1,z)) {
    //ds = std::abs(dy)-std::abs(cy);
    //}

    //TODO: = cx+hr[0] > dx && cx > 0
    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        de = dx - cx;
        E = 0.0;
    }

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        dw = std::abs(dx) - std::abs(cx);
        W = 0.0;
    }

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        dn = dy - cy;
        N = 0.0;
    }

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        ds = std::abs(dy) - std::abs(cy);
        S = 0.0;
    }

    //2/dw*(dw_de)
    W *= -1.0 / (dw * (dw + de));
    E *= -1.0 / (de * (dw + de));
    N *= -1.0 / (dn * (dn + ds));
    S *= -1.0 / (ds * (dn + ds));
    F = -1 / (hr[2] * (hr[2] + hr[2]));
    B = -1 / (hr[2] * (hr[2] + hr[2]));

    //TODO: problem when de,dw,dn,ds == 0
    //is NOT a regular BOUND PT
    C += 1 / de * 1 / dw;
    C += 1 / dn * 1 / ds;
    C += 1 / hr[2] * 1 / hr[2];
    scaleFactor = 0.5;


    //for regular gridpoints no problem with symmetry, just boundary
    //z direction is right
    //implement isLastInside(dir)
    //we have LOCAL x,y coordinates!

    /*
       if(dw != 0 && !wIsB)
       W = -1/dw * (dn+ds) * 2*hr[2];
       else
       W = 0;
       if(de != 0 && !eIsB)
       E = -1/de * (dn+ds) * 2*hr[2];
       else
       E = 0;
       if(dn != 0 && !nIsB)
       N = -1/dn * (dw+de) * 2*hr[2];
       else
       N = 0;
       if(ds != 0 && !sIsB)
       S = -1/ds * (dw+de) * 2*hr[2];
       else
       S = 0;
       F = -(dw+de)*(dn+ds)/hr[2];
       B = -(dw+de)*(dn+ds)/hr[2];
       */

    //if(dw != 0)
    //W = -2*hr[2]*(dn+ds)/dw;
    //else
    //W = 0;
    //if(de != 0)
    //E = -2*hr[2]*(dn+ds)/de;
    //else
    //E = 0;
    //if(dn != 0)
    //N = -2*hr[2]*(dw+de)/dn;
    //else
    //N = 0;
    //if(ds != 0)
    //S = -2*hr[2]*(dw+de)/ds;
    //else
    //S = 0;
    //F = -(dw+de)*(dn+ds)/hr[2];
    //B = -(dw+de)*(dn+ds)/hr[2];

    //// RHS scaleFactor for current 3D index
    //// Factor 0.5 results from the SW/quadratic extrapolation
    //scaleFactor = 0.5*(dw+de)*(dn+ds)*(2*hr[2]);

    // catch the case where a point lies on the boundary
    //FIXME: do this more elegant!
    //double m1 = dw*de;
    //double m2 = dn*ds;
    //if(de == 0)
    //m1 = dw;
    //if(dw == 0)
    //m1 = de;
    //if(dn == 0)
    //m2 = ds;
    //if(ds == 0)
    //m2 = dn;
    ////XXX: dn+ds || dw+de can be 0
    ////C = 2*(dn+ds)*(dw+de)/hr[2];
    //C = 2/hr[2];
    //if(dw != 0 || de != 0)
    //C *= (dw+de);
    //if(dn != 0 || ds != 0)
    //C *= (dn+ds);
    //if(dw != 0 || de != 0)
    //C += (2*hr[2])*(dn+ds)*(dw+de)/m1;
    //if(dn != 0 || ds != 0)
    //C += (2*hr[2])*(dw+de)*(dn+ds)/m2;

    //handle Neumann case
    //if(z == 0 || z == nr[2]-1) {

    //if(z == 0)
    //F = 0.0;
    //else
    //B = 0.0;

    ////neumann stuff
    //W = W/2.0;
    //E = E/2.0;
    //N = N/2.0;
    //S = S/2.0;
    //C /= 2.0;

    //scaleFactor /= 2.0;
    //}

    // handle boundary condition in z direction
    if(z == 0 || z == nr[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        //C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);
        //hr[2]*(nr2[2]-1)/2 = radius
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);

        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor /= 2.0;

    }
}


#endif //#ifdef HAVE_ML_SOLVER
