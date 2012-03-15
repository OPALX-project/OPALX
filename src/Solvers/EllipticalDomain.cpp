#ifdef HAVE_ML_SOLVER

//TODO: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)

#include "EllipticalDomain.h"

EllipticalDomain::EllipticalDomain(Vector_t nr, Vector_t hr)
{
    setNr(nr);
    setHr(hr);
}

EllipticalDomain::EllipticalDomain(double semimajor, double semiminor, Vector_t nr, Vector_t hr, std::string interpl)
{
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

//for this boundary we only have to calculate the intersection on
//one x-y-plane
//nr = new gridsize with ellipse around
void EllipticalDomain::Compute(Vector_t hr)
{

    setHr(hr);
    nxy_m = 0;

    register int x,y;
    register double pos,smajsq,sminsq;

    //clear previous coordinate maps
    IdxMap.clear();
    CoordMap.clear();
    register int idx = 0;

    //XXX: SINCE WE COUNTING idx++ LOOP ORDER MATTERS!
    for(x=0; x < nr[0]; x++) {
        for(y=0; y < nr[1]; y++) {

            if(isInside(x,y,0)) {
            //if(isInside(x+1,y,0) && isInside(x-1,y,0) && isInside(x,y+1,0) && isInside(x,y-1,0)) {
                IdxMap[toCoordIdx(x,y)] = idx++;
                nxy_m++;
            }

        }
    }
    idx = 0;
    for(x=0; x < nr[0]; x++) {
        for(y=0; y < nr[1]; y++) {
            if(isInside(x,y,0)) 
            //if(isInside(x+1,y,0) && isInside(x-1,y,0) && isInside(x,y+1,0) && isInside(x,y-1,0)) 
                CoordMap[idx++] = toCoordIdx(x,y);
        }
    }

    switch(interpolationMethod) {
    
        case CONSTANT: break;
        case LINEAR:
        case QUADRATIC:

            smajsq = SemiMajor*SemiMajor;
            sminsq = SemiMinor*SemiMinor;

            //clear previous intersection points
            IntersectYDir.clear();
            IntersectXDir.clear();
            //IntersectYDir.count(2) == 2!
            double yd,xd;

            for(x=0; x<nr[0]; x++) {
                pos = (x-floor(nr[0]/2))*hr[0]; 
                yd = abs(sqrt(sminsq - sminsq*pos*pos/smajsq));
                IntersectYDir.insert(pair<int, double>(x, yd));
                IntersectYDir.insert(pair<int, double>(x, -yd));
            }

            for(y=0; y<nr[1]; y++) {
                pos = (y-floor(nr[1]/2))*hr[1]; 
                xd = abs(sqrt(smajsq - smajsq*pos*pos/sminsq));
                IntersectXDir.insert(pair<int, double>(y, xd));
                IntersectXDir.insert(pair<int, double>(y, -xd));
            }

    }

}

int EllipticalDomain::getNumXY(int z) {

    return nxy_m;

}

// ignore z since its the same for every z
std::vector<double> EllipticalDomain::getYDirIntersect(int x, int z) {
    std::vector<double> ret;
    for(std::multimap<int, double>::iterator it = IntersectYDir.find(x); it != IntersectYDir.end(); it++)
        ret.push_back(it->second);
    return ret;
}

void EllipticalDomain::getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor) {

    //determine which interpolation method we use for points near the boundary
    switch(interpolationMethod) {

        case CONSTANT: 
            ConstantInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); break;

        case LINEAR: 
            LinearInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); 
            break;

        case QUADRATIC: 
            QuadraticInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); 
            break;
    }

    //simple check if center value of stencil is positive
#ifdef DEBUG
    if(C <= 0)
        cout << "Stencil C is <= 0! This should not case should never occure!" << endl;
#endif
}

//TODO: remove isInside()
//TODO: if constant interpolation is used we dont need to calculate intersections
void EllipticalDomain::ConstantInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    //if(!isInside(x+1,y,z) || !isInside(x-1,y,z) || !isInside(x,y+1,z) || !isInside(x,y-1,z)) {
    //    cout << x << "/" << y << endl;
    //    E=0; W=0; S=0; N=0; F=0; B=0; return;
    //}

    //scaleFactor should be ok. the mesh sizes are _in_ the stencil
    scaleFactor = 1.0;

    W = -1/hr[0]*1/hr[0];
    E = -1/hr[0]*1/hr[0];
    N = -1/hr[1]*1/hr[1];
    S = -1/hr[1]*1/hr[1];
    F = -1/hr[2]*1/hr[2];
    B = -1/hr[2]*1/hr[2];
    C = 2/hr[0]*1/hr[0] + 2/hr[1]*1/hr[1] + 2/hr[2]*1/hr[2];

    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0 /*&& x*hr[0] < dx*/) {
    if(!isInside(x+1,y,z))
        E = 0.0; //we are a right boundary point

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0 /*&& x*hr[0] > dx*/) {
    if(!isInside(x-1,y,z)) 
        W = 0.0; //we are a left boundary point

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0 /*&& y*hr[1] < dy*/) {
    if(!isInside(x,y+1,z)) 
        N = 0.0; //we are a upper boundary point

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0 /*&& y*hr[1] > dy*/) {
    if(!isInside(x,y-1,z)) 
        S = 0.0; //we are a lower boundary point

    if(false /*z == 0 || z == nr[2]-1*/) {

        //case where we are on the Robin BC in Z-direction
        //where we distinguish two cases  
        //IFF: this values should not matter because they 
        //never make it into the discretization matrix
        if(z == 0) 
            F = 0.0;
        else 
            B = 0.0;
       
        //add contribution of Robin discretization to center point
        //d the distance between the center of the bunch and the boundary
        //double cx = (x-(nr[0]-1)/2)*hr[0];
        //double cy = (y-(nr[1]-1)/2)*hr[1];
        //double cz = hr[2]*(nr[2]-1);
        //double d = sqrt(cx*cx+cy*cy+cz*cz);
        double d = hr[2]*(nr[2]-1)/2;
        //cout << cx << " " << cy << " " << cz << " " << 2/(d*hr[2]) << endl;
        C += 2/(d * hr[2]);
        //C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);

        //scale all stencil-points in z-plane with 0.5 (Neumann discretization)
        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;
    }

}

//TODO: remove isInside()
//TODO: scaleFactor -- check
void EllipticalDomain::LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    scaleFactor = 1.0;

    double cx = (x-floor((double)(nr[0]/2.0)))*hr[0];
    double cy = (y-floor((double)(nr[1]/2.0)))*hr[1];

    //since every vector for elliptic domains has ALWAYS size 2 we 
    //can catch all cases manually
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


    double dw=hr[0];
    double de=hr[0];
    double dn=hr[1];
    double ds=hr[1];
    F = -1/hr[2]*1/hr[2];
    B = -1/hr[2]*1/hr[2];
    C = 2*1/hr[2]*1/hr[2];

    //TODO: cx > 0 --> xc <= dx

    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0 /*&& x*hr[0] < dx*/) {
    if(!isInside(x+1,y,z)) {
        //we are a right boundary point
        C += -1/(de*de)*(1-de/(dx-cx)); //1/(dx-cx)*1/(de*de);
        de = 0.0;
        E = 0.0;
    } else {
        C += 1/(de*de);
        E = -1/de*1/de;
    }

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0 /*&& x*hr[0] > dx*/) {
    if(!isInside(x-1,y,z)) {
        //we are a left boundary point
        C += -1/(dw*dw)*(1-dw/(abs(dx)-abs(cx))); //;1/(abs(dx)-abs(cx))*1/(dw*dw);
        dw = 0.0;
        W = 0.0;
    } else {
        C += 1/(dw*dw);
        W = -1/dw*1/dw;
    }

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0 /*&& y*hr[1] < dy*/) {
    if(!isInside(x,y+1,z)) {
        //we are a upper boundary point
        C += -1/(dn*dn)*(1-dn/(dy-cy)); //1/(dy-cy)*1/(dn*dn);
        dn = 0.0;
        N = 0.0;
    } else {
        C += 1/(dn*dn);
        N = -1/dn*1/dn;
    }

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0 /*&& y*hr[1] > dy*/) {
    if(!isInside(x,y-1,z)) {
        //we are a lower boundary point
        C += -1/(ds*ds)*(1-ds/(abs(dy)-abs(cy))); //1/(abs(dy)-abs(cy))*1/(ds*ds);
        ds = 0.0;
        S = 0.0;
    } else {
        C += 1/(ds*ds);
        S = -1/ds*1/ds;
    }

    if(z == 0 || z == nr[2]-1) {

        //case where we are on the NEUMAN BC in Z-direction
        //where we distinguish two cases  

        if(z == 0)
            F = 0.0;
        else
            B = 0.0;
        
        //hr[2]*(nr2[2]-1)/2 = radius
        C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);
        //C += 2*hr[2]/(hr[2]*nr[2]/2.0);

        //neumann stuff
        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;

    }

}

void EllipticalDomain::QuadraticInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    double cx = (x-floor((double)(nr[0]/2.0)))*hr[0];
    double cy = (y-floor((double)(nr[1]/2.0)))*hr[1];

    //since every vector for elliptic domains has ALWAYS size 2 we 
    //can catch all cases manually
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

    double dw=hr[0];
    double de=hr[0];
    double dn=hr[1];
    double ds=hr[1];
    C = 0.0;

    //TODO: = cx+hr[0] > dx && cx > 0
    if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
        //we are a right boundary point
        //if(!isInside(x+1,y,z)) {
        de = dx-cx;
    }

    if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
        //we are a left boundary point
        //if(!isInside(x-1,y,z)) {
        dw = abs(dx)-abs(cx);
    }

    if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
        //we are a upper boundary point
        //if(!isInside(x,y+1,z)) {
        dn = dy-cy;
    }

    if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
        //we are a lower boundary point
        //if(!isInside(x,y-1,z)) {
        ds = abs(dy)-abs(cy);
    }

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

    if(dw != 0)
        W = -2*hr[2]*(dn+ds)/dw;
    else
        W = 0;
    if(de != 0)
        E = -2*hr[2]*(dn+ds)/de;
    else
        E = 0;
    if(dn != 0)
        N = -2*hr[2]*(dw+de)/dn;
    else
        N = 0;
    if(ds != 0)
        S = -2*hr[2]*(dw+de)/ds;
    else
        S = 0;
    F = -(dw+de)*(dn+ds)/hr[2];
    B = -(dw+de)*(dn+ds)/hr[2];

    //RHS scaleFactor for current 3D index
    //Factor 0.5 results from the SW/quadratic extrapolation
    scaleFactor = 0.5*(dw+de)*(dn+ds)*(2*hr[2]);

    //catch the case where a point lies on the boundary
    //TODO: do this more elegant!
    double m1 = dw*de;
    double m2 = dn*ds;
    if(de == 0)
        m1 = dw;
    if(dw == 0)
        m1 = de;
    if(dn == 0)
        m2 = ds;
    if(ds == 0)
        m2 = dn;
    //XXX: dn+ds || dw+de can be 0
    //C = 2*(dn+ds)*(dw+de)/hr[2]; 
    C = 2/hr[2]; 
    if(dw != 0 || de != 0)
        C *= (dw+de);
    if(dn != 0 || ds != 0)
        C *= (dn+ds);
    if(dw != 0 || de != 0)
        C += (2*hr[2])*(dn+ds)*(dw+de)/m1;
    if(dn != 0 || ds != 0)
        C += (2*hr[2])*(dw+de)*(dn+ds)/m2;

    //handle Neumann case
    if(z == 0 || z == nr[2]-1) {

        if(z == 0) 
            F = 0.0;
        else 
            B = 0.0;


        //neumann stuff
        W = W/2.0;
        E = E/2.0;
        N = N/2.0;
        S = S/2.0;
        C /= 2.0;

        scaleFactor /= 2.0;

    }
}

void EllipticalDomain::getBoundaryStencil(int idx, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor)
{

    //TODO: reverse map search?
    //double mem or O(n) search in map to get x,y,z from idx

    int x=0,y=0,z=0;

    getCoord(idx,x,y,z);
    getBoundaryStencil(x,y,z,W,E,S,N,F,B,C,scaleFactor);

}

void EllipticalDomain::getNeighbours(int idx, double& W, double& E, double& S, double& N, double& F, double& B) 
{

    int x=0,y=0,z=0;

    getCoord(idx,x,y,z);
    getNeighbours(x,y,z,W,E,S,N,F,B);

}

void EllipticalDomain::getNeighbours(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B)
{

    if(x > 0)
        W = getIdx(x-1,y,z);
    else
        W = -1;
    if(x < nr[0]-1)
        E = getIdx(x+1,y,z);
    else
        E = -1;

    if(y < nr[1]-1)
        N = getIdx(x,y+1,z);
    else
        N = -1;
    if(y > 0)
        S = getIdx(x,y-1,z);
    else
        S = -1;

    if(z > 0)
        F = getIdx(x,y,z-1);
    else
        F = -1;
    if(z < nr[2]-1)
        B = getIdx(x,y,z+1);
    else
        B = -1;
}

#endif //#ifdef HAVE_ML_SOLVER
