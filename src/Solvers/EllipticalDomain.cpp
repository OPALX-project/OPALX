#ifdef HAVE_ML_SOLVER

//TODO: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)

#include "EllipticalDomain.h"

EllipticalDomain::EllipticalDomain(Vector_t nr, Vector_t hr)
{
  setNr(nr);
  setHr(hr);
}

EllipticalDomain::EllipticalDomain(double semimajor, double semiminor, Vector_t nr, Vector_t hr)
{
  SemiMajor = semimajor;
  SemiMinor = semiminor;
  setNr(nr);
  setHr(hr);
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

  //build up index map
  IdxMap.clear();
  CoordMap.clear();
  register int idx = 0;
  
  //XXX: SINCE WE COUNTING idx++ LOOP ORDER MATTERS!
  for(x=0; x < nr[0]; x++) {
    for(y=0; y < nr[1]; y++) {

        if(isInside(x,y,0)) {
          IdxMap[toCoordIdx(x,y)] = idx++;
          nxy_m++;
        }
        
    }
  }
  idx = 0;
  for(x=0; x < nr[0]; x++) {
    for(y=0; y < nr[1]; y++) {
        if(isInside(x,y,0)) 
          CoordMap[idx++] = toCoordIdx(x,y);
    }
  }

  smajsq = SemiMajor*SemiMajor;
  sminsq = SemiMinor*SemiMinor;

  //clear map
  IntersectYDir.clear();
  IntersectXDir.clear();
  //IntersectYDir.count(2) == 2!
  double yd,xd;

  for(x=0; x<nr[0]; x++) {
    pos = (x-floor(nr[0]/2))*hr[0]; 
    yd = abs(sqrt(sminsq - sminsq*pos*pos/smajsq)); // + 0.5*nr[1]*hr[1]);
    IntersectYDir.insert(pair<int, double>(x, yd));
    IntersectYDir.insert(pair<int, double>(x, -yd));
  }

  for(y=0; y<nr[1]; y++) {
    pos = (y-floor(nr[1]/2))*hr[1]; 
    xd = abs(sqrt(smajsq - smajsq*pos*pos/sminsq)); // + 0.5*nr[0]*hr[0]);
    IntersectXDir.insert(pair<int, double>(y, xd));
    IntersectXDir.insert(pair<int, double>(y, -xd));
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

void EllipticalDomain::getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {


  double cx = (x-floor(nr[0]/2.0))*hr[0];
  double cy = (y-floor(nr[1]/2.0))*hr[1];

  //since every vector here is only an array[2] we 
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

  //symmetric: scale RHS too
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
  //0.5* comes from discretiztaion
  //scaleFactor = 0.5*(dw+de)*(dn+ds)*(df+db);
  scaleFactor = 0.5*(dw+de)*(dn+ds)*(2*hr[2]);

  //catch the case where a point lies on the boundary
  //IFF: do this more elegant!
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
  //IFF: dn+ds || dw+de can be 0
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

  //if(C <= 0)
  //  cout << "!!!!!!!!!! C is <= 0: " << C <<  endl;

  //handle Neumann case
  if(z == 0 || z == nr[2]-1) {
  
    if(z == 0) {
      F = 0.0;
      B = -hr[2];
      C = 2*hr[2];
    } else {
      B = 0.0;
      F = -hr[2];
      C = 2*hr[2];
    }
    
    //neumann stuff
    W = W/2.0;
    E = E/2.0;
    N = N/2.0;
    S = S/2.0;
    C /= 2.0;

    scaleFactor /= 2.0;

    //if(C <= 0)
    //  cout << "!!!!!!!!!! C is <= 0" << endl;
  }

}

/*
 * OLD STENCIL O(1)
void EllipticalDomain::getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C) {


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
  if(idx==2 || idx==3)
    cout << "i: " << idx << " " << x << "/" << y << "/" << z << endl;
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
 
  //IFF: 17.6 swapped signs
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
