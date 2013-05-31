#ifdef HAVE_ML_SOLVER

#include "ArbitraryDomain.h"
#include <iostream>

//IFF: simplified case where we have intersect = 2
//IFF: every node counts own number of gridpoints and then sum to global procs before MyPID (including ghost layers in z direction).
//IFF: triangles intersecting in z dir? (cathode?) -> boundary condition there? neumann? then we should take a normal and stuff!?

ArbitraryDomain::ArbitraryDomain(string fname, Vector_t nr_, Vector_t hr_, NDIndex<3> locidx) {

    filename = fname;
    setNr(nr_);
    setHr(hr_);
    localidx = locidx;
    startIdx = 0;

    LoadFile();

    cout << "file loaded" << endl;

    //TODO: FIT GEOMETRY ON GRID! --> ENLARGE GRID!!

}

//IFF: this method should work correctly since it was checked in render_triangle.cpp
void ArbitraryDomain::LoadFile() {

    //FOR ANALYTICAL TEST:
    double xshift = 0.5;
    double yshift = 0.5;
    /////////////////////

    h5_err_t rc;

    h5_file *f = H5OpenFile(filename.c_str(), 0);
    if(f == NULL)
        ERRORMSG("can't open file: "  << filename << endl);

    h5_size_t num_meshes = H5FedGetNumMeshes(f, TRIANGLE_MESH);
    if(num_meshes != 1)
        ERRORMSG("can't handle more or less than one mesh!!" << endl);

    h5_id_t mesh_id = 0;
    rc = H5FedOpenMesh(f, mesh_id, TRIANGLE_MESH);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    h5_id_t level_id;
    h5_size_t num_levels = H5FedGetNumLevels(f);
    //if(num_levels != 1)
    //  *gmsg << "cannot handle more refinement levels!!" << endl;

    for(level_id = 0; level_id < num_levels; level_id++) {

        rc = H5FedSetLevel(f, level_id);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        h5_id_t id, local_id, parent_id;
        h5_size_t real_num = 0;
        vertex_t tmp;

        h5_id_t level_id = H5FedGetLevel(f);
        h5_size_t num = H5FedGetNumVertices(f);
        rc = H5FedStartTraverseVertices(f);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        //forall vertices
        while((real_num < num) && ((local_id = H5FedTraverseVertices(f, &id, tmp.P)) >= 0)) {
            //FOR ANALYTICAL TEST:
            tmp.P[0] += xshift;
            tmp.P[1] += yshift;
            /////////////////////
            vertices[id] = tmp;
            real_num++;
        }

        num = H5FedGetNumTriangles(f);
        rc = H5FedStartTraverseTriangles(f);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        entity_t triangle;
        h5_id_t vids[3];
        real_num = 0;

        //for all triangles
        while((real_num < num) && ((local_id = H5FedTraverseTriangles(f, &id, &parent_id, vids)) >= 0)) {
            triangle.global_id = id;
            triangle.parent_id = parent_id;
            triangle.vertexIDs[0] = vids[0];
            triangle.vertexIDs[1] = vids[1];
            triangle.vertexIDs[2] = vids[2];
            entities[id] = triangle;
            real_num++;
        }

    }

}

void ArbitraryDomain::Compute(Vector_t hr) {

    setHr(hr);

    //h5_float64_t edge1, edge2, t, q, p, P0;
    double edge1[3], edge2[3], t[3], q[3], p[3], P0[3];
    double det, invDet, u, v, tt;
    double dir[3], origin[3];
    int zGhostOffsetLeft = (localidx[2].first() == 0) ? 0 : 1;
    int zGhostOffsetRight = (localidx[2].last() == nr[2] - 1) ? 0 : 1;
    multimap< pair<int, int>, double >::iterator it;
    pair< multimap< pair<int, int>, double>::iterator, multimap< pair<int, int>, double>::iterator > ret;
    bool hit;

    std::map<h5_id_t, entity_t>::iterator triangleitr;
    //forall triangles
    for(triangleitr = entities.begin(); triangleitr != entities.end(); triangleitr++) {


        //P0 = vertices[triangleitr->second.vertexIDs[0]];
        //edge1 = vertices[triangleitr->second.vertexIDs[1]];
        //edge2 = vertices[triangleitr->second.vertexIDs[2]];

        vertex_t tmp;
        tmp = vertices[triangleitr->second.vertexIDs[0]];
        P0[0] = tmp.P[0];
        P0[1] = tmp.P[1];
        P0[2] = tmp.P[2];

        tmp = vertices[triangleitr->second.vertexIDs[1]];
        edge1[0] = tmp.P[0] - P0[0];
        edge1[1] = tmp.P[1] - P0[1];
        edge1[2] = tmp.P[2] - P0[2];

        tmp = vertices[triangleitr->second.vertexIDs[2]];
        edge2[0] = tmp.P[0] - P0[0];
        edge2[1] = tmp.P[1] - P0[1];
        edge2[2] = tmp.P[2] - P0[2];

        //IFF: use normed dir -> t = intersection (?)
        //x dir rays = forall points in origin y,z with dir (1,0,0) (Y,Z -> COORD)
        dir[0] = 1;
        dir[1] = 0;
        dir[2] = 0;
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        for(int y = 0; y < nr[1]; y++) {
            //ONE LAYER MORE = GHOST LAYER
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].last() + zGhostOffsetRight; z++) {

                origin[2] = z * hr[2];

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
            origin[1] += hr[1];
        }

        //y dir rays = forall points in origin x,z with dir (0,1,0) (X,Z -> COORD)
        dir[0] = 0;
        dir[1] = 1;
        dir[2] = 0;
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        for(int x = 0; x < nr[0]; x++) {
            //ONE LAYER MORE = GHOST LAYER
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].first() + zGhostOffsetRight; z++) {

                origin[2] = z * hr[2];

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
            origin[0] += hr[0];
        }

        //z dir rays = forall points in origin x,y with dir (0,0,1) (X,Y -> COORD)
        dir[0] = 0;
        dir[1] = 0;
        dir[2] = 1;
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        for(int x = 0; x < nr[0]; x++) {
            for(int y = 0; y < nr[1]; y++) {

                origin[1] = y * hr[1];

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
            origin[0] += hr[0];
        }

    }

    cout << "number of intersections in x-dir: " << IntersectXDir.size() << endl;
    cout << "number of intersections in y-dir: " << IntersectYDir.size() << endl;
    cout << "number of intersections in z-dir: " << IntersectZDir.size() << endl;

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

    cout << "ghost nodes left: " << numGhostNodesLeft << endl;

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

    cout << "number of gridpoints: " << numtotal << endl;

    startIdx = 0;
    MPI_Scan(&numtotal, &startIdx, 1, MPI_INTEGER, MPI_SUM, Ippl::getComm());
    startIdx -= numtotal;

    cout << "start idx: " << startIdx << endl;

    //build up index and coord map
    IdxMap.clear();
    CoordMap.clear();
    register int idx = startIdx - numGhostNodesLeft;

    for(int x = 0; x < nr[0]; x++) {
        for(int y = 0; y < nr[1]; y++) {
            for(int z = localidx[2].first() - zGhostOffsetLeft; z <= localidx[2].last() + zGhostOffsetRight; z++) {

                if(isInside(x, y, z)) {
                    IdxMap[toCoordIdx(x, y, z)] = idx;
                    CoordMap[idx] = toCoordIdx(x, y, z);
                    idx++;
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
    double cx = x * hr[0];
    double cy = y * hr[1];
    double cz = z * hr[2];
    double val1 = 0;
    double val2 = 0;
    multimap < std::pair<int, int>, double >::iterator itr;

    //check if x is inside with y,z coords
    std::pair<int, int> coordyz(y, z);
    itr = IntersectXDir.find(coordyz);
    int count = IntersectXDir.count(coordyz);
    if(count == 1)
        ret = (cx == itr->second);
    else {
        val1 = itr->second;
        val2 = (++itr)->second;
        ret = (cx <= max(val1, val2)) && (cx >= min(val1, val2));
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
        ret = ret && (cy <= max(val1, val2)) && (cy >= min(val1, val2));
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
      ret = ret && (cz <= max(val1, val2)) && (cz >= min(val1, val2));
    }*/

    //ret = ret && (cz >= 0 && cz <= nr[2]);

    return ret;

}

int ArbitraryDomain::getNumXY(int z) {

    return numXY[z];

}

void ArbitraryDomain::getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    double cx = x * hr[0];
    double cy = y * hr[1];
    double cz = z * hr[2];

    multimap< pair<int, int>, double >::iterator it;
    pair< multimap< pair<int, int>, double>::iterator, multimap< pair<int, int>, double>::iterator > ret;

    //since every vector here is only an array[2] we
    //can catch all cases manually
    double dx = -1.0, dy = -1.0, dz = -1.0;
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
      if(abs(it->second - z) < hr[2]) {
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
    if(!isInside(x,y,z+1))
      df = dz;

    //we are a lower boundary point
    if(!isInside(x,y,z-1))
      db = dz;
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
    if(dw + de != 0)
        scaleFactor *= (dw + de);
    if(dn + ds != 0)
        scaleFactor *= (dn + ds);
    if(df + db != 0)
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

    //if(C <= 0)
    //  cout << "!!!!!!!!!! C is <= 0: " << C <<  endl;

    //handle Neumann case
    if(z == 0 || z == nr[2] - 1) {

        if(z == 0) {
            F = 0.0;
            B = -hr[2];
            C = 2 * hr[2];
        } else {
            B = 0.0;
            F = -hr[2];
            C = 2 * hr[2];
        }

        //neumann stuff
        W = W / 2.0;
        E = E / 2.0;
        N = N / 2.0;
        S = S / 2.0;
        C /= 2.0;

        scaleFactor /= 2.0;

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

    //IFF: reverse map search?
    //double mem or O(n) search in map to get x,y,z from idx

    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);

}

void ArbitraryDomain::getNeighbours(int idx, double &W, double &E, double &S, double &N, double &F, double &B) {

    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

//IFF:should be ok
//getIdx returns -1 when not inside domain
void ArbitraryDomain::getNeighbours(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B) {

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

#endif //#ifdef HAVE_ML_SOLVER
