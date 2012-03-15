#ifndef TAPER_DOMAIN_H
#define TAPER_DOMAIN_H
#ifdef HAVE_ML_SOLVER

#include <vector>
#include "IrregularDomain.h"

/* Tomas's request for a Taper simulation to compare space charge effects
 * this is a very specific simulation!!
 */

/// EllipticPointList maps from an (x) resp. (y) to a list of double values (=intersections with boundary)
// int encodes point and double intersection value
// FIXME: replace MultiMap with Vector! (performance)
typedef std::multimap<int, double> TaperPointList;

class TaperDomain : public IrregularDomain {

public:

    /// constructor
    TaperDomain(Vector_t nr, Vector_t hr);
    /// constructor
    TaperDomain(double rb, double rs, Vector_t nr, Vector_t hr, std::string interpl, double zmax_big);
    ~TaperDomain();

    /// calculates intersection with the elliptic beam pipe
    void Compute(Vector_t hr);
    //TODO: do we need to export this function??
    std::vector<double> getYDirIntersect(int x, int z);
    /// returns number of nodes in xy plane (here independent of z coordinate)
    int getNumXY(int z);
    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns discretization at 3D index
    void getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B);
    /// returns index of neighbours at 3D index
    void getNeighbours(int idx, double &W, double &E, double &S, double &N, double &F, double &B);
    /// returns type of boundary condition
    string getType() {return "Taper";}
    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int z) {

        double xx = (x - (nr[0] - 1) / 2.0) * hr[0];
        double yy = (y - (nr[1] - 1) / 2.0) * hr[1];
        double zz = z * hr[2] + z_min;

        if(zz <= z_max_big)
            return (xx * xx + yy * yy < radius_big * radius_big);
        else
            return (xx * xx + yy * yy < radius_small * radius_small);
    }

    /// set semi-minor
    void setRadiusBig(double sm) {radius_big = sm;}
    /// set semi-major
    void setRadiusSmall(double sm) {radius_small = sm;}
    /// setters for bunch zmin and zmax
    void setZmin(double zm) {z_min = zm;}
    void setZmax(double zm) {z_max = zm;}

    double getXRangeMin() { return -radius_big; }
    double getXRangeMax() { return radius_big;  }
    double getYRangeMin() { return -radius_big; }
    double getYRangeMax() { return radius_big;  }

    //TODO: ?
    int getStartIdx() {return 0;}

private:

    /// all intersection points with gridlines in Y direction
    TaperPointList IntersectYDir;
    TaperPointList Intersectydir;
    /// all intersection points with gridlines in X direction
    TaperPointList IntersectXDir;
    TaperPointList Intersectxdir;
    /// mapping (x,y,z) -> idx
    std::map<int, int> IdxMap;
    /// mapping idx -> (x,y,z)
    std::map<int, int> CoordMap;

    /// semi-major of the ellipse
    double radius_big;
    /// semi-minor of the ellipse
    double radius_small;
    double z_max_big;
    double z_min;
    double z_max;

    /// number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;
    /// interpolation type
    int interpolationMethod;

    inline int toCoordIdx(int x, int y, int z) { return (z * nr[1] + y) * nr[0] + x; }
    int getIdx(int x, int y, int z) {
        if(isInside(x, y, z) && x >= 0 && y >= 0 && z >= 0)
            return IdxMap[toCoordIdx(x, y, z)];
        else
            return -1;
    }
    inline void getCoord(int idx, int &x, int &y, int &z) {

        int idxx = CoordMap[idx];

        x = idxx % (int)nr[0];
        idxx /= nr[0];
        y = idxx % (int)nr[1];
        idxx /= nr[1];
        z = idxx;

    }
    /*
    /// conversion from (x,y) to index in xy plane
    inline int toCoordIdx(int x, int y) { return y*nr[0] + x; }
    /// conversion from (x,y,z) to index on the 3D grid
    inline int getIdx(int x, int y, int z) {
        if(isInside(x,y,z) && x >= 0 && y >= 0 && z >= 0)
            return IdxMap[toCoordIdx(x,y)]+z*nxy_m;
        else
            return -1;
    }
    /// conversion from a 3D index to (x,y,z)
    inline void getCoord(int idx, int& x, int& y, int& z) {
        int ixy = idx % nxy_m;
        int xy = CoordMap[ixy];
        int inr = (int)nr[0];
        x = xy % inr;
        y = (xy-x)/nr[0];
        z = (idx-ixy)/nxy_m;
    }
    */

    /// different interpolation methods for boundary points
    void ConstantInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void LinearInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void QuadraticInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);

};

#endif //#ifdef HAVE_ML_SOLVER
#endif //#ifdef ELLIPTICAL_DOMAIN_H
