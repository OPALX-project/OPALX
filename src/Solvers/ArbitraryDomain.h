#ifndef ARBITRARY_DOMAIN_H
#define ARBITRARY_DOMAIN_H

#ifdef HAVE_SAAMG_SOLVER

#include <mpi.h>
#include <hdf5.h>
#include "H5hut.h"

#include <vector>
#include <map>
#include <string>
#include <math.h>
#include <cmath>
#include "IrregularDomain.h"
#include "Structure/BoundaryGeometry.h"

class ArbitraryDomain : public IrregularDomain {

public:

    ArbitraryDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr, std::string interpl);
    ArbitraryDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr, Vector_t globalMeanR, Vektor<double, 4> globalToLocalQuaternion, std::string interpl);

    ~ArbitraryDomain();

    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns discretization at 3D index
    void getBoundaryStencil(int idxyz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int idx, int idy, int idz, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns index of neighbours at 3D index
    void getNeighbours(int idxyz, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns type of boundary condition
    std::string getType() {return "Geometric";}
    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int idx, int idy, int idz); 

    /// returns number of nodes in xy plane
    int getNumXY(int idz);
    /// calculates intersection 
    void Compute(Vector_t hr);
    void Compute(Vector_t hr, NDIndex<3> localId);
    // calculates intersection with rotated and shifted geometry 
    void Compute(Vector_t hr, NDIndex<3> localId, Vector_t globalMeanR, Vektor<double, 4> globalToLocalQuaternion);

    int getStartId() {return startId;}

    double getXRangeMin() { return Geo_mincoords_m(0); }
    double getXRangeMax() { return Geo_maxcoords_m(0); }
    double getYRangeMin() { return Geo_mincoords_m(1); }
    double getYRangeMax() { return Geo_maxcoords_m(1); }
    double getZRangeMin() { return Geo_mincoords_m(2); }
    double getZRangeMax() { return Geo_maxcoords_m(2); }

    bool hasGeometryChanged() { return hasGeometryChanged_m; }

private:
    BoundaryGeometry *bgeom_m;

    /// PointList maps from an (x,z) resp. (y,z) pair to double values (=intersections with boundary)
    typedef std::multimap< std::tuple<int, int, int>, double > PointList;

    /// all intersection points with gridlines in X direction
    PointList IntersectHiX, IntersectLoX;

    /// all intersection points with gridlines in Y direction
    PointList IntersectHiY, IntersectLoY;

    /// all intersection points with gridlines in Z direction
    PointList IntersectHiZ, IntersectLoZ;
    
    // meanR to shift from global to local frame 
    Vector_t globalMeanR_m;
    Vektor<double, 4> globalToLocalQuaternion_m;
    Vektor<double, 4> localToGlobalQuaternion_m;

    int startId;

    // Here we store the number of nodes in a xy layer for a given z coordinate
    std::map<int, int> numXY;
    std::map<int, int> numYZ;
    std::map<int, int> numXZ;

    // Number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;
    // mapping (x,y,z) -> idxyz
    std::map<int, int> IdxMap;
    // Mapping idxyz -> (x,y,z)
    std::map<int, int> CoordMap;
    // Interpolation type
    int interpolationMethod;
 
    // Flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    Vector_t Geo_nr_m;
    Vector_t Geo_hr_m;
    Vector_t Geo_mincoords_m;
    Vector_t Geo_maxcoords_m;

    // Conversion from (x,y,z) to index in xyz plane
    inline int toCoordIdx(int idx, int idy, int idz);
    // Conversion from (x,y,z) to index on the 3D grid
    int getIdx(int idx, int idy, int idz);
    // Conversion from a 3D index to (x,y,z)
    inline void getCoord(int idxyz, int &x, int &y, int &z);

    inline void crossProduct(double A[], double B[], double C[]);
    inline double dotProduct(double v1[], double v2[]) { return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]); }

    // Different interpolation methods for boundary points
    void ConstantInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void LinearInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void QuadraticInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);

    // Rotate positive axes with quaternion -DW
    inline void rotateWithQuaternion(Vector_t &v, Vektor<double, 4> const quaternion);
 	
    inline void rotateXAxisWithQuaternion(Vector_t &v, Vektor<double, 4> const quaternion);
    inline void rotateYAxisWithQuaternion(Vector_t &v, Vektor<double, 4> const quaternion);
    inline void rotateZAxisWithQuaternion(Vector_t &v, Vektor<double, 4> const quaternion);
};

#endif //#ifdef HAVE_SAAMG_SOLVER
#endif //#ifdef ARBITRARY_DOMAIN
