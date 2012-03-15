#ifndef ELLIPTICAL_DOMAIN_H
#define ELLIPTICAL_DOMAIN_H
#ifdef HAVE_ML_SOLVER

#include <vector>
#include "IrregularDomain.h"

/// EllipticPointList maps from an (x) resp. (y) to a list of double values (=intersections with boundary) 
// int encodes point and double intersection value
// FIXME: replace MultiMap with Vector! (performance):
//typedef std::vector< vector<double> > EllipticPointList;
typedef std::multimap<int, double> EllipticPointList;

class EllipticalDomain : public IrregularDomain {

    public:

        /// constructor
        EllipticalDomain(Vector_t nr, Vector_t hr);
        /// constructor
        EllipticalDomain(double semimajor, double semiminor, Vector_t nr, Vector_t hr, std::string interpl);
        ~EllipticalDomain();

        /// calculates intersection with the elliptic beam pipe
        void Compute(Vector_t hr);
        //TODO: do we need to export this function??
        std::vector<double> getYDirIntersect(int x, int z);
        /// returns number of nodes in xy plane (here independent of z coordinate)
        int getNumXY(int z);
        /// returns discretization at (x,y,z)
        void getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor);
        /// returns discretization at 3D index
        void getBoundaryStencil(int idx, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor);
        /// returns index of neighbours at (x,y,z)
        void getNeighbours(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B);
        /// returns index of neighbours at 3D index
        void getNeighbours(int idx, double& W, double& E, double& S, double& N, double& F, double& B);
        /// returns type of boundary condition
        string getType() {return "Elliptic";}
        /// queries if a given (x,y,z) coordinate lies inside the domain
        inline bool isInside(int x, int y, int z) {
            double xx = (x-(nr[0]-1)/2.0)*hr[0];
            double yy = (y-(nr[1]-1)/2.0)*hr[1];
            return (xx*xx/(SemiMajor*SemiMajor) + yy*yy/(SemiMinor*SemiMinor) < 1); //CHANGED: <= to <
        }

        /// set semi-minor
        void setSemiMinor(double sm) {SemiMinor = sm;}
        /// set semi-major
        void setSemiMajor(double sm) {SemiMajor = sm;}
  
        double getXRangeMin() { return -SemiMajor; }
        double getXRangeMax() { return SemiMajor;  }
        double getYRangeMin() { return -SemiMinor; }
        double getYRangeMax() { return SemiMinor;  }

        //TODO: ?
        int getStartIdx() {return 0;}

    private:

        /// all intersection points with gridlines in Y direction
        EllipticPointList IntersectYDir;
        /// all intersection points with gridlines in X direction
        EllipticPointList IntersectXDir;
        /// mapping (x,y,z) -> idx
        std::map<int, int> IdxMap;
        /// mapping idx -> (x,y,z)
        std::map<int, int> CoordMap;

        /// semi-major of the ellipse
        double SemiMajor;
        /// semi-minor of the ellipse
        double SemiMinor;
        /// number of nodes in the xy plane (for this case: independent of the z coordinate)
        int nxy_m;
        /// interpolation type
        int interpolationMethod;

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

        /// different interpolation methods for boundary points 
        void ConstantInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor);
        void LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor);
        void QuadraticInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor);

};

#endif //#ifdef HAVE_ML_SOLVER
#endif //#ifdef ELLIPTICAL_DOMAIN_H
