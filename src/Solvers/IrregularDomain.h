#ifndef IRREGULAR_DOMAIN_H
#define IRREGULAR_DOMAIN_H
#ifdef HAVE_ML_SOLVER

#include <vector>
#include <string>
#include "Ippl.h"

typedef UniformCartesian<3, double> Mesh_t;
typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

/// enumeration corresponding to different interpolation methods at the boundary
enum {
    CONSTANT,
    LINEAR,
    QUADRATIC
};

/// this class defines a common abstract interface for different types of boundaries
class IrregularDomain {

public:

    /** method to compute the intersection points with the boundary geometry (stored in some appropriate data structure)
     * \param hr updated mesh spacings
     */
    virtual void Compute(Vector_t hr) = 0;

    /** method to get the number of gridpoints in a given z plane
     * \param z coordinate of the z plane
     * \return int number of grid nodes in the given z plane
     */
    virtual int getNumXY(int z) = 0;

    /// method to calculate the stencil at a boundary points
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \param W stencil value of the element in the west of idx: (x-1)
    /// \param E stencil value of the element in the east of idx: (x+1)
    /// \param S stencil value of the element in the south of idx: (y-1)
    /// \param N stencil value of the element in the north of idx: (y+1)
    /// \param F stencil value of the element in front of idx: (z-1)
    /// \param B stencil value of the element in the back of idx: (z+1)
    /// \param C stencil value of the element in the center
    virtual void getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) = 0;

    /// method to calculate the stencil at a boundary points
    /// \param idx index of the current element in the matrix
    /// \param W stencil value of the element in the west of idx: (x-1)
    /// \param E stencil value of the element in the east of idx: (x+1)
    /// \param S stencil value of the element in the south of idx: (y-1)
    /// \param N stencil value of the element in the north of idx: (y+1)
    /// \param F stencil value of the element in front of idx: (z-1)
    /// \param B stencil value of the element in the back of idx: (z+1)
    /// \param C stencil value of the element in the center
    virtual void getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) = 0;

    /// method to calculate the neighbours in the matrix of the current index (x,y,z)
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \param W stencil index of the element in the west of idx: (x-1)
    /// \param E stencil index of the element in the east of idx: (x+1)
    /// \param S stencil index of the element in the south of idx: (y-1)
    /// \param N stencil index of the element in the north of idx: (y+1)
    /// \param F stencil index of the element in front of idx: (z-1)
    /// \param B stencil index of the element in the back of idx: (z+1)
    virtual void getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B) = 0;
    virtual void getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B) = 0;

    /// method that identifies a specialized boundary geometry
    /// \return std::string containing a description of the boundary geometry used
    virtual std::string getType() = 0;

    /// method that checks if a given point lies inside the boundary
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \return boolean indicating if the point lies inside the boundary
    virtual bool isInside(int x, int y, int z) = 0;

    Vector_t getNr() { return nr; }
    Vector_t getHr() { return hr; }
    void setNr(Vector_t nri) { nr = nri; }
    void setHr(Vector_t hri) { hr = hri; }

    void setMinMaxZ(double minz, double maxz) { zMin_m=minz; zMax_m=maxz; }
    double getMinZ() { return zMin_m; }
    double getMaxZ() { return zMax_m; }

    virtual double getXRangeMin() = 0;
    virtual double getXRangeMax() = 0;
    virtual double getYRangeMin() = 0;
    virtual double getYRangeMax() = 0;

    virtual int getIdx(int x, int y, int z) = 0;
    virtual bool hasGeometryChanged() = 0;

protected:

    // a irregular domain is always defined on a grid
    /// number of mesh points in each direction
    Vector_t nr;
    /// mesh-spacings in each direction
    Vector_t hr;

    /// min/max of bunch in floor coordinates
    double zMin_m;
    double zMax_m;
};

#endif //#ifdef HAVE_ML_SOLVER
#endif //#ifndef IRREGULAR_DOMAIN_H
