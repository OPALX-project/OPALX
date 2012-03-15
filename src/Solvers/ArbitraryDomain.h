/* can handle arbitrary domains
 * at the moment implementation for validation against solution with zylinder on [0,1]
 */

#ifndef ARBITRARY_DOMAIN
#define ARBITRARY_DOMAIN
#ifdef HAVE_ML_SOLVER

#include "IrregularDomain.h"
#include <mpi.h>
#include <hdf5.h>
#include "H5Part.h"
#include "H5Fed.h"
#include <vector>
#include <map>
#include <math.h>

struct vertex {
  h5_float64_t P[3];
};

typedef struct vertex vertex_t; 

struct entity {
  h5_id_t global_id;
  h5_id_t parent_id;
  h5_id_t vertexIDs[3];
};
typedef struct entity entity_t;

///PointList maps from an (x,z) resp. (y,z) pair to double values (=intersections with boundary) 
typedef std::multimap< std::pair<int,int>, double > PointList;

class ArbitraryDomain : public IrregularDomain {

public:

  ArbitraryDomain(string fname, Vector_t nr_, Vector_t hr_, NDIndex<3> locidx);

  /// load geometry file
  void LoadFile();
  /// calculates intersection with the elliptic beam pipe
  void Compute(Vector_t hr);
  /// returns number of nodes in xy plane
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
  string getType() {return "Geometric";}
  /// queries if a given (x,y,z) coordinate lies inside the domain
  inline bool isInside(int x, int y, int z);

	/// set a filename where the geometry residues
  void setFilename(string fname) {filename = fname;}

  int getStartIdx() {return startIdx;}
  /// conversion from (x,y,z) to index on the 3D grid
  int getIdx(int x, int y, int z); 

private:

  /// all intersection points with gridlines in Y direction
  PointList IntersectYDir;
  /// all intersection points with gridlines in X direction
  PointList IntersectXDir;
  /// all intersection points with gridlines in Z direction
  PointList IntersectZDir;
  /// filepath to the mesh
  string filename;
  /// vector storing the triangles
  std::map<h5_id_t, vertex_t> vertices;
  std::map<h5_id_t, entity_t> entities;

  int startIdx;
  NDIndex<3> localidx;

	/// here we store the number of nodes in a xy layer for a given z coordinate
	std::map<int, int> numXY;

  /// mapping (x,y,z) -> idx
  std::map<int, int> IdxMap;
  /// mapping idx -> (x,y,z)
  std::map<int, int> CoordMap;
	
  /// conversion from (x,y,z) to index in xyz plane
  inline int toCoordIdx(int x, int y, int z); 
  // conversion from (x,y,z) to index on the 3D grid
  //inline int getIdx(int x, int y, int z);
  /// conversion from a 3D index to (x,y,z)
  inline void getCoord(int idx, int& x, int& y, int& z);
  
  inline void crossProduct(double A[], double B[], double C[]);
  inline double dotProduct(double v1[], double v2[]) { return  (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]); }

};

#endif //#ifdef HAVE_ML_SOLVER
#endif //#ifdef ARBITRARY_DOMAIN
