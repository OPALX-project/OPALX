/*
  Implementation of the class BoundaryGeometry.

  Copyright & License: See Copyright.readme in src directory
 */

#include <fstream>

#include "H5hut.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/PriEmissionPhysics.h"
#include "Expressions/SRefExpr.h"
#include "Elements/OpalBeamline.h"

extern Inform* gmsg;

#define SQR(x) ((x)*(x))
#define PointID(triangle_id, vertex_id) allbfaces_m[4 * (triangle_id) + (vertex_id)]
#define Point(triangle_id, vertex_id)   geo3Dcoords_m[allbfaces_m[4 * (triangle_id) + (vertex_id)]]

#define FAST_VOXELIZATION

/*

  Some
   _   _      _                 
  | | | | ___| |_ __   ___ _ __ 
  | |_| |/ _ \ | '_ \ / _ \ '__|
  |  _  |  __/ | |_) |  __/ |   
  |_| |_|\___|_| .__/ \___|_|   
             |_|  

  functions
 */
struct VectorLessX {
    bool operator() (Vector_t x1, Vector_t x2) {
        return x1 (0) < x2 (0);
    }
};

struct VectorLessY {
    bool operator() (Vector_t x1, Vector_t x2) {
        return x1 (1) < x2 (1);
    }
};

struct VectorLessZ {
    bool operator() (Vector_t x1, Vector_t x2) {
        return x1 (2) < x2 (2);
    }
};

/**
   Calculate the maximum of coordinates of geometry,i.e the maximum of X,Y,Z
 */
Vector_t get_max_extend (std::vector<Vector_t>& coords) {
    const Vector_t x = *max_element (
        coords.begin (), coords.end (), VectorLessX ());
    const Vector_t y = *max_element (
        coords.begin (), coords.end (), VectorLessY ());
    const Vector_t z = *max_element (
        coords.begin (), coords.end (), VectorLessZ ());
    return Vector_t (x (0), y (1), z (2));
}

/*
   Compute the minimum of coordinates of geometry, i.e the minimum of X,Y,Z
 */
Vector_t get_min_extend (std::vector<Vector_t>& coords) {
    const Vector_t x = *min_element (
        coords.begin (), coords.end (), VectorLessX ());
    const Vector_t y = *min_element (
        coords.begin (), coords.end (), VectorLessY ());
    const Vector_t z = *min_element (
        coords.begin (), coords.end (), VectorLessZ ());
    return Vector_t (x (0), y (1), z (2));
}

/*
  "ULP compare" for double precision floating point numbers.
  See:
    http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

  Note:
    An updated version of this document with improved code is here:
    http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/

 */
static int64_t fcmp (
    const double A,
    const double B,
    const int maxUlps ) {
                    
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    assert (maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    assert (sizeof (long long) == sizeof (int64_t) );
    assert (sizeof (long long) == sizeof (double) );
                    
    // Make [ab]Int lexicographically ordered as a twos-complement int
    const double* pa = &A;
    int64_t aInt = *(int64_t*)pa;
    if (aInt < 0)
        aInt = 0x8000000000000000LL - aInt;
                    
    const double* pb = &B;
    int64_t bInt = *(int64_t*)pb;
    if (bInt < 0)
        bInt = 0x8000000000000000LL - bInt;
                    
    const int64_t intDiff = aInt - bInt;
    if (llabs(intDiff) <= maxUlps)
        return 0;
    return intDiff;
}

static inline bool is_in_voxel (Vector_t& point, Vector_t& min, Vector_t& max) {
    if (fcmp (point [0], min [0], 10) < 0) return false;
    if (fcmp (point [1], min [1], 10) < 0) return false;
    if (fcmp (point [2], min [2], 10) < 0) return false;
    if (fcmp (point [0], max [0], 10) > 0) return false;
    if (fcmp (point [1], max [1], 10) > 0) return false;
    if (fcmp (point [2], max [2], 10) > 0) return false;
    return true;
}

/*
  write legacy VTK file of voxel mesh
*/
static void write_voxel_mesh (
    const std::unordered_set<int> ids,
    const Vector_t hr_m,
    const Vektor<int,3> nr,
    const Vector_t origin
    ) {
    /*----------------------------------------------------------------------*/
    const size_t numpoints = 8 * ids.size ();
    std::ofstream of;
    of.open (string ("data/testBBox.vtk").c_str ());
    assert (of.is_open ());
    of.precision (6);
    
    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated using BoundaryGeometry::computeMeshVoxelization"
       << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
    of << "POINTS " << numpoints << " float" << std::endl;

    const auto end_it = ids.end();
    const auto nr0_times_nr1 = nr[0] * nr[1];
    for (auto id = ids.begin (); id != end_it; id++) {
        int k = (*id - 1) / nr0_times_nr1;
        int rest = (*id - 1) % nr0_times_nr1;
        int j = rest / nr[0];
        int i = rest % nr[0]; 

        Vector_t P;
        P[0] = i * hr_m[0] + origin[0];
        P[1] = j * hr_m[1] + origin[1];
        P[2] = k * hr_m[2] + origin[2];

        of << P[0]           << " " << P[1]           << " " << P[2]           << std::endl;
        of << P[0] + hr_m[0] << " " << P[1]           << " " << P[2]           << std::endl;
        of << P[0]           << " " << P[1] + hr_m[1] << " " << P[2]           << std::endl;
        of << P[0] + hr_m[0] << " " << P[1] + hr_m[1] << " " << P[2]           << std::endl;
        of << P[0]           << " " << P[1]           << " " << P[2] + hr_m[2] << std::endl;
        of << P[0] + hr_m[0] << " " << P[1]           << " " << P[2] + hr_m[2] << std::endl;
        of << P[0]           << " " << P[1] + hr_m[1] << " " << P[2] + hr_m[2] << std::endl;
        of << P[0] + hr_m[0] << " " << P[1] + hr_m[1] << " " << P[2] + hr_m[2] << std::endl;
    }
    of << std::endl;
    const auto num_cells = ids.size ();
    of << "CELLS " << num_cells << " " << 9 * num_cells << std::endl;
    for (size_t i = 0; i < num_cells; i++)
        of << "8 "
           << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3 << " "
           << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << std::endl;
    of << "CELL_TYPES " << num_cells << std::endl;
    for (size_t i = 0; i <  num_cells; i++)
        of << "11" << std::endl;
    of << "CELL_DATA " << num_cells << std::endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
    of << "LOOKUP_TABLE " << "default" << std::endl;
    for (size_t i = 0; i <  num_cells; i++)
        of << (float)(i) << std::endl;
    of << std::endl;
    of << "COLOR_SCALARS " << "BBoxColor " << 4 << std::endl;
    for (size_t i = 0; i < num_cells; i++) {
        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << std::endl;
    }
    of << std::endl;
}

/*___________________________________________________________________________

  Triangle-cube intersection test.

  See:
  http://tog.acm.org/resources/GraphicsGems/gemsiii/triangleCube.c

 */

#include <math.h>


#define LERP( A, B, C) ((B)+(A)*((C)-(B)))
#define MIN3(a,b,c) ((((a)<(b))&&((a)<(c))) ? (a) : (((b)<(c)) ? (b) : (c)))
#define MAX3(a,b,c) ((((a)>(b))&&((a)>(c))) ? (a) : (((b)>(c)) ? (b) : (c)))
#define INSIDE 0
#define OUTSIDE 1

typedef struct {
    Vector_t v1;                 /* Vertex1 */
    Vector_t v2;                 /* Vertex2 */
    Vector_t v3;                 /* Vertex3 */
} Triangle; 

typedef struct {
    Vector_t v1;
    Vector_t v2;
} Cube;

/*___________________________________________________________________________*/

/* Which of the six face-plane(s) is point P outside of? */

static inline int
face_plane (
    const Vector_t p
    ) {
    int outcode;

    outcode = 0;
    if (p[0] >  .5) outcode |= 0x01;
    if (p[0] < -.5) outcode |= 0x02;
    if (p[1] >  .5) outcode |= 0x04;
    if (p[1] < -.5) outcode |= 0x08;
    if (p[2] >  .5) outcode |= 0x10;
    if (p[2] < -.5) outcode |= 0x20;
    return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Which of the twelve edge plane(s) is point P outside of? */

static inline int
bevel_2d (
    const Vector_t p
    ) {
    int outcode;

    outcode = 0;
    if ( p[0] + p[1] > 1.0) outcode |= 0x001;
    if ( p[0] - p[1] > 1.0) outcode |= 0x002;
    if (-p[0] + p[1] > 1.0) outcode |= 0x004;
    if (-p[0] - p[1] > 1.0) outcode |= 0x008;
    if ( p[0] + p[2] > 1.0) outcode |= 0x010;
    if ( p[0] - p[2] > 1.0) outcode |= 0x020;
    if (-p[0] + p[2] > 1.0) outcode |= 0x040;
    if (-p[0] - p[2] > 1.0) outcode |= 0x080;
    if ( p[1] + p[2] > 1.0) outcode |= 0x100;
    if ( p[1] - p[2] > 1.0) outcode |= 0x200;
    if (-p[1] + p[2] > 1.0) outcode |= 0x400;
    if (-p[1] - p[2] > 1.0) outcode |= 0x800;
    return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

  Which of the eight corner plane(s) is point P outside of?
*/
static inline int
bevel_3d (
    const Vector_t p
    ) {
    int outcode;

    outcode = 0;
    if (( p[0] + p[1] + p[2]) > 1.5) outcode |= 0x01;
    if (( p[0] + p[1] - p[2]) > 1.5) outcode |= 0x02;
    if (( p[0] - p[1] + p[2]) > 1.5) outcode |= 0x04;
    if (( p[0] - p[1] - p[2]) > 1.5) outcode |= 0x08;
    if ((-p[0] + p[1] + p[2]) > 1.5) outcode |= 0x10;
    if ((-p[0] + p[1] - p[2]) > 1.5) outcode |= 0x20;
    if ((-p[0] - p[1] + p[2]) > 1.5) outcode |= 0x40;
    if ((-p[0] - p[1] - p[2]) > 1.5) outcode |= 0x80;
    return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Test the point "alpha" of the way from P1 to P2
  See if it is on a face of the cube
  Consider only faces in "mask"
*/

static inline int
check_point (
    const Vector_t p1,
    const Vector_t p2,
    const double alpha,
    const int mask
    ) {
    Vector_t plane_point;

    plane_point[0] = LERP(alpha, p1[0], p2[0]);
    plane_point[1] = LERP(alpha, p1[1], p2[1]);
    plane_point[2] = LERP(alpha, p1[2], p2[2]);
    return(face_plane(plane_point) & mask);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Compute intersection of P1 --> P2 line segment with face planes
  Then test intersection point to see if it is on cube face
  Consider only face planes in "outcode_diff"
  Note: Zero bits in "outcode_diff" means face line is outside of
*/
static inline int
check_line (
    const Vector_t p1,
    const Vector_t p2,
    const int outcode_diff
    ) {
    if ((0x01 & outcode_diff) != 0)
        if (check_point(p1,p2,( .5-p1[0])/(p2[0]-p1[0]),0x3e) == INSIDE) return(INSIDE);
    if ((0x02 & outcode_diff) != 0)
        if (check_point(p1,p2,(-.5-p1[0])/(p2[0]-p1[0]),0x3d) == INSIDE) return(INSIDE);
    if ((0x04 & outcode_diff) != 0) 
        if (check_point(p1,p2,( .5-p1[1])/(p2[1]-p1[1]),0x3b) == INSIDE) return(INSIDE);
    if ((0x08 & outcode_diff) != 0) 
        if (check_point(p1,p2,(-.5-p1[1])/(p2[1]-p1[1]),0x37) == INSIDE) return(INSIDE);
    if ((0x10 & outcode_diff) != 0) 
        if (check_point(p1,p2,( .5-p1[2])/(p2[2]-p1[2]),0x2f) == INSIDE) return(INSIDE);
    if ((0x20 & outcode_diff) != 0) 
        if (check_point(p1,p2,(-.5-p1[2])/(p2[2]-p1[2]),0x1f) == INSIDE) return(INSIDE);
    return(OUTSIDE);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Test if 3D point is inside 3D triangle
*/

#define EPS 10e-5

static inline int
SIGN3 (
    Vector_t A
    ) {
    return ((A[0] < EPS) ? 4 : 0 | (A[0] > -EPS) ? 32 : 0 |
            (A[1] < EPS) ? 2 : 0 | (A[1] > -EPS) ? 16 : 0 |
            (A[2] < EPS) ? 1 : 0 | (A[2] > -EPS) ? 8 : 0);
}

static int
point_triangle_intersection (
    const Vector_t p,
    const Triangle t
    ) {
    /*
      First, a quick bounding-box test:
      If P is outside triangle bbox, there cannot be an intersection.
    */
    if (p[0] > MAX3(t.v1[0], t.v2[0], t.v3[0])) return(OUTSIDE);  
    if (p[1] > MAX3(t.v1[1], t.v2[1], t.v3[1])) return(OUTSIDE);
    if (p[2] > MAX3(t.v1[2], t.v2[2], t.v3[2])) return(OUTSIDE);
    if (p[0] < MIN3(t.v1[0], t.v2[0], t.v3[0])) return(OUTSIDE);
    if (p[1] < MIN3(t.v1[1], t.v2[1], t.v3[1])) return(OUTSIDE);
    if (p[2] < MIN3(t.v1[2], t.v2[2], t.v3[2])) return(OUTSIDE);
    
    /*
      For each triangle side, make a vector out of it by subtracting vertexes;
      make another vector from one vertex to point P.
      The crossproduct of these two vectors is orthogonal to both and the
      signs of its X,Y,Z components indicate whether P was to the inside or
      to the outside of this triangle side.                                
    */
    const Vector_t vect12 = t.v1 - t.v2;
    const Vector_t vect1h = t.v1 - p;
    const Vector_t cross12_1p = cross (vect12, vect1h);
    const int sign12 = SIGN3(cross12_1p);      /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

    const Vector_t vect23 = t.v2 - t.v3;
    const Vector_t vect2h = t.v2 - p;
    const Vector_t cross23_2p = cross (vect23, vect2h);
    const int sign23 = SIGN3(cross23_2p);

    const Vector_t vect31 = t.v3 - t.v1;
    const Vector_t vect3h = t.v3 - p;
    const Vector_t cross31_3p = cross (vect31, vect3h);
    const int sign31 = SIGN3(cross31_3p);

    /*
      If all three crossproduct vectors agree in their component signs,
      then the point must be inside all three.
      P cannot be OUTSIDE all three sides simultaneously.
    */
    return ((sign12 & sign23 & sign31) == 0) ? OUTSIDE : INSIDE;
}


/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  This is the main algorithm procedure.
  Triangle t is compared with a unit cube,
  centered on the origin.
  It returns INSIDE (0) or OUTSIDE(1) if t
  intersects or does not intersect the cube.
*/
static int
t_c_intersection (
    Triangle t
    ) {
    int v1_test;
    int v2_test;
    int v3_test;

    /*
      First compare all three vertexes with all six face-planes
      If any vertex is inside the cube, return immediately!
    */
   if ((v1_test = face_plane(t.v1)) == INSIDE) return(INSIDE);
   if ((v2_test = face_plane(t.v2)) == INSIDE) return(INSIDE);
   if ((v3_test = face_plane(t.v3)) == INSIDE) return(INSIDE);

   /*
     If all three vertexes were outside of one or more face-planes,
     return immediately with a trivial rejection!
   */
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

   /*
     Now do the same trivial rejection test for the 12 edge planes
   */
   v1_test |= bevel_2d(t.v1) << 8; 
   v2_test |= bevel_2d(t.v2) << 8; 
   v3_test |= bevel_2d(t.v3) << 8;
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);  

   /*
     Now do the same trivial rejection test for the 8 corner planes
   */
   v1_test |= bevel_3d(t.v1) << 24; 
   v2_test |= bevel_3d(t.v2) << 24; 
   v3_test |= bevel_3d(t.v3) << 24; 
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);   

   /*
     If vertex 1 and 2, as a pair, cannot be trivially rejected
     by the above tests, then see if the v1-->v2 triangle edge
     intersects the cube.  Do the same for v1-->v3 and v2-->v3./
     Pass to the intersection algorithm the "OR" of the outcode
     bits, so that only those cube faces which are spanned by
     each triangle edge need be tested.
   */
   if ((v1_test & v2_test) == 0)
      if (check_line(t.v1,t.v2,v1_test|v2_test) == INSIDE) return(INSIDE);
   if ((v1_test & v3_test) == 0)
      if (check_line(t.v1,t.v3,v1_test|v3_test) == INSIDE) return(INSIDE);
   if ((v2_test & v3_test) == 0)
      if (check_line(t.v2,t.v3,v2_test|v3_test) == INSIDE) return(INSIDE);

   /*
     By now, we know that the triangle is not off to any side,
     and that its sides do not penetrate the cube.  We must now
     test for the cube intersecting the interior of the triangle.
     We do this by looking for intersections between the cube
     diagonals and the triangle...first finding the intersection
     of the four diagonals with the plane of the triangle, and
     then if that intersection is inside the cube, pursuing
     whether the intersection point is inside the triangle itself.

     To find plane of the triangle, first perform crossproduct on 
     two triangle side vectors to compute the normal vector.
   */  
   Vector_t vect12 = t.v1 - t.v2;
   Vector_t vect13 = t.v1 - t.v3;
   Vector_t norm = cross (vect12, vect13);

   /*
     The normal vector "norm" X,Y,Z components are the coefficients
     of the triangles AX + BY + CZ + D = 0 plane equation.  If we
     solve the plane equation for X=Y=Z (a diagonal), we get
     -D/(A+B+C) as a metric of the distance from cube center to the
     diagonal/plane intersection.  If this is between -0.5 and 0.5,
     the intersection is inside the cube.  If so, we continue by
     doing a point/triangle intersection.
     Do this for all four diagonals.
   */
   double d = norm[0] * t.v1[0] + norm[1] * t.v1[1] + norm[2] * t.v1[2];

   /*
     if one of the diagonals is parallel to the plane, the other will
     intersect the plane
   */
   double denom;
   if(fabs(denom=(norm[0] + norm[1] + norm[2]))>EPS) {
       /* skip parallel diagonals to the plane; division by 0 can occur */
       Vector_t hitpp = d / denom;
       if (fabs(hitpp[0]) <= 0.5)
           if (point_triangle_intersection(hitpp,t) == INSIDE) return(INSIDE);
   }
   if(fabs(denom=(norm[0] + norm[1] - norm[2]))>EPS) {
       Vector_t hitpn;
       hitpn[2] = -(hitpn[0] = hitpn[1] = d / denom);
       if (fabs(hitpn[0]) <= 0.5)
           if (point_triangle_intersection(hitpn,t) == INSIDE) return(INSIDE);
   }       
   if(fabs(denom=(norm[0] - norm[1] + norm[2]))>EPS) {       
       Vector_t hitnp;
       hitnp[1] = -(hitnp[0] = hitnp[2] = d / denom);
       if (fabs(hitnp[0]) <= 0.5)
           if (point_triangle_intersection(hitnp,t) == INSIDE) return(INSIDE);
   }
   if(fabs(denom=(norm[0] - norm[1] - norm[2]))>EPS) {
       Vector_t hitnn;
       hitnn[1] = hitnn[2] = -(hitnn[0] = d / denom);
       if (fabs(hitnn[0]) <= 0.5)
           if (point_triangle_intersection(hitnn,t) == INSIDE) return(INSIDE);
   }
   
   /*
     No edge touched the cube; no cube diagonal touched the triangle.
     We're done...there was no intersection.
   */
   return(OUTSIDE);
}

static inline void
scaleCube (
    Cube& c,
    const Vector_t& scale
    ) {
    c.v1[0] *= scale[0];
    c.v1[1] *= scale[1];
    c.v1[2] *= scale[2];
    c.v2[0] *= scale[0];
    c.v2[1] *= scale[1];
    c.v2[2] *= scale[2];
}

static inline void
scaleTriangle (
    Triangle& t,
    const Vector_t& scale
    ) {
    t.v1[0] *= scale[0];
    t.v1[1] *= scale[1];
    t.v1[2] *= scale[2];
    t.v2[0] *= scale[0];
    t.v2[1] *= scale[1];
    t.v2[2] *= scale[2];
    t.v3[0] *= scale[0];
    t.v3[1] *= scale[1];
    t.v3[2] *= scale[2];
}

static inline void
shiftTriangle (
    Triangle& t,
    const Vector_t& shift
    ) {
    t.v1 -= shift;
    t.v2 -= shift;
    t.v3 -= shift;
}

int
intersect3dTriangleCube (
    const Triangle& t,
    const Cube& c
    ) {
    Cube c_ = c;
    Triangle t_ = t;
    Vector_t extend = c_.v2 - c_.v1;
    const Vector_t scale = 1.0 / extend; 
    scaleTriangle (t_, scale);
    scaleCube (c_, scale);

    shiftTriangle (t_, c_.v1 - 0.5);
    
    return t_c_intersection (t);
}

/*
  ____                           _              
 / ___| ___  ___  _ __ ___   ___| |_ _ __ _   _ 
| |  _ / _ \/ _ \| '_ ` _ \ / _ \ __| '__| | | |
| |_| |  __/ (_) | | | | | |  __/ |_| |  | |_| |
 \____|\___|\___/|_| |_| |_|\___|\__|_|   \__, |
                                          |___/
*/

BoundaryGeometry::BoundaryGeometry() :
    Definition (
        SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry."),
    Tribarycent_m (NULL),
    TriPrPartloss_m (NULL),
    TriSePartloss_m (NULL),
    TriFEPartloss_m (NULL),
    allbfaces_m (NULL) {

    itsAttr[FGEOM] = Attributes::makeString
        ("FGEOM",
         "Specifies the geometry file [h5fed]",
         "");

    itsAttr[TOPO] = Attributes::makeString
        ("TOPO",
         "BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written ",
         "ELLIPTIC");

    itsAttr[LENGTH] = Attributes::makeReal
        ("LENGTH",
         "Specifies the length of a tube shaped elliptic beam pipe [m]",
         1.0);

    itsAttr[S] = Attributes::makeReal
        ("S",
         "Specifies the start of a tube shaped elliptic beam pipe [m]",
         0.0);

    itsAttr[A] = Attributes::makeReal
        ("A",
         "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]",
         0.025);

    itsAttr[B] = Attributes::makeReal
        ("B",
         "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]",
         0.025);

    itsAttr[L1] = Attributes::makeReal
        ("L1",
         "In case of BOXCORNER Specifies first part with hight == B [m]",
         0.5);

    itsAttr[L2] = Attributes::makeReal
        ("L2",
         "In case of BOXCORNER Specifies first second with hight == B-C [m]",
         0.2);

    itsAttr[C] = Attributes::makeReal
        ("C",
         "In case of BOXCORNER Specifies hight of corner C [m]",
         0.01);

    itsAttr[DISTR] = Attributes::makeString
        ("DISTR",
         "Distribution to be generated on the surface",
         "");
    itsAttr[DISTRS] = Attributes::makeStringArray
        ("DISTRS",
         "Distribution array to be generated on the surface");

    itsAttr[XYZSCALE] = Attributes::makeReal
        ("XYZSCALE",
         "Multiplicative scaling factor for coordinates ",
         1.0);

    itsAttr[XSCALE] = Attributes::makeReal
        ("XSCALE",
         "Multiplicative scaling factor for X coordinates ",
         1.0);

    itsAttr[YSCALE] = Attributes::makeReal
        ("YSCALE",
         "Multiplicative scaling factor for Y coordinates ",
         1.0);

    itsAttr[ZSCALE] = Attributes::makeReal
        ("ZSCALE",
         "Multiplicative scaling factor for Z coordinates ",
         1.0);

    itsAttr[ZSHIFT] = Attributes::makeReal
        ("ZSHIFT",
         "Shift in z direction",
         0.0);

    itsAttr[APERTURE]  = Attributes::makeRealArray
        ("APERTURE", "The element aperture");

    BoundaryGeometry* defGeometry = clone ("UNNAMED_GEOMETRY");
    defGeometry->builtin = true;

    TPInside_m = IpplTimings::getTimer ("Particle Inside");
    TPreProc_m = IpplTimings::getTimer ("Pre Processing");
    TRayTrace_m = IpplTimings::getTimer ("Ray tracing");
    h5FileName_m = Attributes::getString (itsAttr[FGEOM]);
    Tinward_m = IpplTimings::getTimer ("Check inward");
    try {
        defGeometry->update ();
        OpalData::getInstance ()->define (defGeometry);
    } catch (...) {
        delete defGeometry;
    }
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);    

    if (!h5FileName_m.empty ())
        initialize ();

}

BoundaryGeometry::BoundaryGeometry(
    const string& name,
    BoundaryGeometry* parent
    ) :
    Definition (name, parent),
    Tribarycent_m (NULL),
    TriPrPartloss_m (NULL),
    TriSePartloss_m (NULL),
    TriFEPartloss_m (NULL),
    allbfaces_m (NULL) {
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);    

    h5FileName_m = Attributes::getString (itsAttr[FGEOM]);
    if (!h5FileName_m.empty ())
        initialize ();
    TPInside_m = IpplTimings::getTimer ("Particle Inside");
    TPreProc_m = IpplTimings::getTimer ("zshiftPre Processing");
    TRayTrace_m = IpplTimings::getTimer ("Ray tracing");
    Tinward_m = IpplTimings::getTimer ("Check inward");

 }

BoundaryGeometry::~BoundaryGeometry() {
    /*
       if (allbfaces_m)
          delete allbfaces_m;
       if (Tribarycent_m)
          delete Tribarycent_m;
       if (TriPrPartloss_m )
       delete TriPrPartloss_m ;
       if (TriFEPartloss_m)
       delete TriFEPartloss_m;
       if (TriSePartloss_m)
       delete TriSePartloss_m;
     */

    gsl_rng_free(randGen_m);
}

bool BoundaryGeometry::canReplaceBy (Object* object) {
    // Can replace only by another GEOMETRY.
    return dynamic_cast<BGeometryBase*>(object) != 0;
}

BoundaryGeometry* BoundaryGeometry::clone (const string& name) {
    return new BoundaryGeometry (name, this);
}

void BoundaryGeometry::update () {
    if (getOpalName ().empty ()) setOpalName ("UNNAMED_GEOMETRY");
}


void BoundaryGeometry::execute () {
    update ();
    TPInside_m = IpplTimings::getTimer ("Particle Inside");
    TPreProc_m = IpplTimings::getTimer ("Pre Processing");
    TRayTrace_m = IpplTimings::getTimer ("Ray tracing");
    Tinward_m = IpplTimings::getTimer ("Check inward");
}

BoundaryGeometry* BoundaryGeometry::find (const string& name) {
    BoundaryGeometry* geom = dynamic_cast<BoundaryGeometry*>(
        OpalData::getInstance ()->find (name));

    if (geom == 0)
        throw OpalException ("BoundaryGeometry::find()", "Geometry \""
                             + name + "\" not found.");
    return geom;
}

void BoundaryGeometry::updateElement (ElementBase* element) {
}

int
BoundaryGeometry::intersectTriangleVoxel (
    const int triangle_id,
    const int i,
    const int j,
    const int k
    ) {
    const Triangle t = {getPoint (triangle_id, 1),
                        getPoint (triangle_id, 2),
                        getPoint (triangle_id, 3)};

    //const int k = (voxel_id - 1) / (nr_m[0] * nr_m[1]);
    //const int rest = (voxel_id - 1) % (nr_m[0] * nr_m[1]);
    //const int j = rest / nr_m[0];
    //const int i = rest % nr_m[0]; 
    
    Cube c;
    c.v1 = {
        i * hr_m[0] + voxelMesh_m.minExtend[0],
        j * hr_m[1] + voxelMesh_m.minExtend[1],
        k * hr_m[2] + voxelMesh_m.minExtend[2]
    };

    c.v2 = c.v1 + hr_m;

    return intersect3dTriangleCube (t, c);
}

/*
  Find the 3D intersection of a line segment, ray or line with a triangle.

  Input:
        kind: type of test: SEGMENT, RAY or LINE
        P0, P0: defining 
            a line segment from P0 to P1 or
            a ray starting at P0 with directional vector P1-P0 or
            a line through P0 and P1
        V0, V1, V2: the triangle vertices
        
  Output:
        I: intersection point (when it exists)

  Return values for line segment and ray test :
        -1 = triangle is degenerated (a segment or point)
        0 =  disjoint (no intersect)
        1 =  are in the same plane
        2 =  intersect in unique point I1

  Return values for line intersection test :
        -1: triangle is degenerated (a segment or point)
        0:  disjoint (no intersect)
        1:  are in the same plane
        2:  intersect in unique point I1, with r < 0.0
        3:  intersect in unique point I1, with 0.0 <= r <= 1.0
        4:  intersect in unique point I1, with 1.0 < r

  For algorithm and implementation see:
  http://geomalgorithms.com/a06-_intersect-2.html

  Copyright 2001 softSurfer, 2012 Dan Sunday
  This code may be freely used and modified for any purpose
  providing that this copyright notice is included with it.
  SoftSurfer makes no warranty for this code, and cannot be held
  liable for any real or imagined damage resulting from its use.
  Users of this code must verify correctness for their application.
 */

int
BoundaryGeometry::intersectLineTriangle (
    const enum INTERSECTION_TESTS kind,
    const Vector_t& P0,
    const Vector_t& P1,
    const int triangle_id,
    Vector_t& I
    ) {
    const Vector_t V0 = getPoint (triangle_id, 1);
    const Vector_t V1 = getPoint (triangle_id, 2);
    const Vector_t V2 = getPoint (triangle_id, 3);

    // get triangle edge vectors and plane normal
    const Vector_t u = V1 - V0;         // triangle vectors
    const Vector_t v = V2 - V0;
    const Vector_t n = cross (u, v);
    if (n == (Vector_t)0)               // triangle is degenerate
        return -1;                      // do not deal with this case
    
    const Vector_t dir = P1 - P0;       // ray direction vector
    const Vector_t w0 = P0 - V0;
    const double a = -dot(n,w0);
    const double b = dot(n,dir);
    if (fcmp (b, 0.0, 10) == 0) {       // ray is  parallel to triangle plane
        if (a == 0) {                   // ray lies in triangle plane
            return 1;
        } else {                        // ray disjoint from plane
            return 0;
        }
    }
    
    // get intersect point of ray with triangle plane
    const double r = a / b;
    switch (kind) {
    case RAY:
        if (r < 0.0) {                  // ray goes away from triangle
            return 0;                   // => no intersect
        }
    case SEGMENT:
        if (r < 0 || 1.0 < r) {         // intersection on extended
            return 0;                   // segment
        }
    case LINE:
        break;
    };
    I = P0 + r * dir;                   // intersect point of ray and plane
    
    // is I inside T?
    const double uu = dot(u,u);
    const double uv = dot(u,v);
    const double vv = dot(v,v);
    const Vector_t w = I - V0;
    const double wu = dot(w,u);
    const double wv = dot(w,v);
    const double D = uv * uv - uu * vv;
    
    // get and test parametric coords
    const double s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) {           // I is outside T
        return 0;
    }
    const double t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0) {     // I is outside T
        return 0;
    }
    // intersection point is in triangle
    if (r < 0.0) {                      // in extended segment in opposite
        return 2;                       // direction of ray
    } else if ((0.0 <= r) && (r <= 1.0)) { // in segment
        return 3;
    } else {                            // in extended segment in
        return 4;                       // direction of ray 
    }
}

static inline double magnitude (
    const Vector_t& v
    ) {
    return sqrt (dot (v,v));
}

/*
  Game plan:
  Count number of intersection of the line segment defined by P and a reference
  pt with the boundary. If the reference pt is inside the boundary and the number
  of intersections is even, then P is inside the geometry. Otherwise P is outside.
  To count the number of intersection, we divide the line segment in N segments
  and run the line-segment boundary intersection test for all these segments.
  N must be choosen carefully. It shouldn't be to large to avoid needless test.
 */
int BoundaryGeometry::fastIsInside (
    const Vector_t reference_pt,        // [in] a reference point which must be in the boundary
    const Vector_t P                    // [in] point to test
    ) {
    const Vector_t v = P - reference_pt;

    const int N = ceil (magnitude (v) / MIN3 (hr_m[0], hr_m[1], hr_m[2]));
    const Vector_t v_ = v / N;
    Vector_t P0 = P;
    Vector_t P1 = P + v_;
    Vector_t I;
    int triangle_id = -1;
    int result = 0;
    for (int i = 0; i < N; i++) {
        result += intersectLineSegmentBoundary (P0, P1, I, triangle_id) == 3 ? 1 : 0;
        P0 = P1;
        P1 += v_;
    }
    return result;
}

/*
  P must be *inside* the boundary geometry!
 */
int BoundaryGeometry::intersectRayBoundary (
    const Vector_t& P,
    const Vector_t& v,
    Vector_t& I) {

    const int N = ceil (magnitude (v) / MIN3 (hr_m[0], hr_m[1], hr_m[2]));
    const Vector_t v_ = v / N;
    Vector_t P0 = P;
    Vector_t P1 = P + v_;
    int triangle_id = -1;
    for (int i = 0; i < N; i++) {
        if (3 == intersectLineSegmentBoundary (P0, P1, I, triangle_id))
            return 1;
        P0 = P1;
        P1 += v_;
    }
    return 0;
}

/*
  Map point to unique voxel ID.

  Remember:
  * hr_m:  is the  mesh size
  * nr_m:  number of mesh points
  */
inline int
BoundaryGeometry::mapVoxelIndices2ID (
    const int i,
    const int j,
    const int k
    ) {
    if (i < 0 || i >= nr_m[0] ||
        j < 0 || j >= nr_m[1] ||
        k < 0 || k >= nr_m[2]) {
        return 0;
    }
    return 1 + k * nr_m[0] * nr_m[1] + j * nr_m[0] + i;
}

inline bool
BoundaryGeometry::mapPoint2VoxelIndices(
    const Vector_t pt,
    int& i,
    int& j,
    int& k
    ) {
    i = floor ((pt[0] - voxelMesh_m.minExtend [0]) / hr_m[0]);
    j = floor ((pt[1] - voxelMesh_m.minExtend [1]) / hr_m[1]);
    k = floor ((pt[2] - voxelMesh_m.minExtend [2]) / hr_m[2]);
    if (0 <= i && i < nr_m[0] &&
        0 <= j && j < nr_m[1] &&
        0 <= k && k < nr_m[2]) {
        return true;
    }
    return false;
}
    
inline int
BoundaryGeometry::mapPoint2VoxelID (
    const Vector_t x
    ) {
    int i, j, k;
    if (mapPoint2VoxelIndices (x, i, j, k)) {
        return mapVoxelIndices2ID (i, j, k);
    }
    return -1;
}

inline Vector_t&
BoundaryGeometry::mapPoint2Voxel (
    const Vector_t& pt
    ) {
    Vector_t r;
    const int i = floor ((pt[0] - voxelMesh_m.minExtend [0]) / hr_m[0]);
    const int j = floor ((pt[1] - voxelMesh_m.minExtend [1]) / hr_m[1]);
    const int k = floor ((pt[2] - voxelMesh_m.minExtend [2]) / hr_m[2]);

    r = {
        i * hr_m[0] + voxelMesh_m.minExtend[0],
        j * hr_m[1] + voxelMesh_m.minExtend[1],
        k * hr_m[2] + voxelMesh_m.minExtend[2]
    };
    Vector_t &result = r;
    return result;
}

void BoundaryGeometry::initialize () {

    class Local {

    public:

        static void computeGeometryInterval (BoundaryGeometry* bg) {

            bg->mincoords_m = get_min_extend (bg->geo3Dcoords_m);
            bg->maxcoords_m = get_max_extend (bg->geo3Dcoords_m);
            bg->len_m = bg->maxcoords_m - bg->mincoords_m;

            /*
              Calculate the maximum dimension of triangles. This value will be used to
              define the cubic box size
            */
            bg->longest_side_max_m = 0.0;
            bg->longest_side_min_m = 0.01;
            for (int i = 0; i < bg->num_triangles_m; i++) {
                // compute length of longest edge
                const Vector_t x1 = bg->getPoint (i, 1);
                const Vector_t x2 = bg->getPoint (i, 2);
                const Vector_t x3 = bg->getPoint (i, 3);
                const double length_edge1 = sqrt (
                    SQR (x1[0] - x2[0]) + SQR (x1[1] - x2[1]) + SQR (x1[2] - x2[2]));
                const double length_edge2 = sqrt (
                    SQR (x3[0] - x2[0]) + SQR (x3[1] - x2[1]) + SQR (x3[2] - x2[2]));
                const double length_edge3 = sqrt (
                    SQR (x3[0] - x1[0]) + SQR (x3[1] - x1[1]) + SQR (x3[2] - x1[2]));
        
                double max = length_edge1;
                if (length_edge2 > max) max = length_edge2;
                if (length_edge3 > max) max = length_edge3;
        
                // save min and max of length of longest edge
                if (bg->longest_side_max_m < max) bg->longest_side_max_m = max;
                if (bg->longest_side_min_m > max) bg->longest_side_min_m = max;
            }

            /*
              In principal the number of discretization nr_m is the extend of
              the geometry divided by the extend of the largest triangle. Whereby
              the extend of a triangle is defined as the lenght of its longest
              edge. Thus the largest triangle is the triangle with the longest edge.

              But if the hot spot, i.e., the multipacting/field emission zone is
              too small that the normal bounding box covers the whole hot spot, the
              expensive triangle-line intersection tests will be frequently called.
              In these cases, we have to use smaller bounding box size to speed up
              simulation. 

              Todo:
              The relation between bounding box size and simulation time step &
              geometry shape maybe need to be summarized and modeled in a more
              flexible manner and could be adjusted in input file.
            */
            bg->nr_m (0) = (int)floor (bg->len_m (0) / bg->longest_side_max_m * 8.0);
            bg->nr_m (1) = (int)floor (bg->len_m (1) / bg->longest_side_max_m * 8.0);
            bg->nr_m (2) = (int)floor (bg->len_m (2) / bg->longest_side_max_m * 8.0);

            bg->hr_m = bg->len_m / bg->nr_m;
            bg->voxelMesh_m.minExtend = bg->mincoords_m - 0.5 * bg->hr_m;
            bg->voxelMesh_m.maxExtend = bg->maxcoords_m + 0.5 * bg->hr_m;
            bg->voxelMesh_m.extend = bg->voxelMesh_m.maxExtend - bg->voxelMesh_m.minExtend;
            bg->nr_m += 1;

            bg->outside_point_m = bg->maxcoords_m + bg->hr_m;
            *gmsg << "* Geometry interval built done." << endl;
        }


#define surrounding_voxels( voxel_id, nx, ny ) {                        \
            (voxel_id) - (nx) * (ny) - (nx) - 1,                        \
                (voxel_id) - (nx) * (ny) - (nx),                        \
                (voxel_id) - (nx) * (ny) - (nx) + 1,                    \
                (voxel_id) - (nx) * (ny) - 1,                           \
                (voxel_id) - (nx) * (ny),                               \
                (voxel_id) - (nx) * (ny) + 1,                           \
                (voxel_id) - (nx) * (ny) + (nx) - 1,                    \
                (voxel_id) - (nx) * (ny) + (nx),                        \
                (voxel_id) - (nx) * (ny) + (nx) + 1,                    \
                (voxel_id) - (nx) - 1,                                  \
                (voxel_id) - (nx),                                      \
                (voxel_id) - (nx) + 1,                                  \
                (voxel_id) - 1,                                         \
                (voxel_id),                                             \
                (voxel_id) + 1,                                         \
                (voxel_id) + (nx) - 1,                                  \
                (voxel_id) + (nx),                                      \
                (voxel_id) + (nx) + 1,                                  \
                (voxel_id) + (nx) * (ny) - (nx) - 1,                    \
                (voxel_id) + (nx) * (ny) - (nx),                        \
                (voxel_id) + (nx) * (ny) - (nx) + 1,                    \
                (voxel_id) + (nx) * (ny) - 1,                           \
                (voxel_id) + (nx) * (ny),                               \
                (voxel_id) + (nx) * (ny) + 1,                           \
                (voxel_id) + (nx) * (ny) + (nx) - 1,                    \
                (voxel_id) + (nx) * (ny) + (nx),                        \
                (voxel_id) + (nx) * (ny) + (nx) + 1,                    \
                }
        /*
          add voxels to integer array
          sort this array
          get unique set 
         */

        static inline void computeLineVoxelization (
            BoundaryGeometry* bg,
            const Vector_t& P,
            const Vector_t& v,
            const int num_segments,
            std::unordered_set<int>& voxel_ids
            ) {
            Vector_t P1 = P;
            const Vector_t v_ = v / num_segments;
            for (int j = 0; j < num_segments; j++, P1 += v_) {
                const int voxel_id = bg->mapPoint2VoxelID (P + j*v_);
                assert (voxel_id > 0);
                const int idc[27] = 
                    surrounding_voxels (voxel_id, bg->nr_m[0], bg->nr_m[1]);
                int offset = -1;
                while (idc[++offset] <= 0);
                voxel_ids.insert (idc+offset, idc+27-offset);
            }
        }

        static inline void computeTriangleVoxelization (
            BoundaryGeometry* bg,
            const int triangle_id,
            std::unordered_set<int>& voxel_ids
            ) {
            const int num_segments = 16;

            /*
              Discretize the three central lines and the three triangle edges
              to get a more complete boundary index set.
            */
            const Vector_t V0 = bg->getPoint (triangle_id, 1);
            const Vector_t V1 = bg->getPoint (triangle_id, 2);
            const Vector_t V2 = bg->getPoint (triangle_id, 3);

            computeLineVoxelization (bg, V0, 0.5 * (V1 + V2) - V0, num_segments, voxel_ids);
            computeLineVoxelization (bg, V1, 0.5 * (V2 + V0) - V1, num_segments, voxel_ids);
            computeLineVoxelization (bg, V2, 0.5 * (V0 + V1) - V2, num_segments, voxel_ids);
            computeLineVoxelization (bg, V0, V1 - V0,              num_segments, voxel_ids);
            computeLineVoxelization (bg, V1, V2 - V1,              num_segments, voxel_ids);
            computeLineVoxelization (bg, V2, V0 - V2,              num_segments, voxel_ids);
        }

        static inline void computeTriangleVoxelizationViaBBox (
            BoundaryGeometry* bg,
            const int triangle_id,
            std::unordered_set<int>& voxel_ids
            ) {
            Vector_t v1 = bg->getPoint (triangle_id, 1);
            Vector_t v2 = bg->getPoint (triangle_id, 2);
            Vector_t v3 = bg->getPoint (triangle_id, 3);
            Vector_t bbox_min = {
                MIN3 (v1[0], v2[0], v3[0]),
                MIN3 (v1[1], v2[1], v3[1]),
                MIN3 (v1[2], v2[2], v3[2]) };
            Vector_t bbox_max = {
                MAX3 (v1[0], v2[0], v3[0]),
                MAX3 (v1[1], v2[1], v3[1]),
                MAX3 (v1[2], v2[2], v3[2]) };
            int i_min, j_min, k_min;
            int i_max, j_max, k_max;
            bg->mapPoint2VoxelIndices (bbox_min, i_min, j_min, k_min);
            bg->mapPoint2VoxelIndices (bbox_max, i_max, j_max, k_max);
            for (int i = i_min; i <= i_max; i++) {
                for (int j = j_min; j <= j_max; j++) {
                    for (int k = k_min; k <= k_max; k++) {
                        // test if voxel (i,j,k) has an intersection with triangle
                        if (bg->intersectTriangleVoxel (triangle_id, i, j, k) == INSIDE) {
                            int voxel_id = bg->mapVoxelIndices2ID (i, j, k);
                            voxel_ids.insert (voxel_id);
                        }
                    }
                }
            } 
        }

        static void computeMeshVoxelization (BoundaryGeometry* bg) {

            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                std::unordered_set<int> voxel_ids;
#if defined(FAST_VOXELIZATION)
                computeTriangleVoxelizationViaBBox (bg, triangle_id, voxel_ids);
#else
                computeTriangleVoxelization (bg, triangle_id, voxel_ids);
#endif
#if 0
                *gmsg << "* Triangle: " << triangle_id << " voxels:";
                for (auto it = voxel_ids.begin(); it != voxel_ids.end(); it++) {
                    *gmsg << " " << *it;
                }
                *gmsg << endl;
#endif
                // add voxeliziation of triangle to voxelization of mesh
                bg->boundaryVoxelIDs_m.insert (voxel_ids.begin(), voxel_ids.end());

                // 
                std::unordered_map< int, std::unordered_set<int> > map_of_voxel_ids_to_intersect_triangles;
                for (auto voxelIt = voxel_ids.begin ();
                     voxelIt != voxel_ids.end ();
                     voxelIt++) {
                    int voxel_id = *voxelIt;

                    auto it = map_of_voxel_ids_to_intersect_triangles.find (voxel_id);
                    if (it == map_of_voxel_ids_to_intersect_triangles.end ()) {
                        map_of_voxel_ids_to_intersect_triangles [voxel_id] =
                            std::unordered_set<int> (&triangle_id, &triangle_id+1);
                    } else
                        it->second.insert (triangle_id);
                }

                // add map for given triangle to map for mesh
                for (auto mapIt = map_of_voxel_ids_to_intersect_triangles.begin ();
                     mapIt != map_of_voxel_ids_to_intersect_triangles.end ();
                     mapIt++) {
                    auto it = bg->CubicLookupTable_m.find (mapIt->first);
                    if (it == bg->CubicLookupTable_m.end ()) {
                        bg->CubicLookupTable_m [mapIt->first] = mapIt->second;
                    } else {
                        it->second.insert (mapIt->second.begin(), mapIt->second.end());
                    }
                }
                if (triangle_id > 0 && (triangle_id % 1000) == 0)
                    *gmsg << "* Triangle ID: " << triangle_id << endl;
            } // for_each triangle
            *gmsg << "* Boundary index set built done." << endl;
            if(Ippl::myNode() == 0) {
                write_voxel_mesh (bg->boundaryVoxelIDs_m, bg->hr_m, bg->nr_m, bg->voxelMesh_m.minExtend);
            }
        }

/*
  
  Following combinations are possible:
              1,1 && 2,2   1,2 && 2,1   1,3 && 2,1
              1,1 && 2,3   1,2 && 2,3   1,3 && 2,2
              1,1 && 3,2   1,2 && 3,1   1,3 && 3,1
              1,1 && 3,3   1,2 && 3,3   1,3 && 3,2

             (2,1 && 1,2) (2,2 && 1,1) (2,3 && 1,1)
             (2,1 && 1,3) (2,2 && 1,3) (2,3 && 1,2)
              2,1 && 3,2   2,2 && 3,1   2,3 && 3,1
              2,1 && 3,3   2,2 && 3,3   2,3 && 3,2

             (3,1 && 1,2) (3,2 && 1,1) (3,3 && 1,1)
             (3,1 && 1,3) (3,2 && 1,3) (3,3 && 1,2)
             (3,1 && 2,2) (3,2 && 2,1) (3,3 && 2,1)
             (3,1 && 2,3) (3,2 && 2,3) (3,3 && 2,2)

  Note:
     Since we find vertices with lower enumeration first, we
     can ignore combinations in ()

                  2 2           2 3           3 2           3 3    
                   *             *             *             *     
                  /|\           /|\           /|\           /|\    
                 / | \         / | \         / | \         / | \   
                /  |  \       /  |  \       /  |  \       /  |  \  
               /   |   \     /   |   \     /   |   \     /   |   \ 
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 1   3   3   1 1   2   2   1 1   3   2   1 1   2      
diff:            (1,1)         (1,2)         (2,1)         (2,2)
change orient.:   yes           no            no            yes


                  2 1           2 3           3 1           3 3    
                   *             *             *             *     
                  /|\           /|\           /|\           /|\    
                 / | \         / | \         / | \         / | \   
                /  |  \       /  |  \       /  |  \       /  |  \  
               /   |   \     /   |   \     /   |   \     /   |   \ 
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 2   3   3   1 2   1   2   1 2   3   2   1 2   1      
diff:            (1,-1)        (1,1)         (2,-1)        (2,1)
change orient.:   no            yes           yes           no


                  2 1           2 2           3 1           3 2    
                   *             *             *             *     
                  /|\           /|\           /|\           /|\    
                 / | \         / | \         / | \         / | \   
                /  |  \       /  |  \       /  |  \       /  |  \  
               /   |   \     /   |   \     /   |   \     /   |   \ 
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 3   2   3   1 3   1   2   1 3   2   2   1 3   1      
diff:            (1,-2)        (1,-1)        (2,-2)        (2,-1)
change orient.:   yes           no            no            yes

                                              3 2           3 3    
                                               *             *     
                                              /|\           /|\    
                                             / | \         / | \   
                                            /  |  \       /  |  \  
                                           /   |   \     /   |   \ 
                                          *----*----*   *----*----*
                                          1   2 1   3   1   2 1   2
diff:                                        (1,1)         (1,2)
change orient.:                               yes           no

                                              3 1           3 3
                                               *             *     
                                              /|\           /|\    
                                             / | \         / | \   
                                            /  |  \       /  |  \  
                                           /   |   \     /   |   \ 
                                          *----*----*   *----*----*
                                          1   2 2   3   1   2 2   1
diff:                                        (1,-1)        (1,1)
change orient.:                               no            yes

                                              3 1           3 2
                                               *             *     
                                              /|\           /|\    
                                             / | \         / | \   
                                            /  |  \       /  |  \  
                                           /   |   \     /   |   \ 
                                          *----*----*   *----*----*
                                          1   2 3   2   1   2 3   1
diff:                                        (1,-2)        (1,-1)
change orient.:                               yes           no


Change orientation if diff is:
(1,1) || (1,-2) || (2,2) || (2,-1) || (2,-1)

*/

        static void orientTriangles (BoundaryGeometry* bg, int ref_id, int triangle_id) {

            // find pts of common edge
            int ic[2];
            int id[2];
            int n = 0;

            for (int i = 1; i <= 3; i++) {
                for (int j = 1; j <= 3; j++) {
                    if (bg->PointID (triangle_id, j) == bg->PointID (ref_id, i)) {
                        id[n] = j;
                        ic[n] = i;
                        n++;
                        if (n == 2) goto edge_found;
                    }
                }
            }
            assert (n == 2);
        edge_found:
            int diff = id[1] - id[0];
            if ((((ic[1] - ic[0]) == 1) && ((diff == 1) || (diff == -2))) || 
                (((ic[1] - ic[0]) == 2) && ((diff == -1) || (diff == 2)))) {
                bg->PointID (triangle_id, id[0]) = bg->PointID (ref_id, ic[1]);
                bg->PointID (triangle_id, id[1]) = bg->PointID (ref_id, ic[0]);
            }
            bg->isOriented_m [triangle_id] = true;
            std::unordered_set<int> neighbors = bg->triangleNeighbors_m[triangle_id];

            const auto endIt = neighbors.end ();
            for (auto triangleIt = neighbors.begin(); triangleIt != endIt; triangleIt++) {
                if (!bg->isOriented_m [*triangleIt])
                    orientTriangles (bg, triangle_id, *triangleIt);
            }
        }

        /*
          Determine if a point x is outside or inside the geometry or just on
          the boundary. Return true if point is inside geometry or on the
          boundary, false otherwise

          The basic idea here is:
          If a line segment from the point to test to a random point outside
          the geometry has has an even number of intersections with the
          boundary, the point is outside the geometry.

          Note:
          If the point is on the boundary, the number of intersections is 1.
          Points on the boundary are handled as inside.

          A random selection of the reference point outside the boundary avoids
          some specific issues, like line parallel to boundary.
         */
        static inline bool isInside (BoundaryGeometry* bg, const Vector_t x) {
            IpplTimings::startTimer (bg->Tinward_m);

            Vector_t y = Vector_t (
                //x[0],
                bg->maxcoords_m[0] * (1.1 + gsl_rng_uniform(bg->randGen_m)),
                bg->maxcoords_m[1] * (1.1 + gsl_rng_uniform(bg->randGen_m)),
                bg->maxcoords_m[2] * (1.1 + gsl_rng_uniform(bg->randGen_m)));

            std::vector<Vector_t> intersection_points;
            int num_intersections = 0;

            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                Vector_t result;
                if (bg->intersectLineTriangle (SEGMENT, x, y, triangle_id, result)) {
                    intersection_points.push_back (result);
                    num_intersections++;
                }
            }
            IpplTimings::stopTimer (bg->Tinward_m);
            return ((intersection_points.size () % 2) == 1);
        }


        static void computeTriangleNeighbors (
            BoundaryGeometry* bg
            ) {
            std::unordered_set<int> adjacencies_to_pt  [bg->num_points_m];

            // for each triangles find adjacent triangles for each vertex
            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                for (int j = 1; j <= 3; j++) {
                    int pt_id = bg->PointID (triangle_id, j);
                    assert (pt_id < bg->num_points_m);
                    adjacencies_to_pt [pt_id].insert (triangle_id);
                }
            }

            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                std::unordered_set<int>  to_A = adjacencies_to_pt [bg->PointID (triangle_id, 1)];
                std::unordered_set<int>  to_B = adjacencies_to_pt [bg->PointID (triangle_id, 2)];
                std::unordered_set<int>  to_C = adjacencies_to_pt [bg->PointID (triangle_id, 3)];

                std::unordered_set<int>  intersect;
                std::set_intersection (
                    to_A.begin(), to_A.end(),
                    to_B.begin(), to_B.end(),
                    std::inserter(intersect,intersect.begin()));
                std::set_intersection(
                    to_B.begin(), to_B.end(),
                    to_C.begin(), to_C.end(),
                    std::inserter(intersect,intersect.begin()));
                std::set_intersection(
                    to_C.begin(), to_C.end(),
                    to_A.begin(), to_A.end(),
                    std::inserter(intersect,intersect.begin()));

                bg->triangleNeighbors_m.insert (std::pair <int, std::unordered_set<int>> (triangle_id, intersect));
                                
            }
        }

        static inline Vector_t normalVector (
            BoundaryGeometry* const bg,
            const int triangle_id
            ) {
            const Vector_t A = bg->getPoint (triangle_id, 1);
            const Vector_t B = bg->getPoint (triangle_id, 2);
            const Vector_t C = bg->getPoint (triangle_id, 3);

            const Vector_t N = cross (B - A, C - A);
            const double magnitude = sqrt (SQR (N (0)) + SQR (N (1)) + SQR (N (2)));
            assert (fcmp (magnitude, 0.0, 10) > 0); // in case we have degenerted triangles
            return N / magnitude;
        }

        static bool hasInwardPointingNormal (
            BoundaryGeometry* const bg,
            const int triangle_id
            ) {
            const Vector_t A = bg->getPoint (triangle_id, 1);
            const Vector_t B = bg->getPoint (triangle_id, 2);
            const Vector_t C = bg->getPoint (triangle_id, 3);
            
            // choose a point P close to the center of the triangle
            const Vector_t P = (A+B+C)/3 + bg->TriNormal_m[triangle_id] * 0.1 * bg->longest_side_min_m;

            /*
              The triangle normal points inward, if P is
              - outside the geometry and the dot product is negativ
              - or inside the geometry and the dot product is positiv

              Remember:
                The dot product is positiv only if both vectors are
                pointing in the same direction.
            */
            const bool is_inside = isInside (bg, P);
            const double dotPA_N = dot (P - A, bg->TriNormal_m[0]);
            return (is_inside && dotPA_N >= 0) || (!is_inside && dotPA_N < 0);
        }

        /*
          Recursively get inward-pointing normal of all surface triangles.
          
          The basic idea is as follow:
          -  get the inward-pointing normal of the first triangle by determine
             whether a nearby point is inside or outside the boundary geometry.
             (using ray-triangle intersection and even/odd intersection number).
          -  Then use a recursion method to switch the vertex order of adjacent
             triangles. The inward normal is stored in TriNormal_m.
        */
        static void makeTriangleNormalInwardPointing (BoundaryGeometry* bg) {
            bg->isOriented_m = new bool[bg->num_triangles_m];
            memset (bg->isOriented_m, 0, sizeof (bg->isOriented_m[0])*bg->num_triangles_m);
            bg->TriNormal_m.reserve (bg->num_triangles_m);

            bg->TriNormal_m[0] = normalVector (bg, 0);
            if (!hasInwardPointingNormal (bg, 0)) {
                bg->TriNormal_m[0] = -bg->TriNormal_m[0];
                int id = bg->PointID (0, 2);
                bg->PointID (0, 2) = bg->PointID (0, 3);
                bg->PointID (0, 3) = id;
            }
            bg->isOriented_m [0] = true;

            /*
              compute edge-neightbors for each triangle
              then orient all triangles, starting with the neighbors of the first
            */
            computeTriangleNeighbors (bg);

            const auto neighbors = bg->triangleNeighbors_m[0];
            const auto endIt = neighbors.end ();
            for (auto triangleIt = neighbors.begin(); triangleIt != endIt; triangleIt++) {
                if (!bg->isOriented_m [*triangleIt])
                    orientTriangles (bg, 0, *triangleIt);
            }

            // compute inward-normals
            for (int triangle_id = 1; triangle_id < bg->num_triangles_m; triangle_id++) {
                bg->TriNormal_m.push_back (normalVector (bg, triangle_id));
            }
            *gmsg << "* Triangle Normal built done." << endl;
        }

        // Calculate the area of triangle given by id.
        static inline double computeArea (BoundaryGeometry* bg, int id) {
            Vector_t AB = bg->getPoint (id, 2) - bg->getPoint (id, 1);
            Vector_t AC = bg->getPoint (id, 3) - bg->getPoint (id, 1);
            return(0.5 * sqrt (dot (AB, AB) * dot (AC, AC) - dot (AB, AC) * dot (AB, AC)));
        }

        /*
          We define some tags in namespace BGphysics for each surface triangle to
          identify the physical reactions for each triangle when amplitude of
          electrostatic field exceeds some threshold or particles incident the surface.
        */
        static void setBGphysicstag (BoundaryGeometry* bg) {
            for (int i = 0; i < bg->num_triangles_m; i++) {
                bg->TriBGphysicstag_m.push_back (
                    BGphysics::Absorption
                    | BGphysics::FNEmission
                    | BGphysics::SecondaryEmission);
            }
        }

    };

    h5_int64_t rc;

    *gmsg << "* Initializing Boundary Geometry..." << endl;
    IpplTimings::startTimer (TPreProc_m);

    apert_m = Attributes::getRealArray(itsAttr[APERTURE]);
 
    if (hasApperture()) {
        *gmsg << "* Found additional aperture." << endl;
        for (unsigned int i=0; i<apert_m.size(); i=i+3)
            *gmsg << "* zmin = " << apert_m[i]
                  << " zmax = " << apert_m[i+1]
                  << " r= " << apert_m[i+2] << endl;
    }

    *gmsg << "* Filename: " << h5FileName_m.c_str() << endl;

    double xscale = Attributes::getReal(itsAttr[XSCALE]); 
    double yscale = Attributes::getReal(itsAttr[YSCALE]); 
    double zscale = Attributes::getReal(itsAttr[ZSCALE]); 

    double xyzscale = Attributes::getReal(itsAttr[XYZSCALE]); 

    *gmsg << "* X-scale all points of geometry by " << xscale << endl;
    *gmsg << "* Y-scale all points of geometry by " << yscale << endl;
    *gmsg << "* Z-scale all points of geometry by " << zscale << endl;
    *gmsg << "* Scale all points of geometry by " << xyzscale << endl;

    rc = H5SetErrorHandler (H5AbortErrorhandler);
    if (rc != H5_SUCCESS)
        ERRORMSG ("H5 rc = " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    H5SetVerbosityLevel (1);
    h5_file_t* f = H5OpenFile (h5FileName_m.c_str (), H5_O_RDONLY, Ippl::getComm());
    h5t_mesh_t* m = NULL;
    H5FedOpenTriangleMesh (f, "0", &m);
    H5FedSetLevel (m, 0);

    num_triangles_m = H5FedGetNumElementsTotal (m);
    allbfaces_m = new int[num_triangles_m * 4];

    // iterate over all co-dim 0 entities, i.e. elements
    h5_loc_id_t local_id;
    int i = 0;
    h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
    while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
        h5_loc_id_t local_vids[4];
        H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
        PointID (i, 0) = 0;
        PointID (i, 1) = local_vids[0];
        PointID (i, 2) = local_vids[1];
        PointID (i, 3) = local_vids[2];
        i++;
    }
    H5FedEndTraverseEntities (iter);

    // loop over all vertices
    num_points_m = H5FedGetNumVerticesTotal (m);
    double* point_coords = new double[3 * num_points_m];
    for (i = 0; i < num_points_m; i++) {
        h5_float64_t P[3];
        H5FedGetVertexCoordsByIndex (m, i, P);
        point_coords[i * 3]     = P[0] * xyzscale * xscale;
        point_coords[i * 3 + 1] = P[1] * xyzscale * yscale;
        point_coords[i * 3 + 2] = P[2] * xyzscale * zscale;
    }
    H5FedCloseMesh (m);
    H5CloseFile (f);

    double zshift = (double)(Attributes::getReal (itsAttr[ZSHIFT]));

    for (int i = 0; i < num_points_m; i++) {
        geo3Dcoords_m.push_back (
            Vector_t (
                point_coords[3 * i],
                point_coords[3 * i + 1],
                point_coords[3 * i + 2] + zshift));
    }
    delete point_coords;
    *gmsg << "* Vertex built done." << endl;

    Local::computeGeometryInterval (this);

    Local::computeMeshVoxelization (this);
    Local::makeTriangleNormalInwardPointing (this);
    Local::setBGphysicstag (this);

    Tribarycent_m = new Vector_t[num_triangles_m];
    TriPrPartloss_m = new double[num_triangles_m];
    TriFEPartloss_m = new double[num_triangles_m];
    TriSePartloss_m = new double[num_triangles_m];
    for (int i = 0; i < num_triangles_m; i++) {
        Tribarycent_m[i] = (getPoint (i, 1) + getPoint (i, 2) + getPoint (i, 3)) / 3.0;
        Triarea_m.push_back (Local::computeArea (this, i));

        TriPrPartloss_m[i] = 0.0;
        TriFEPartloss_m[i] = 0.0;
        TriSePartloss_m[i] = 0.0;
    }
    *gmsg << "* Triangle barycent built done." << endl;

    *gmsg << *this << endl;
    IpplTimings::stopTimer (TPreProc_m);
    exit(42);
}


/*
  Compute intersection between line-segment given by the endpoints P0 and P1
  and boundary geometry.

  result:
    0   no intersection
    1   line-segment is on boundary
    2   unique intersection, but both Pts are outside
    3   unique intersection in segment

    WARNING:
    This is *not* a general porpose method! If neither P0 nor P1 is in the
    voxelized boundary geometry, no intersection will be found!
*/
int
BoundaryGeometry::intersectLineSegmentBoundary (
    const Vector_t P0,                  // [in] starting point of ray
    const Vector_t P1,                  // [in] end point of ray
    Vector_t& intersection_pt,          // [out] intersection with boundary
    int& triangle_id                    // [out] triangle the line segment intersects with
    ) {
    triangle_id = -1;
    /*
      Test if P0 or P1 is inside the voxelized boundary
      If not, we are done ...
    */
    int voxel_id;
    if (boundaryVoxelIDs_m.find (mapPoint2VoxelID (P0)) != boundaryVoxelIDs_m.end ()) {
        // particle is in geometry at timestep n
        voxel_id = mapPoint2VoxelID (P0);
    } else if (boundaryVoxelIDs_m.find (mapPoint2VoxelID (P1)) != boundaryVoxelIDs_m.end ()) {
        // particle is in geometry at timestep n+1
        voxel_id = mapPoint2VoxelID (P1);
    } else {
        return 0;
    }

    /*
      Build an array containing the IDs of 27(3*3*3) voxels. The ID
      of the voxel containing the particle, is the center of these
      voxels. Just for the situation that even the line segment
      crosses more than one voxel.
    */
    const int idc[27] = surrounding_voxels (voxel_id, nr_m[0], nr_m[1]);
    
    /* Test all the 27 voxels to find if the line segment has
       intersection with the triangles in those voxels. */
    const Vector_t v = P1 - P0;
    int intersect_result = 0;
    for (int k = 0; k < 27; k++) {
        const auto triangles_overlaping_with_voxel = CubicLookupTable_m.find (idc[k]);
        if (triangles_overlaping_with_voxel == CubicLookupTable_m.end ())
            continue;                   // not in voxelization of boundary
        
        // for each triangle in this voxel
        for (auto it = triangles_overlaping_with_voxel->second.begin ();
             it != triangles_overlaping_with_voxel->second.end ();
             it++) {
            if (dot (v, TriNormal_m[*it]) >= 0.0) 
                continue;               // particle moves away from triangle
            
            // particle moves towards triangle
            intersect_result = intersectLineTriangle (
                LINE,
                P0, P1,
                *it,
                intersection_pt);
            switch (intersect_result) {
            case -1:                    // triangle is degenerated
                assert (intersect_result != -1);
                exit (42);              // terminate even if NDEBUG is set
            case 0:                     // no intersection
            case 4:                     // both points are inside
                intersect_result = 0;
                break;              
            case 1:                     // line and triangle are in same plane
            case 2:                     // both points are outside
            case 3:                     // unique intersection in segment
                *gmsg << "* Intersection test returned: " << intersect_result << endl;
                triangle_id = (*it);
                goto done;
            };
        }                               // end for all triangles
    }                                   // end for all voxels
done:
    return intersect_result;
}

/**
   Determine whether a particle with position @param r, momenta @param v,
   and time step @param dt will hit the boundary.

   Basic algorithms are as follows:
   1) Determine if the particle is near the boundary by checking r is in
      boundary bounding box index set.
      if r is in, then do the following checking step, else return -1 to
      indicate that particle is far from boundary and to be integrated.
   2) Traversal all the triangles in the bounding cubic box which the
       particle is in, as well as triangles in the adjacent 26 bounding
       cubic boxes
       if the momenta has oppsite direction with those triangles' normals,
       then check if the particle has intersection with those triangles. If
       intersection exsists, then return 0.
 */
int
BoundaryGeometry::PartInside (
    const Vector_t r,                   // [in] particle position
    const Vector_t v,                   // [in] momentum
    const double dt,                    // [in]
    const int Parttype,                 // [in] type of particle
    const double Qloss,                 // [in]
    Vector_t& intecoords,               // [out] intersection with boundary
    int& triangle_id                    // [out] appropriate triangle 
    ) {
    int ret = -1;
    if (v == (Vector_t)0)               // nothing to do if momenta == 0
        return ret;

    const double p_sq = dot (v, v);
    const double betaP = 1.0 / sqrt (1.0 + p_sq);

    const Vector_t P0 = r;              //particle position in timestep n;
    const Vector_t P1 = r + (Physics::c * betaP * v * dt); //particle position in tstep n+1

    IpplTimings::startTimer (TPInside_m);
    Vector_t intersection_pt = 0.0;
    intersectLineSegmentBoundary (P0, P1, intersection_pt, triangle_id);
    if (triangle_id >= 0) {
        intecoords = intersection_pt;
        if (Parttype == 0)
            TriPrPartloss_m[triangle_id] += Qloss;
        else if (Parttype == 1)
            TriFEPartloss_m[triangle_id] += Qloss;
        else
            TriSePartloss_m[triangle_id] += Qloss;
        ret = 0;
    }
    IpplTimings::stopTimer (TPInside_m);
    return ret;
}

void
BoundaryGeometry::writeGeomToVtk (string fn) {
    std::ofstream of;
    of.open (fn.c_str ());
    assert (of.is_open ());
    of.precision (6);
    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated using DataSink::writeGeoToVtk" << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
    of << "POINTS " << num_points_m << " float" << std::endl;
    for (int i = 0; i < num_points_m; i++)
        of << geo3Dcoords_m[i](0) << " "
	   << geo3Dcoords_m[i](1) << " "
	   << geo3Dcoords_m[i](2) << std::endl;
    of << std::endl;

    of << "CELLS "
       << num_triangles_m << " "
       << 4 * num_triangles_m << std::endl;
    for (int i = 0; i < num_triangles_m; i++)
        of << "3 "
	   << PointID (i, 1) << " "
	   << PointID (i, 2) << " "
	   << PointID (i, 3) << std::endl;
    of << "CELL_TYPES " << num_triangles_m << std::endl;
    for (int i = 0; i < num_triangles_m; i++)
	of << "5" << std::endl;
    of << "CELL_DATA " << num_triangles_m << std::endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
    of << "LOOKUP_TABLE " << "default" << std::endl;
    for (int i = 0; i < num_triangles_m; i++)
	of << (float)(i) << std::endl;
    of << std::endl;
}

Inform&
BoundaryGeometry::printInfo (Inform& os) const {
    os << "* *************Boundary Geometry Info*********************************************** " << endl;
    os << "* GEOMETRY                   " << getOpalName () << '\n'
       << "* FGEOM                      " << Attributes::getString (itsAttr[FGEOM]) << '\n'
       << "* TOPO                       " << Attributes::getString (itsAttr[TOPO]) << '\n'
       << "* LENGTH                     " << Attributes::getReal (itsAttr[LENGTH]) << '\n'
       << "* S                          " << Attributes::getReal (itsAttr[S]) << '\n'
       << "* A                          " << Attributes::getReal (itsAttr[A]) << '\n'
       << "* B                          " << Attributes::getReal (itsAttr[B]) << '\n';
    if (getTopology () == string ("BOXCORNER")) {
        os << "* C                          " << Attributes::getReal (itsAttr[C]) << '\n'
           << "* L1                         " << Attributes::getReal (itsAttr[L1]) << '\n'
           << "* L1                         " << Attributes::getReal (itsAttr[L2]) << '\n';
    }
    os << "* Total triangle num         " << num_triangles_m << '\n'
       << "* Total points num           " << num_points_m << '\n'
       << "* Triangle side(m)   Max=    " << longest_side_max_m << '\n'
       << "*                    Min=    " << longest_side_min_m << '\n'
       << "* Geometry bounds(m) Max=    " << maxcoords_m << '\n'
       << "*                    Min=    " << mincoords_m << '\n'
       << "* Geometry length(m)         " << len_m << '\n'
       << "* Boundary box grid num      " << nr_m << '\n'
       << "* Boundary box size(m)       " << hr_m << '\n'
       << "* Size of boundary index set " << boundaryVoxelIDs_m.size () << '\n'
       << "* Number of all boxes        " << nr_m (0) * nr_m (1) * nr_m (2) << '\n'
        << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}

/*
   ____  _               _          
  |  _ \| |__  _   _ ___(_) ___ ___ 
  | |_) | '_ \| | | / __| |/ __/ __|
  |  __/| | | | |_| \__ \ | (__\__ \
  |_|   |_| |_|\__, |___/_|\___|___/
                |___/          

  start here ...
*/

/**
   Determine physical behaviour when particle hits the boundary triangle,
   non secondary emission version.
 */
int BoundaryGeometry::doBGphysics (
    const Vector_t& intecoords,
    const int& triId
    ) {
    short BGtag = TriBGphysicstag_m[triId];
    if ((BGtag & BGphysics::Nop) == BGphysics::Nop) {
        return -1;
    } else if (((BGtag & BGphysics::Absorption) == BGphysics::Absorption)
               && ((BGtag & BGphysics::FNEmission) != BGphysics::FNEmission)) {
        return 0;
    } else {
        return 1;
    }
}

/**
   Determine physical behaviour when particle hits the boundary triangle,
   call Furman-Pivi's secondary emission model.
 */
int BoundaryGeometry::doBGphysics (
    const Vector_t& intecoords,
    const int& triId,
    const double& incQ,
    const Vector_t& incMomentum,
    PartBunch* itsBunch,
    double& seyNum
    ) {
    Inform msg ("BGphyDebug");

    const double p_sq = dot (incMomentum, incMomentum);    
    const double incEnergy = Physics::m_e * (sqrt (1.0 + p_sq) - 1.0) * 1.0e9;   // energy in eV

    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if ((BGtag & BGphysics::Nop) == BGphysics::Nop) {
        ret = - 1;
    } else if ((BGtag & BGphysics::Absorption) == BGphysics::Absorption &&
               (BGtag & BGphysics::FNEmission) != BGphysics::FNEmission &&
               (BGtag & BGphysics::SecondaryEmission) != BGphysics::SecondaryEmission) {
        ret = 0;
    } else {
        // Secondary Emission;
        ret = 1;
        int se_Num = 0;
        int seType = 0;
        if ((BGtag & BGphysics::SecondaryEmission) == BGphysics::SecondaryEmission) {
            double cosTheta = - dot (incMomentum, TriNormal_m[triId]) /
                              sqrt (dot (incMomentum, incMomentum));
            if (cosTheta < 0) {
                //cosTheta should be positive
                std::cout << "cosTheta = " << cosTheta
                          << " intecoords " << intecoords (0) << " " << intecoords (1)
                          << " " << intecoords (2) << " "
                          << std::endl;
                std::cout << "incident momentum=("
                          << incMomentum (0) << "," << incMomentum (1) << "," << incMomentum (2) << ")"
                          << " triNormal=("
                          << TriNormal_m[triId](0) << ","
                          << TriNormal_m[triId](1) << "," << TriNormal_m[triId](2) << ") "
                          << std::endl;
            }
            assert(cosTheta>=0);
            int idx = 0;
            if (intecoords != Point (triId, 1)) {
                idx = 1; // intersection is not the 1st vertex
            } else {
                idx = 2; // intersection is the 1st vertex
            }
            sec_phys_m.nSec (incEnergy,
                             cosTheta,
                             seBoundaryMatType_m,
                             se_Num,
                             seType,
                             incQ,
                             TriNormal_m[triId],
                             intecoords,
                             Point (triId, idx),
                             itsBunch,
                             seyNum,
                             ppVw_m,
                             vVThermal_m,
                             nEmissionMode_m);
        }
    }
    return ret;
}

/**
   Determine physical behaviour when particle hits the boundary triangle,
   call Vaughan's secondary emission model.
 */
int BoundaryGeometry::doBGphysics (
    const Vector_t& intecoords,
    const int& triId,
    const double& incQ,
    const Vector_t& incMomentum,
    PartBunch* itsBunch,
    double& seyNum,
    const int& para_null
    ) {
    const double p_sq = dot (incMomentum, incMomentum);    
    const double incEnergy = Physics::m_e * (sqrt (1.0 + p_sq) - 1.0) * 1.0e9;   // energy in eV

    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if ((BGtag & BGphysics::Nop) == BGphysics::Nop) {
        ret = - 1;
    } else if (((BGtag & BGphysics::Absorption) == BGphysics::Absorption) &&
               ((BGtag & BGphysics::FNEmission) != BGphysics::FNEmission) &&
               ((BGtag & BGphysics::SecondaryEmission) != BGphysics::SecondaryEmission)) {
        ret = 0;
    } else {
        // Secondary Emission;
        int se_Num = 0;
        int seType = 0;
        if ((BGtag & BGphysics::SecondaryEmission) == BGphysics::SecondaryEmission) {
            double cosTheta = - dot (incMomentum, TriNormal_m[triId]) /
                              sqrt (dot (incMomentum, incMomentum));
            //cosTheta should be positive
            if (cosTheta < 0) {
                std::cout << "cosTheta = " << cosTheta << std::endl;
                INFOMSG ("incident momentum=" << incMomentum
                                              << " triNormal=" << TriNormal_m[triId]
                                              << " dot=" << dot (incMomentum, TriNormal_m[triId])
                                              << endl);
            }
            //assert(cosTheta>=0);
            int idx = 0;
            if (intecoords != Point (triId, 1)) {
                // intersection is not the 1st vertex
                idx = 1;
            } else {
                // intersection is the 1st vertex
                idx = 2;
            }
            sec_phys_m.nSec (incEnergy,
                             cosTheta,
                             se_Num,
                             seType,
                             incQ,
                             TriNormal_m[triId],
                             intecoords,
                             Point (triId, idx),
                             itsBunch,
                             seyNum,
                             ppVw_m,
                             vSeyZero_m,
                             vEzero_m,
                             vSeyMax_m,
                             vEmax_m,
                             vKenergy_m,
                             vKtheta_m,
                             vVThermal_m,
                             nEmissionMode_m);
        }
    }
    return ret;
}

/**
   Here we call field emission model.
 */
/// \returns size_t
///     - number of emitted electrons at the surface
size_t BoundaryGeometry::doFNemission (
    OpalBeamline& itsOpalBeamline,
    PartBunch* itsBunch,
    const double t
    ) {
    // Self-field is not considered at moment. Only 1D Child-Langmuir law is
    // implemented for space charge limited current density.
    const double fa = parameterFNA_m / workFunction_m * fieldEnhancement_m * fieldEnhancement_m;
    /*  int node_num = Ippl::getNodes();

       size_t *count = new size_t [node_num];
       // itsBunch->getLocalNum();
       for(int i = 0; i < node_num; i++) {

       count[i] = 0;

       }*/
    size_t Nstp = 0;
    for (int i = 0; i < num_triangles_m; i++) {
        if ((TriBGphysicstag_m[i] & BGphysics::FNEmission) == BGphysics::FNEmission) {
            Vector_t E (0.0), B (0.0);
            Vector_t centroid (0.0);
            itsOpalBeamline.getFieldAt (Tribarycent_m[i], centroid, t, E, B);
            double Enormal = dot (TriNormal_m[i], E);
            /* Enormal should be negative as E field direction should be
               opposite to inward normal of surface */
            if (Enormal < fieldFNthreshold_m) {
                std::vector<Vector_t> vertex;
                vertex.push_back (Point (i, 1));
                vertex.push_back (Point (i, 2));
                vertex.push_back (Point (i, 3));
                PriEmissionPhysics::Fieldemission (itsBunch, fa, Enormal,
                                                   parameterFNB_m,
                                                   workFunction_m,
                                                   parameterFNVYZe_m,
                                                   parameterFNVYSe_m,
                                                   parameterFNY_m,
                                                   fieldEnhancement_m,
                                                   maxFNemission_m,
                                                   Triarea_m[i],
                                                   vertex,
                                                   TriNormal_m[i],
                                                   Nstp);
            }
        }
    }
    *gmsg << "* Emit " << Nstp << " field emission particles at the surfaces" << endl;
    return Nstp;
}

/**
   Initialize some darkcurrent particles near the surface with inward
      momenta.
 */
void BoundaryGeometry::createParticlesOnSurface (
    size_t n,
    double darkinward,
    OpalBeamline& itsOpalBeamline,
    PartBunch& itsBunch
    ) {
    int tag = 1002;
    int Parent = 0;
    if (Ippl::myNode () == 0) {
        for (size_t i = 0; i < n; i++) {
            short BGtag = BGphysics::Absorption;
            int k = 0;
            Vector_t E (0.0), B (0.0);
            while (((BGtag & BGphysics::Absorption) == BGphysics::Absorption &&
                    (BGtag & BGphysics::FNEmission) != BGphysics::FNEmission &&
                    (BGtag & BGphysics::SecondaryEmission) != BGphysics::SecondaryEmission)
                   ||
                   (fabs (E (0)) < eInitThreshold_m &&
                    fabs (E (1)) < eInitThreshold_m &&
                    fabs (E (2)) < eInitThreshold_m)) {
                E = Vector_t (0.0);
                B = Vector_t (0.0);
                int tmp = (int)(IpplRandom () * num_triangles_m);
                BGtag = TriBGphysicstag_m[tmp];
                k = tmp;
                Vector_t centroid (0.0);
                itsOpalBeamline.getFieldAt (Tribarycent_m[k] + darkinward * TriNormal_m[k],
                                            centroid, itsBunch.getdT (), E, B);
            }
            partsr_m.push_back (Tribarycent_m[k] + darkinward * TriNormal_m[k]);

        }
        Message* mess = new Message ();
        putMessage (*mess, partsr_m.size ());
        for (std::vector<Vector_t>::iterator myIt = partsr_m.begin (); myIt != partsr_m.end (); myIt++) {
            putMessage (*mess, *myIt);

        }
        Ippl::Comm->broadcast_all (mess, tag);
    } else {
        // receive particle position message
        size_t nData = 0;
        Message* mess = Ippl::Comm->receive_block (Parent, tag);
        getMessage (*mess, nData);
        for (size_t i = 0; i < nData; i++) {
            Vector_t tmp = Vector_t (0.0);
            getMessage (*mess, tmp);
            partsr_m.push_back (tmp);
        }

    }
}

/**
   Initialize primary particles near the surface with inward momenta.
 */
void BoundaryGeometry::createPriPart (
    size_t n,
    double darkinward,
    OpalBeamline& itsOpalBeamline,
    PartBunch* itsBunch
    ) {
    if (Options::ppdebug) {
        int tag = 1001;
        int Parent = 0;
        if (Ippl::myNode () == 0) {
            /* limit the initial particle in the center of the lower
               parallel plate. There is a distance of 0.01*length in
               x direction as margin. */
            double x_low = mincoords_m (0) + 0.5 * len_m (0) - 0.49 * len_m (0);

            /* limit the initial particle in the center of the upper
               parallel
               plate. There is a distance of 0.01*length in x direction as
               margin. */
            double x_up = mincoords_m (0) + 0.5 * len_m (0) + 0.49 * len_m (0);

            /* limit the initial particle in the center of the lower
               parallel
               plate. There is a distance of 0.01*length in y direction as
               margin. */
            double y_low = mincoords_m (1) + 0.5 * len_m (1) - 0.49 * len_m (1);

            /* limit the initial particle in the center of the upper
               parallel
               plate. There is a distance of 0.01*length in y direction as
               margin. */
            double y_up = mincoords_m (1) + 0.5 * len_m (1) + 0.49 * len_m (1);

            for (size_t i = 0; i < n / 2; i++) {
                double zCoord = maxcoords_m (2);
                double xCoord = maxcoords_m (0);
                double yCoord = maxcoords_m (1);
                while (zCoord > 0.000001 ||
                       zCoord < - 0.000001 ||
                       xCoord > x_up ||
                       xCoord < x_low ||
                       yCoord > y_up ||
                       yCoord < y_low) {

                    int k = (int)(IpplRandom () * num_triangles_m);
                    zCoord = Tribarycent_m[k](2);
                    xCoord = Tribarycent_m[k](0);
                    yCoord = Tribarycent_m[k](1);
                    if (Tribarycent_m[k](2) < 0.000001 &&
                        Tribarycent_m[k](2) > - 0.000001 &&
                        Tribarycent_m[k](0) < x_up &&
                        Tribarycent_m[k](0) > x_low &&
                        Tribarycent_m[k](1) < y_up &&
                        Tribarycent_m[k](1) > y_low) {
                        partsr_m.push_back (Tribarycent_m[k] + darkinward * TriNormal_m[k]);
                        partsp_m.push_back (TriNormal_m[k]);
                    }
                }
            }
            for (size_t i = 0; i < n / 2; i++) {
                double zCoord = maxcoords_m (2);
                double xCoord = maxcoords_m (0);
                double yCoord = maxcoords_m (1);
                while (zCoord > (maxcoords_m (2) + 0.000001) ||
                       (zCoord < (maxcoords_m (2) - 0.00000)) ||
                       xCoord > x_up ||
                       xCoord < x_low ||
                                yCoord > y_up ||
                       yCoord < y_low) {
                    int k = (int)(IpplRandom () * num_triangles_m);
                    zCoord = Tribarycent_m[k](2);
                    xCoord = Tribarycent_m[k](0);
                    yCoord = Tribarycent_m[k](1);
                    if ((Tribarycent_m[k](2) < maxcoords_m (2) + 0.000001) &&
                        (Tribarycent_m[k](2) > maxcoords_m (2) - 0.000001) &&
                        Tribarycent_m[k](0) < x_up &&
                        Tribarycent_m[k](0) > x_low &&
                        Tribarycent_m[k](1) < y_up &&
                        Tribarycent_m[k](1) > y_low) {
                        partsr_m.push_back (Tribarycent_m[k] + darkinward * TriNormal_m[k]);
                        partsp_m.push_back (TriNormal_m[k]);
                    }
                }
            }

            Message* mess = new Message ();
            putMessage (*mess, partsr_m.size ());
            for (std::vector<Vector_t>::iterator myIt = partsr_m.begin (),
                     myItp = partsp_m.begin ();
                 myIt != partsr_m.end ();
                 myIt++, myItp++) {
                putMessage (*mess, *myIt);
                putMessage (*mess, *myItp);
            }
            Ippl::Comm->broadcast_all (mess, tag);
        } else {
            // receive particle position message
            size_t nData = 0;
            Message* mess = Ippl::Comm->receive_block (Parent, tag);
            getMessage (*mess, nData);
            for (size_t i = 0; i < nData; i++) {
                Vector_t tmpr = Vector_t (0.0);
                Vector_t tmpp = Vector_t (0.0);
                getMessage (*mess, tmpr);
                getMessage (*mess, tmpp);
                partsr_m.push_back (tmpr);
                partsp_m.push_back (tmpp);
            }
        }
    } else {
        int tag = 1001;
        int Parent = 0;
        if (Ippl::myNode () == 0) {
            for (size_t i = 0; i < n; i++) {
                short BGtag = BGphysics::Absorption;
                int k = 0;
                Vector_t E (0.0), B (0.0);
                while ((((BGtag & BGphysics::Absorption) == BGphysics::Absorption) &&
                        ((BGtag & BGphysics::FNEmission) != BGphysics::FNEmission) &&
                        ((BGtag & BGphysics::SecondaryEmission) != BGphysics::SecondaryEmission))
                       ||
                       (fabs (E (0)) < eInitThreshold_m &&
                        fabs (E (1)) < eInitThreshold_m &&
                        fabs (E (2)) < eInitThreshold_m)) {
                    E = Vector_t (0.0);
                    B = Vector_t (0.0);
                    int tmp = (int)(IpplRandom () * num_triangles_m);
                    BGtag = TriBGphysicstag_m[tmp];
                    k = tmp;
                    Vector_t centroid (0.0);
                    itsOpalBeamline.getFieldAt (
                        Tribarycent_m[k] + darkinward * TriNormal_m[k],
                        centroid,
                        itsBunch->getdT (),
                        E,
                        B);
                }
                partsr_m.push_back (Tribarycent_m[k] + darkinward * TriNormal_m[k]);
            }
            Message* mess = new Message ();
            putMessage (*mess, partsr_m.size ());
            for (std::vector<Vector_t>::iterator myIt = partsr_m.begin ();
                 myIt != partsr_m.end ();
                 myIt++) {
                putMessage (*mess, *myIt);

            }
            Ippl::Comm->broadcast_all (mess, tag);
        } else {
            // receive particle position message
            size_t nData = 0;
            Message* mess = Ippl::Comm->receive_block (Parent, tag);
            getMessage (*mess, nData);
            for (size_t i = 0; i < nData; i++) {
                Vector_t tmp = Vector_t (0.0);
                getMessage (*mess, tmp);
                partsr_m.push_back (tmp);
            }
        }
    }
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
