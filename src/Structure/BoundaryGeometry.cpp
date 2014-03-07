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

using namespace Expressions;
using namespace Physics;

extern Inform* gmsg;

#define SQR(x) ((x)*(x))
#define PointID(triangle_id, vertex_id) allbfaces_m[4 * (triangle_id) + (vertex_id)]
#define Point(triangle_id, vertex_id)   geo3Dcoords_m[allbfaces_m[4 * (triangle_id) + (vertex_id)]]

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
    const double& incEnergy,
    const double& incQ,
    const Vector_t& incMomentum,
    PartBunch* itsBunch,
    double& seyNum
    ) {
    Inform msg ("BGphyDebug");
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
    const double& incEnergy,
    const double& incQ,
    const Vector_t& incMomentum,
    PartBunch* itsBunch,
    double& seyNum,
    const int& para_null
    ) {
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

/*
   helper functions for STL max/min_element
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
    double A,
    double B,
    int maxUlps ) {
                    
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    assert (maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    assert (sizeof (long long) == sizeof (int64_t) );
    assert (sizeof (long long) == sizeof (double) );
                    
    // Make [ab]Int lexicographically ordered as a twos-complement int
    double* pa = &A;
    int64_t aInt = *(int64_t*)pa;
    if (aInt < 0)
        aInt = 0x8000000000000000LL - aInt;
                    
    double* pb = &B;
    int64_t bInt = *(int64_t*)pb;
    if (bInt < 0)
        bInt = 0x8000000000000000LL - bInt;
                    
    int64_t intDiff = aInt - bInt;
    if (llabs(intDiff) <= maxUlps)
        return 0;
    return intDiff;
}


/*
  Map point to unique voxel ID.

  Remember:
  * hr_m:  is the  mesh size
  * nr_m:  number of mesh points
  */
int BoundaryGeometry::map_point_to_voxel_id (Vector_t x) {
    int id_tx = floor ((x[0] - mincoords_m[0]) / hr_m[0]);
    int id_ty = floor ((x[1] - mincoords_m[1]) / hr_m[1]);
    int id_tz = floor ((x[2] - mincoords_m[2]) / hr_m[2]);
    
    if (id_tx == -1) id_tx = 0;
    if (id_ty == -1) id_ty = 0;
    if (id_tz == -1) id_tz = 0;
    
    if (id_tx < 0 || id_ty < 0 || id_tz < 0) {
        return 0;
    }
    return 1 + id_tz * nr_m[0] * nr_m[1] + id_ty * nr_m[0] + id_tx;
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
                Vector_t x1 = bg->getPoint (i, 1);
                Vector_t x2 = bg->getPoint (i, 2);
                Vector_t x3 = bg->getPoint (i, 3);
                double length_edge1 = sqrt (
                    SQR (x1[0] - x2[0]) + SQR (x1[1] - x2[1]) + SQR (x1[2] - x2[2]));
                double length_edge2 = sqrt (
                    SQR (x3[0] - x2[0]) + SQR (x3[1] - x2[1]) + SQR (x3[2] - x2[2]));
                double length_edge3 = sqrt (
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
            bg->outside_point_m = bg->maxcoords_m + bg->hr_m;
            *gmsg << "* Geometry interval built done." << endl;
        }

        /*
          Make the boundary set by using
          * triangle vertices
          * bounding box vertices
          * several points in triangle central lines
          * several points in triangle edges
          */
        static void makeBoundaryIndexSet (BoundaryGeometry* bg) {
            class Local {
            public:
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
                  Get the smallest bounding box of triangle given by ID. The
                  smallest bounding box is used to make sure that all part of a triangle
                  is known by boundary bounding box.
                */
                static std::vector<Vector_t> getBBoxOfTriangle (BoundaryGeometry* bg, size_t id) {
                    Vector_t min = bg->getPoint (id, 1);
                    Vector_t max = min;
                    for (int i = 2; i <= 3; i++) {
                        Vector_t P = bg->getPoint (id, i);
                        if (P(0) < min[0]) min[0] = P(0);
                        if (P(1) < min[1]) min[1] = P(1);
                        if (P(2) < min[2]) min[2] = P(2);
                        if (P(0) > max[0]) max[0] = P(0);
                        if (P(1) > max[1]) max[1] = P(1);
                        if (P(2) > max[2]) max[2] = P(2);
                    }
                    std::vector<Vector_t> ret;
                    ret.push_back (min);
                    ret.push_back (max);
                    PAssert (ret.size () != 0);
                    return ret;
                }

                static void write_bbox_mesh (
                    std::set<size_t> ids,
                    Vector_t hr_m,
                    Vektor<int,3> nr,
                    Vector_t origin
                    ) {
                    /*-------------------------------------------------------------------------*/
                    size_t numpoints = 8 * ids.size ();
                    std::set<size_t>::iterator id;
                    std::ofstream of;
                    of.open (string ("data/testBBox.vtk").c_str ());
                    assert (of.is_open ());
                    of.precision (6);
    
                    of << "# vtk DataFile Version 2.0" << std::endl;
                    of << "generated using BoundaryGeometry::makeBoundaryIndexSet" << std::endl;
                    of << "ASCII" << std::endl << std::endl;
                    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
                    of << "POINTS " << numpoints << " float" << std::endl;
    
                    for (id = ids.begin (); id != ids.end (); id++) {
                        size_t k = (*id - 1) / (nr[0] * nr[1]);
                        size_t rest = (*id - 1) % (nr[0] * nr[1]);
                        size_t j = rest / nr[0];
                        size_t i = rest % nr[0]; 

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
    
                    of << "CELLS " << ids.size () << " " << 9 * ids.size () << std::endl;
                    for (size_t i = 0; i < ids.size (); i++)
                        of << "8 "
                           << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3 << " "
                           << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << std::endl;
                    of << "CELL_TYPES " << ids.size () << std::endl;
                    for (size_t i = 0; i <  ids.size (); i++)
                        of << "11" << std::endl;
                    of << "CELL_DATA " << ids.size () << std::endl;
                    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
                    of << "LOOKUP_TABLE " << "default" << std::endl;
                    for (size_t i = 0; i <  ids.size (); i++)
                        of << (float)(i) << std::endl;
                    of << std::endl;
                    of << "COLOR_SCALARS " << "BBoxColor " << 4 << std::endl;
                    for (size_t i = 0; i < ids.size (); i++) {
                        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << std::endl;
                    }
                    of << std::endl;
                }


            };


            for (int i = 0; i < bg->num_triangles_m; i++) {
                std::vector<Vector_t> coords;
                /* Discretize the three central lines and the three triangle edges to
                   200 segments to get a more complete boundary index set. */
                const Vector_t c1 = 0.5 * (bg->getPoint (i, 2) + bg->getPoint (i, 3)) - bg->getPoint (i, 1);
                const Vector_t c2 = 0.5 * (bg->getPoint (i, 3) + bg->getPoint (i, 1)) - bg->getPoint (i, 2);
                const Vector_t c3 = 0.5 * (bg->getPoint (i, 1) + bg->getPoint (i, 2)) - bg->getPoint (i, 3);

                const Vector_t e1 = bg->getPoint (i, 2) - bg->getPoint (i, 1);
                const Vector_t e2 = bg->getPoint (i, 3) - bg->getPoint (i, 2);
                const Vector_t e3 = bg->getPoint (i, 1) - bg->getPoint (i, 3);
                const int num_segments = 200;
                const double x = 1.0 / num_segments;
                for (int j = 1; j <= num_segments; j++) {
                    // discretize the three central lines.
                    coords.push_back (bg->getPoint (i, 1) + x * j * c1);
                    coords.push_back (bg->getPoint (i, 2) + x * j * c2);
                    coords.push_back (bg->getPoint (i, 3) + x * j * c3);
                    // discretize the three  triangle edges:
                    coords.push_back (bg->getPoint (i, 1) + x * j * e1);
                    coords.push_back (bg->getPoint (i, 2) + x * j * e2);
                    coords.push_back (bg->getPoint (i, 3) + x * j * e3);
                }

#if 0
                std::vector<Vector_t> ret = Local::getBBoxOfTriangle (bg, i);

                Vector_t min = ret[0];
                Vector_t max = ret[1];

                coords.push_back (Vector_t (min[0], min[1], min[2]));
                coords.push_back (Vector_t (min[0], min[1], max[2]));
                coords.push_back (Vector_t (min[0], max[1], max[2]));
                coords.push_back (Vector_t (min[0], max[1], min[2]));
                coords.push_back (Vector_t (max[0], min[1], min[2]));
                coords.push_back (Vector_t (max[0], min[1], max[2]));
                coords.push_back (Vector_t (max[0], max[1], max[2]));
                coords.push_back (Vector_t (max[0], max[1], min[2]));

                /*
                  :TODO: better description why we are doing this. Actually I (Achim)
                  don't understand it.

                  add some margins to make sure that boundary bounding box has both
                  positive and negtive margins w.r.t boundary. Just make sure that the
                  size of bounding box for triangle is smaller than the boundary bounding
                  box.
                */
                min -= bg->hr_m;
                max += bg->hr_m;
                Vector_t P;
                P = Vector_t (min[0], min[1], min[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (min[0], min[1], min[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (min[0], min[1], max[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (min[0], max[1], max[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (min[0], max[1], min[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->mincoords_m)) coords.push_back (P);
                P = Vector_t (max[0], min[1], min[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (max[0], min[1], max[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (max[0], max[1], max[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
                P = Vector_t (max[0], max[1], min[2]); if (Local::is_in_voxel (P, bg->mincoords_m, bg->maxcoords_m)) coords.push_back (P);
#endif
                std::vector<Vector_t>::iterator point;
                for (point = coords.begin (); point != coords.end (); point++) {
                    /*
                      if (!is_in_bbox (*point, mincoords_m, maxcoords_m))
                      continue;
                    */
                    int id = bg->map_point_to_voxel_id (*point);
                    assert (id > 0);
                    // insert ID to std::set! 
                    bg->boundary_ids_m.insert (id);
            
                    std::map< size_t, std::set<size_t> >::iterator It;
                    It = bg->CubicLookupTable_m.find (id);
                    if (It == bg->CubicLookupTable_m.end ()) {
                        std::set<size_t> tmp;
                        tmp.insert (i);
                        bg->CubicLookupTable_m.insert (std::pair <size_t, std::set<size_t> > (id, tmp));
                    } else
                        (*It).second.insert (i);
                }
            }
            if(Ippl::myNode() == 0) {
                Local::write_bbox_mesh (bg->boundary_ids_m, bg->hr_m, bg->nr_m, bg->mincoords_m);
            }
            *gmsg << "* Boundary index set built done." << endl;
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
            std::set<int> neighbors = bg->triangleNeighbors_m[triangle_id];

            for (std::set<int>::iterator triangle_iter = neighbors.begin();
                 triangle_iter != neighbors.end();
                 triangle_iter++) {
                if (!bg->isOriented_m [*triangle_iter])
                    orientTriangles (bg, triangle_id, *triangle_iter);
            }
        }

        /*
          Compute intersection point of line segment given by x and y 
          and triangle given by its ID.

          Returns true if the value returned in argument 'result' is an 
          intersection point. 

          Used in: isInside() 
        */

        static bool computeLineTriangleIntersectionPoint (
            BoundaryGeometry* bg,
            const Vector_t& x,
            const Vector_t& y,
            const size_t& triangle_id,
            Vector_t& result
            ) {
            const Vector_t t0 = bg->getPoint (triangle_id, 1);
            const Vector_t u = bg->getPoint (triangle_id, 2) - t0;
            const Vector_t v = bg->getPoint (triangle_id, 3) - t0;
            const Vector_t n = cross (u, v);

            const Vector_t lseg = y - x;
            if (fabs (dot (n, lseg)) < 1.0e-10) {
                // triangle is parallel to line segment
                return false;
            }

            const double rI = dot(n, t0-x) / dot(n, lseg);
            if((rI < 0) || (rI > 1)) {
                // intersection point is on the extended line
                return false;
            }
            result = x + rI * lseg;
            const Vector_t w = result - t0;
            const double tmp1 = dot (u, v);
            const double tmp2 = dot (w, v);
            const double tmp3 = dot (u, w);
            const double tmp4 = dot (u, u);
            const double tmp5 = dot (v, v);
            const double tmp6 = tmp1 * tmp1 - tmp4 * tmp5;
            const double sI = (tmp1 * tmp2 - tmp5 * tmp3) / tmp6;
            const double tI = (tmp1 * tmp3 - tmp4 * tmp2) / tmp6;
            return (sI >= 0) && (tI >= 0) && ((sI + tI) <= 1);
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
                if (computeLineTriangleIntersectionPoint (bg, x, y, triangle_id, result)) {
                    intersection_points.push_back (result);
                    num_intersections++;
                }
            }
            IpplTimings::stopTimer (bg->Tinward_m);
            return ((intersection_points.size () % 2) == 1);
        }


        static void computeTriangleNeighbors (BoundaryGeometry* bg) {
            std::set<int> adjacencies_to_pt  [bg->num_points_m];

            // for each triangles find adjacent triangles for each vertex
            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                for (int j = 1; j <= 3; j++) {
                    int pt_id = bg->PointID (triangle_id, j);
                    assert (pt_id < bg->num_points_m);
                    adjacencies_to_pt [pt_id].insert (triangle_id);
                }
            }

            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                std::set<int>  to_A = adjacencies_to_pt [bg->PointID (triangle_id, 1)];
                std::set<int>  to_B = adjacencies_to_pt [bg->PointID (triangle_id, 2)];
                std::set<int>  to_C = adjacencies_to_pt [bg->PointID (triangle_id, 3)];

                std::set<int>  intersect;
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

                bg->triangleNeighbors_m.insert (std::pair <int, std::set<int>> (triangle_id, intersect));
                                
            }
        }

        static inline Vector_t normalVector (BoundaryGeometry* const bg, const int triangle_id) {
            const Vector_t A = bg->getPoint (triangle_id, 1);
            const Vector_t B = bg->getPoint (triangle_id, 2);
            const Vector_t C = bg->getPoint (triangle_id, 3);

            /*
              compute triangle normal
             */
            const Vector_t N = cross (B - A, C - A);
            const double magnitude = sqrt (SQR (N (0)) + SQR (N (1)) + SQR (N (2)));
            assert (fcmp (magnitude, 0.0, 10) > 0); // in case we have degenerted triangles
            return N / magnitude;
        }

        static bool hasInwardPointingNormal (BoundaryGeometry* const bg, const int triangle_id, Vector_t& N) {
            N = normalVector (bg, triangle_id);

            // choose a point P close to the center of the triangle
            const Vector_t P = (A+B+C)/3 + N * 0.1 * bg->longest_side_min_m;

            /*
              The triangle normal points inward, if P is
              - outside the geometry and the dot product is negativ
              - or inside the geometry and the dot product is positiv

              Remember:
                The dot product is positiv only if both vectors are
                pointing in the same direction.
            */
            const bool is_inside = isInside (bg, P);
            const double dotPA_N = dot (P - A, N);
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

            Vector_t N;
            if (!hasInwardPointingNormal (bg, 0, N)) {
                N = -N;
                int id = bg->PointID (0, 2);
                bg->PointID (0, 2) = bg->PointID (0, 3);
                bg->PointID (0, 3) = id;
            }

            bg->TriNormal_m.push_back (N);

            bg->isOriented_m [0] = true;

            // compute edge-neightbors for each triangle
            computeTriangleNeighbors (bg);

            // orient all triangles, starting with the neighbors of the first
            std::set<int> neighbors = bg->triangleNeighbors_m[0];
            for (std::set<int>::iterator triangle_iter = neighbors.begin();
                 triangle_iter != neighbors.end();
                 triangle_iter++) {
                if (!bg->isOriented_m [*triangle_iter])
                    orientTriangles (bg, 0, *triangle_iter);
            }

#if 0
            // for debugging only!
            for (int triangle_id = 0; triangle_id < bg->num_triangles_m; triangle_id++) {
                if (!hasInwardPointingNormal (bg, triangle_id)) {
                    *gmsg << "* Wrong oriented triangle: " << triangle_id << endl;
                }
            }
#endif
            // compute inward-normals
            bg->TriNormal_m.reserve (bg->num_triangles_m);
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

    Local::makeBoundaryIndexSet (this);
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
int BoundaryGeometry::PartInside (
    const Vector_t r,
    const Vector_t v,
    const double dt,
    int Parttype,
    const double Qloss,
    Vector_t& intecoords,
    int& triId,
    double& Energy
    ) {
    class Local {

    public:
        /*
          Find a intersection between a line segment and triangle,faster by using
          the pre-computed oriented normal.
     
          @param x0         start of line segment.
          @param x1         end of line segment
          @param i          triangle ID

          Algorithms:
          1) find the intersection between line segment \f$\vec{x1-x0}\f$ and plane
          defined by point t0 and triangle normal;
          if the dot product of line segment and plane normal equals to zero and
          x0 is not the barycentric point of the triangle(in the plane), then
          the line segment is parallel to plane
          return no intersection,
          else if particle is really move (\f$ x0 \neq x1 \f$), then
          return  initialized position-triangle barycentric point as intersection.

          The intersection position rI w.r.t x0 is obtained from:

          \f$ rI=\frac{\vec{n} \cdot \vec{(t0-x0)}}{\vec{n} \cdot \vec{(x1-x0)}} \f$,

          where t0 is the first vertex of triangle, n is the normal of triangle.
          The intersection point Itsec is obtained from:
          \f$ Itsec=x0+rI(x1-x0) \f$.

          2) check if the intersection point is inside the triangle by using parametric
          coordinates sI and tI of the intersecion point.
          First calculate sI and tI. The parametric plane equation is given by:
          \f$ t(sI,tI)=t0+sI(t1-t0)+tI(t2-t0)=t0+sI\vec{u}+tI\vec{v} \f$.
          \f$\vec{w}=\vec{Itsec-t0}\f$ is also in the plane, solve equation:
          \f$\vec{w}=t0+sI\vec{u}+tI\vec{v}\f$ , we obtain the sI and tI.
          \f$ sI=\frac{(\vec{u} \cdot \vec{v})(\vec{w} \cdot \vec{v})-(\vec{v} \cdot \vec{v})(\vec{w} \cdot \vec{u})}{(\vec{u} \cdot \vec{v})^2-(\vec{u} \cdot \vec{u})(\vec{v} \cdot \vec{v})} \f$,
          \f$ tI=\frac{(\vec{u} \cdot \vec{v})(\vec{w} \cdot \vec{u})-(\vec{u} \cdot \vec{u})(\vec{w} \cdot \vec{v})}{(\vec{u} \cdot \vec{v})^2-(\vec{u} \cdot \vec{u})(\vec{v} \cdot \vec{v})} \f$.
          If \f$ sI \geq 0 \f$, \f$ tI \geq 0 \f$ and \f$ sI+tI \leq 1 \f$, then
          the intersection is inside the triangle, and return the intersection
          coordinate Itsec.
        */
        static bool FindIntersection (
            BoundaryGeometry* bg,
            const Vector_t& x,          // [in] start of line segment
            const Vector_t& y,          // [in] end of line segment
            const size_t& triangle_id,  // [in] triangle ID
            double& rI,
            Vector_t& intersection_pt   // [out] intersection point
            ) {
            IpplTimings::startTimer (bg->TRayTrace_m);

            bool result = false;
            const Vector_t t0 = bg->getPoint (triangle_id, 1);
            const Vector_t u = bg->getPoint (triangle_id, 2) - t0;
            const Vector_t v = bg->getPoint (triangle_id, 3) - t0;
            const Vector_t lt = t0 - x;
            const Vector_t n = bg->TriNormal_m[triangle_id];

            const Vector_t lseg = y - x; // length and direction of line segment;
            const double dotLT = dot (n, lseg);
            if (fabs (dotLT) < 1.0e-10) {
                if ((x == bg->Tribarycent_m[triangle_id]) && (x != y)) {
                    /*
                      Some initialized particles have momenta parallel to its
                      triangle normal, this kind of particles will lose
                      directly
                    */
                    intersection_pt = bg->Tribarycent_m[triangle_id];
                    return true;
                }
            } else {
                // find intersection position w.r.t x and the unit is (y-x);
                rI = dot (n, lt) / dotLT;

                // find the coordinate of intersection plane.
                const Vector_t ItSec = x + rI * lseg;

                // test if intersection is inside the triangle.
                const Vector_t w = ItSec - t0;
                const double tmp1 = dot (u, v);
                const double tmp2 = dot (w, v);
                const double tmp3 = dot (u, w);
                const double tmp4 = dot (u, u);
                const double tmp5 = dot (v, v);
                const double temp = (tmp1 * tmp1 - tmp4 * tmp5);
                const double sI = (tmp1 * tmp2 - tmp5 * tmp3) / temp;
                const double tI = (tmp1 * tmp3 - tmp4 * tmp2) / temp;
                if ((sI >= 0.0) && (tI >= 0.0) && ((sI + tI) <= 1.0)) {
                    intersection_pt = ItSec;
                    result = true;
                }
            }
            IpplTimings::stopTimer (bg->TRayTrace_m);
            return result;
        }
    };

    int ret = -1;
    const double p_sq = dot (v, v);
    const double betaP = 1.0 / sqrt (1.0 + p_sq);

    const Vector_t temp1 = r; //particle position in timestep n;
    const Vector_t temp = r + (c * betaP * v * dt); //particle position in tstep n+1;
    double rI = 0.0;
    Vector_t intersection_pt = outside_point_m;

    IpplTimings::startTimer (TPInside_m);

    /* test if particle position in timestep n is inside the cubic bounding box.
       If true, do the following tests */
    int id;
    if (boundary_ids_m.find (map_point_to_voxel_id (temp1)) != boundary_ids_m.end ()) {
        // particle is in geometry at timestep n
        id = map_point_to_voxel_id (temp1);
    } else if (boundary_ids_m.find (map_point_to_voxel_id (temp)) != boundary_ids_m.end ()) {
        // particle is in geometry at timestep n+1
        id = map_point_to_voxel_id (temp);
    } else {
        goto out;
    }
    { /* we need this brace! Otherwise the compiler (gcc) complains about
         initializing nx, ny und idc after above goto statement. */

        /* Build an array containing the IDs of 27(3*3*3) voxels. The ID
           of the voxel containing the particle, is the center of these
           voxels. Just for the situation that even the line segment
           cross more than one voxel.*/
        int nx = nr_m[0];
        int ny = nr_m[1];
        int idc[27] = {
            id,                    id + 1,                id - 1,
            id + nx,               id - nx,               id + nx + 1,
            id + nx - 1,           id - nx - 1,           id - nx + 1,
            id + nx * ny,          id + nx * ny + 1,      id + nx * ny - 1,
            id + nx * ny + nx,     id + nx * ny - nx,     id + nx * ny + nx + 1,
            id + nx * ny + nx - 1, id + nx * ny - nx - 1, id + nx * ny - nx + 1,
            id - nx * ny,          id - nx * ny + 1,      id - nx * ny - 1,
            id - nx * ny + nx,     id - nx * ny - nx,     id - nx * ny + nx + 1,
            id - nx * ny + nx - 1, id - nx * ny - nx - 1, id - nx * ny - nx + 1
        };
        
        /* Test all the 27 voxels to find if the line segment has
           intersection with the triangles in those voxels. */
        for (int k = 0; k < 27; k++) {
            std::map< size_t, std::set<size_t> >::iterator It;
            It = CubicLookupTable_m.find (idc[k]);
            if (It == CubicLookupTable_m.end ())
                continue; // not a voxel
            
            // for each triangle in this voxel
            std::set<size_t> ::iterator faceIt;
            for (faceIt = (*It).second.begin ();
                 faceIt != (*It).second.end ();
                 faceIt++) {
                if (v != 0 && dot (v, TriNormal_m[*faceIt]) <= 0.0) {
                    /* If the particle have a momenta greater zero with opposite
                       direction to triangle normal, do the following tests. */
                    if (Local::FindIntersection (
                            this,
                            temp1,      // IN: particle position in tstep n
                            temp,       // IN: particle position in tstep n+1
                            *faceIt,    // IN: triangle id
                            rI,         // OUT: ratio
                            intersection_pt // OUT: intersection points
                            )) {       
                        /* Test if the intersection is between the particle
                           position in tstep n and particle position in
                           tstep n+1 or is in the extension of line segment
                           when particle position in tstep n is already
                           outside the geometry( this may be not accurate
                           and may be the source of problem.) */
                        if ((rI >= -0.00001 && rI <= 1.00001) ||
                            (rI < 0 && dot (temp1 - Tribarycent_m[*faceIt], TriNormal_m[*faceIt]) <= 0.0)) {
                            intecoords = intersection_pt;
                            triId = (*faceIt);
                            assert (dot (TriNormal_m[*faceIt], v) < 0 || intersection_pt == temp1);
                            // energy in eV
                            Energy = Physics::m_e * (sqrt (1.0 + p_sq) - 1.0) * 1.0e9;
                            if (Parttype == 0)
                                TriPrPartloss_m[*faceIt] += Qloss;
                            else if (Parttype == 1)
                                TriFEPartloss_m[*faceIt] += Qloss;
                            else
                                TriSePartloss_m[*faceIt] += Qloss;
                            if (TriPrPartloss_m[*faceIt] > 0 ||
                                TriSePartloss_m[*faceIt] > 0 ||
                                TriFEPartloss_m[*faceIt] > 0) {
                                ;
                            }
                            ret = 0;
                            goto out;
                        }
                    }
                }
            } // end for all triangles
        } // end for all voxels
    }
out:
    IpplTimings::stopTimer (TPInside_m);
    return ret;
}

void BoundaryGeometry::updateElement (ElementBase* element) {
}

void BoundaryGeometry::writeGeomToVtk (string fn) {
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

Inform& BoundaryGeometry::printInfo (Inform& os) const {
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
       << "* Size of boundary index set " << boundary_ids_m.size () << '\n'
       << "* Number of all boxes        " << nr_m (0) * nr_m (1) * nr_m (2) << '\n'
        << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
