//
// vi: set et ts=4 sw=4 sts=4:

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

BoundaryGeometry::BoundaryGeometry() :
    Definition (
        SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry."),
    allbfaces_m (NULL),
    Tribarycent_m (NULL),
    TriPrPartloss_m (NULL),
    TriSePartloss_m (NULL),
    TriFEPartloss_m (NULL),
    TriPrPartlossZ_m (NULL),
    TriSePartlossZ_m (NULL),
    TriFEPartlossZ_m (NULL) {
    itsAttr[FGEOM] = Attributes::makeString
        ("FGEOM",
         "Specifies the geometry file [h5fed]",
         "");

    itsAttr[TOPO] = Attributes::makeString
        ("TOPO",
         "BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written ",
         "ELLIPTIC");

    itsAttr[LENGHT] = Attributes::makeReal
        ("LENGHT",
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

    itsAttr[ZSHIFT] = Attributes::makeReal
        ("ZSHIFT",
         "Shift in z direction",
         0.0);

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
    if (!h5FileName_m.empty ())
        initialize ();
}

BoundaryGeometry::BoundaryGeometry(
    const string& name,
    BoundaryGeometry* parent
    ) :
    Definition (name, parent),
    allbfaces_m (NULL),
    Tribarycent_m (NULL),
    TriPrPartloss_m (NULL),
    TriSePartloss_m (NULL),
    TriFEPartloss_m (NULL),
    TriPrPartlossZ_m (NULL),
    TriSePartlossZ_m (NULL),
    TriFEPartlossZ_m (NULL) {
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
       if (TriPrPartlossZ_m)
          delete TriPrPartlossZ_m;
       if (TriSePartlossZ_m)
          delete TriSePartlossZ_m;
       if (TriFEPartlossZ_m)
          delete TriFEPartlossZ_m;
       if (TriPrPartloss_m )
       delete TriPrPartloss_m ;
       if (TriFEPartloss_m)
       delete TriFEPartloss_m;
       if (TriSePartloss_m)
       delete TriSePartloss_m;
     */
}

bool BoundaryGeometry::canReplaceBy (Object* object) {
    // Can replace only by another GEOMETRY.
    return dynamic_cast<Geometry*>(object) != 0;
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
            if ((intecoords[0] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](0)) ||
                (intecoords[1] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](1)) ||
                (intecoords[2] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](2))) {
                idx = 4 * triId + 1; // intersection is not the 1st vertex
            } else {
                idx = 4 * triId + 2; // intersection is the 1st vertex
            }
            sec_phys_m.nSec (incEnergy,
                             cosTheta,
                             seBoundaryMatType_m,
                             se_Num,
                             seType,
                             incQ,
                             TriNormal_m[triId],
                             intecoords,
                             geo3Dcoords_m[allbfaces_m[idx]],
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
            if ((intecoords[0] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](0)) ||
                (intecoords[1] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](1)) ||
                (intecoords[2] != geo3Dcoords_m[allbfaces_m[4 * triId + 1]](2))) {
                // intersection is not the 1st vertex
                idx = 4 * triId + 1;
            } else {
                // intersection is the 1st vertex
                idx = 4 * triId + 2;
            }
            sec_phys_m.nSec (incEnergy,
                             cosTheta,
                             se_Num,
                             seType,
                             incQ,
                             TriNormal_m[triId],
                             intecoords,
                             geo3Dcoords_m[allbfaces_m[idx]],
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
    for (int i = 0; i < numbfaces_global_m; i++) {
        if ((TriBGphysicstag_m[i] & BGphysics::FNEmission) == BGphysics::FNEmission) {
            Vector_t E (0.0), B (0.0);
            Vector_t centroid (0.0);
            itsOpalBeamline.getFieldAt (Tribarycent_m[i], centroid, t, E, B);
            double Enormal = dot (TriNormal_m[i], E);
            /* Enormal should be negative as E field direction should be
               opposite to inward normal of surface */
            if (Enormal < fieldFNthreshold_m) {
                std::vector<Vector_t> vertex;
                vertex.push_back (geo3Dcoords_m[allbfaces_m[4 * i + 1]]);
                vertex.push_back (geo3Dcoords_m[allbfaces_m[4 * i + 2]]);
                vertex.push_back (geo3Dcoords_m[allbfaces_m[4 * i + 3]]);
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
   Determine if a line segment  has intersects with a triangle.
   Algorithm
      :http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
   @param x0 and @param x1 define the line segment.
   @param i is the id of the triangle.
 */

/**
   Determine if a point is outside, inside or just on the boundary.

   @param x stands for coordinates of the test point.

   @result false point is outside boundary
   @result true point is inside boundary

   The basic idea is if a line segment starting from the test point has
   odd intersects with a closed boundary, then the test point is inside
   the geometry;
   if the intersects have even number, then the test points
   is outside the geometry;
   if the test point is amoung the intersects, then the test point is just
   on the boundary. Makesure the end point of the line
   segment is outside the geometry boundary.
 */
bool BoundaryGeometry::isInside (Vector_t x) {
    Vector_t x0 = x;
    Vector_t x1;
    x1[0] = x0[0];
    //x1[1] = x0[1];
    RANLIB_class* rGen = new RANLIB_class (265314159, 4);
    x1[1] = maxcoords_m[1] * (1.1 + rGen->uniform (0.0, 1.0));
    x1[2] = maxcoords_m[2] * (1.1 + rGen->uniform (0.0, 1.0));
    //x1[2] = x0[2];
    delete rGen;

    /*
       Random number could avoid some specific situation,
       like line parallel to boundary......
       x1 could be any point outside the boundary ;
     */
    IpplTimings::startTimer (Tinward_m);
    std::vector<Vector_t> IntesecNum = PartBoundaryInteNum (x0, x1);
    IpplTimings::stopTimer (Tinward_m);
    if (IntesecNum[0] == x0) {
        return true; // x0 is just on the boundary;
    } else {
        if (((IntesecNum.size () % 2) == 0) || (*IntesecNum.begin () == x1)) {
            return false; // x0 is  outside the boundary;
        } else
            return true;  // x0 is inside the boundary;
    }
}

std::vector<Vector_t>  BoundaryGeometry::GridIntersection (Vector_t x0, Vector_t x1) {

    std::vector<Vector_t> SegDiscrete;
    std::vector<Vector_t> Isp;
    std::vector<Vector_t>::iterator TriIscIt;

    std::vector<int> TriId;
    std::vector<int>::iterator TriIdIt;
    double hr_tmp = hr_m[0];
    if (hr_m[1] < hr_tmp)
        hr_tmp = hr_m[1];
    if (hr_m[2] < hr_tmp)
        hr_tmp = hr_m[2];  //interval should be smaller than the size of
                           // box;
    //hr_tmp*=0.001;
    int Seglen = (((int)floor (sqrt (dot (x0 - x1, x0 - x1)) / hr_tmp)) + 1);
    int count = 0;
    for (int i = 0; i < Seglen; i++) {
        SegDiscrete.push_back (x0 + hr_tmp * i * (x1 - x0) / sqrt (dot (x0 - x1, x0 - x1)));
    }
    SegDiscrete.push_back (x1);

    for (std::vector<Vector_t>::iterator myit = SegDiscrete.begin ();
         myit != SegDiscrete.end ();
         myit++) {
        size_t id = f (*myit);

        std::vector<size_t> ::iterator idIt;
        std::map< size_t, std::vector<size_t> >::iterator It = CubicLookupTable_m.find (id);
        if (It == CubicLookupTable_m.end ())
            continue;

        // cubic box is the boundary bounding box
        for (idIt = (*It).second.begin (); idIt != (*It).second.end (); idIt++) {
            Vector_t tmp =  LineInsTri (x0, x1, *idIt);
            TriIdIt = std::find (TriId.begin (), TriId.end (), *idIt);
            if (tmp != x1 &&
                (TriIdIt == TriId.end () || TriId.size () == 0)) {
                TriIscIt = std::find (Isp.begin (), Isp.end (), tmp);
                if (TriIscIt == Isp.end () || Isp.size () == 0) {
                    TriId.push_back (*idIt);
                    Isp.push_back (tmp);
                    ++count;
                }
            }
        }
    }
    if (count == 0)
        Isp.push_back (x1);
    return Isp;
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

void BoundaryGeometry::initialize () {

    h5_int64_t rc;

    *gmsg << "* Iniitializing Boundary Geometry ... ..." << endl;
    IpplTimings::startTimer (TPreProc_m);

    xyzscale_m = getXYZScale ();
    *gmsg << "* Scale all points of the geometry by " << xyzscale_m << endl;

    rc = H5SetErrorHandler (H5AbortErrorhandler);
    if (rc != H5_SUCCESS)
        ERRORMSG ("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    H5SetVerbosityLevel (1);
    h5_file_t* f = H5OpenFile (h5FileName_m.c_str (), H5_O_RDONLY, Ippl::getComm());
    h5t_mesh_t* m = NULL;
    H5FedOpenTriangleMesh (f, "0", &m);
    H5FedSetLevel (m, 0);

    numbfaces_global_m = H5FedGetNumElementsTotal (m);
    allbfaces_m = new int[numbfaces_global_m * 4];
    Tribarycent_m = new Vector_t[numbfaces_global_m];

    // iterate over all co-dim 0 entities, i.e. elements
    h5_loc_id_t local_id;
    int i = 0;
    h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
    while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
        h5_loc_id_t local_vids[4];
        H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
        allbfaces_m[4 * i]   = 0;
        allbfaces_m[4 * i + 1] = local_vids[0];
        allbfaces_m[4 * i + 2] = local_vids[1];
        allbfaces_m[4 * i + 3] = local_vids[2];
        i++;
    }
    H5FedEndTraverseEntities (iter);

    // loop over all vertices
    numpoints_global_m = H5FedGetNumVerticesTotal (m);
    double* point_coords = new double[3 * numpoints_global_m];
    for (i = 0; i < numpoints_global_m; i++) {
        h5_float64_t P[3];
        H5FedGetVertexCoordsByIndex (m, i, P);
        point_coords[i * 3]   = P[0] * xyzscale_m;
        point_coords[i * 3 + 1] = P[1] * xyzscale_m;
        point_coords[i * 3 + 2] = P[2] * xyzscale_m;
    }
    H5FedCloseMesh (m);
    H5CloseFile (f);

    Vector_t geo3d_tmp;

    double zshift = getZshift ();

    for (int i = 0; i < numpoints_global_m; i++) {
        geo3d_tmp[0] = point_coords[3 * i];
        geo3d_tmp[1] = point_coords[3 * i + 1];
        geo3d_tmp[2] = point_coords[3 * i + 2] + zshift;
        geo3Dcoords_m.push_back (geo3d_tmp);
    }
    *gmsg << "*  Vertex built done." << endl;

    mincoords_m = get_min_extend (geo3Dcoords_m);
    TriPrPartloss_m = new double[numbfaces_global_m];
    TriFEPartloss_m = new double[numbfaces_global_m];
    TriSePartloss_m = new double[numbfaces_global_m];
    for (int i = 0; i < numbfaces_global_m; i++) {
        Tribarycent_m[i] = (geo3Dcoords_m[allbfaces_m[4 * i + 1]]
                            + geo3Dcoords_m[allbfaces_m[4 * i + 2]]
                            + geo3Dcoords_m[allbfaces_m[4 * i + 3]]) / 3.0;
        Triarea_m.push_back (TriangleArea (i));

        TriPrPartloss_m[i] = 0.0;
        TriFEPartloss_m[i] = 0.0;
        TriSePartloss_m[i] = 0.0;
    }
    if (point_coords)
        delete point_coords;
    *gmsg << "*  Triangle barycent built done." << endl;

    maxcoords_m = get_max_extend (geo3Dcoords_m);

    len_m = maxcoords_m - mincoords_m;
    getMaxDimenssion (); // get maximum and minimum of triangle side.

    /*
       In principal, the number of discretization nr_m is maximum lenth in
       each dimension divided by maximum of triangle length. But if the hot
       spot, i.e., the multipacting/field emission zone is too small that
       normal bounding box covers the whole hot spot, the expensive
       triangle-line intersection tests will be frequently called. In these
       cases, we need to use smaller bounding box size to speed up
          simulation
       by setting a larger nr_m, for example use nr_m(0) /= 0.2 in the
          following.

       Todo:
       The relation between bounding box size and simulation time step &
       geometry shape maybe need to be summarized and modeled in a more
       flexible manner and could be adjusted in input file.
     */
    nr_m (0) = (int)floor (len_m (0) / triangle_max_m / 0.2); //2.0
    nr_m (1) = (int)floor (len_m (1) / triangle_max_m / 0.2);
    nr_m (2) = (int)floor (len_m (2) / triangle_max_m / 0.1); //2.0

    hr_m = len_m / nr_m;
    out_m = maxcoords_m + hr_m;
    TriPrPartlossZ_m = new double[nr_m (2)];
    TriSePartlossZ_m = new double[nr_m (2)];
    TriFEPartlossZ_m = new double[nr_m (2)];
    *gmsg << "*  Geometry interval built done." << endl;

    makeBoundaryIndexSet ();
    *gmsg << "*  BoundaryIndexSet built done." << endl;

    makeTriNormal ();
    *gmsg << "*  Triangle Normal built done." << endl;

    setBGphysicstag ();
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
                int tmp = (int)(IpplRandom () * numbfaces_global_m);
                BGtag = TriBGphysicstag_m[tmp];
                k = tmp;
                Vector_t centroid (0.0);
                itsOpalBeamline.getFieldAt (Tribarycent_m[k] + darkinward * TriNormal_m[k], centroid, itsBunch.getdT (), E, B);
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
               parallel
               plate. There is a distance of 0.01*length in x direction as
               margin. */
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

                    int k = (int)(IpplRandom () * numbfaces_global_m);
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
                    int k = (int)(IpplRandom () * numbfaces_global_m);
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
            for (std::vector<Vector_t>::iterator myIt = partsr_m.begin (), myItp = partsp_m.begin ();
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
                    int tmp = (int)(IpplRandom () * numbfaces_global_m);
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
   Make the boundary set by using triangle vertex, bounding box vertex,
   as well as points in triangle central lines.
 */
void BoundaryGeometry::makeBoundaryIndexSet () {
    std::set<size_t>::iterator bbIt;
    std::vector<Vector_t> BBox;
    for (int i = 0; i < numbfaces_global_m; i++) {
        std::vector<Vector_t> discreteTri;
        std::vector<Vector_t>::iterator disIt;
        /* Discrete three central lines and three triangle sides to
           200 segments to get a more complete boundary index set. */
        for (int j = 0; j < 200; j++) {
            // Discrete three central lines.
            discreteTri.push_back (
                getVertexCoord (i, 1)
                + 0.005 * (j + 1) * (0.5 * (getVertexCoord (i, 2) + getVertexCoord (i, 3)) - getVertexCoord (i, 1))
                );
            discreteTri.push_back (
                getVertexCoord (i, 2)
                + 0.005 * (j + 1) * (0.5 * (getVertexCoord (i, 3) + getVertexCoord (i, 1)) - getVertexCoord (i, 2))
                );
            discreteTri.push_back (
                getVertexCoord (i, 3)
                + 0.005 * (j + 1) * (0.5 * (getVertexCoord (i, 1) + getVertexCoord (i, 2)) - getVertexCoord (i, 3))
                );
            // Discrete three  triangle sides.
            discreteTri.push_back (
                getVertexCoord (i, 1)
                + 0.005 * (j + 1) * (getVertexCoord (i, 2) - getVertexCoord (i, 1))
                );
            discreteTri.push_back (
                getVertexCoord (i, 2)
                + 0.005 * (j + 1) * (getVertexCoord (i, 3) - getVertexCoord (i, 2))
                );
            discreteTri.push_back (
                getVertexCoord (i, 3)
                + 0.005 * (j + 1) * (getVertexCoord (i, 1) - getVertexCoord (i, 3))
                );
        }
        discreteTri.push_back (Tribarycent_m[i]);
        for (disIt = discreteTri.begin (); disIt != discreteTri.end (); disIt++) {
            int id = f (*disIt);

            //  assert(id > 0);
            bbIt = boundary_ids_m.find (id);
            if ((bbIt == boundary_ids_m.end ())) {
                Vector_t temp;
                temp[0] = floor (((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];
                temp[1] = floor (((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor (((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back (temp);
                boundary_ids_m.insert (id);
            }

            std::vector<size_t> tmp;
            std::map< size_t, std::vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find (id);
            if (It == CubicLookupTable_m.end ()) {
                tmp.push_back (i);
                CubicLookupTable_m.insert (std::pair <size_t, std::vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back (i);
        }
        std::vector<Vector_t> ret, cubic_coords;
        ret = SetMinMaxBound (i);
        cubic_coords.push_back (ret[0]);
        cubic_coords.push_back ((ret[0](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back ((ret[0](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back ((ret[0](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back ((ret[1](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back ((ret[1](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back ((ret[1](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back ((ret[1](0), ret[0](1), ret[0](2)));

        cubic_coords.push_back (ret[2]);
        cubic_coords.push_back ((ret[2](0), ret[2](1), ret[3](2)));
        cubic_coords.push_back ((ret[2](0), ret[3](1), ret[3](2)));
        cubic_coords.push_back ((ret[2](0), ret[3](1), ret[2](2)));
        cubic_coords.push_back ((ret[3](0), ret[3](1), ret[2](2)));
        cubic_coords.push_back ((ret[3](0), ret[3](1), ret[3](2)));
        cubic_coords.push_back ((ret[3](0), ret[2](1), ret[3](2)));
        cubic_coords.push_back ((ret[3](0), ret[2](1), ret[2](2)));

        for (disIt = cubic_coords.begin (); disIt != cubic_coords.end (); disIt++) {
            int id = f (*disIt);
	    // :FIXME: id == 0 happens, why and when?
            // assert(id > 0);
            bbIt = boundary_ids_m.find (id);
            if ((bbIt == boundary_ids_m.end ())) {
                Vector_t temp;
                // this may cause the increct display.
                temp[0] = floor (((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];
                temp[1] = floor (((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor (((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back (temp);
                boundary_ids_m.insert (id);
            }

            std::vector<size_t> tmp;
            std::map< size_t, std::vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find (id);
            if (It == CubicLookupTable_m.end ()) {
                tmp.push_back (i);
                CubicLookupTable_m.insert (std::pair <size_t, std::vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back (i);
        }
    }

    /*-------------------------------------------------------------------------*/
    size_t numpoints = 8 * BBox.size ();
    std::vector<Vector_t>::iterator bIt;
    std::ofstream of;
    of.open (string ("data/testBBox.vtk").c_str ());
    assert (of.is_open ());
    of.precision (6);

    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated using BoundaryGeometry::makeBoundaryIndexSet" << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
    of << "POINTS " << numpoints << " float" << std::endl;

    for (bIt = BBox.begin (); bIt != BBox.end (); bIt++) {
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << std::endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << std::endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << std::endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << std::endl;
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << std::endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << std::endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] + hr_m[2] << std::endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1] << " " << (*bIt)[2] + hr_m[2] << std::endl;
    }
    of << std::endl;

    of << "CELLS " << BBox.size () << " " << 9 * BBox.size () << std::endl;
    for (size_t i = 0; i < BBox.size (); i++)
        of << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
           << " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << std::endl;
    of << "CELL_TYPES " << BBox.size () << std::endl;
    for (size_t i = 0; i <  BBox.size (); i++)
        of << "11" << std::endl;
    of << "CELL_DATA " << BBox.size () << std::endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
    of << "LOOKUP_TABLE " << "default" << std::endl;
    for (size_t i = 0; i <  BBox.size (); i++)
        of << (float)(i) << std::endl;
    of << std::endl;
    of << "COLOR_SCALARS " << "BBoxColor " << 4 << std::endl;
    for (size_t i = 0; i < BBox.size (); i++) {
        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << std::endl;
    }
    of << std::endl;
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
    int ret;
    const double p_sq = dot (v, v);
    const double betaP = 1.0 / sqrt (1.0 + p_sq);

    const Vector_t temp1 = r; //particle position in tstep n;
    const Vector_t temp = r + (c * betaP * v * dt); //particle position in tstep n+1;
    double rI = 0.0;
    Vector_t Isc = out_m;

    IpplTimings::startTimer (TPInside_m);

    /* test if particle position in tstep n is inside the cubic bounding box.
       If true, do the following tests */
    int flg = 0;
    int id;
    if (isInGeometry (temp1)) {
        id = f (temp1);
    } else if (isInGeometry (temp)) {
        id = f (temp);
    } else {
        // cannot use goto out here (compiler complains)
        IpplTimings::stopTimer (TPInside_m);
        return -1;
    }

    /* Build an array containing the IDs of 27(3*3*3) cubic boxes. The ID
       of the cubic box which contains the particle, is in the center of
       all this boxes. Just for the situation that even the line segment
       cross more than one cubic box.*/
    int idc[27] = {
        id,                                 id + 1,                                 id - 1,
        id + nr_m (0),                       id - nr_m (0),                           id + nr_m (0) + 1,
        id + nr_m (0) - 1,                   id - nr_m (0) - 1,                       id - nr_m (0) + 1,
        id + nr_m (0) * nr_m (1),               id + nr_m (0) * nr_m (1) + 1,               id + nr_m (0) * nr_m (1) - 1,
        id + nr_m (0) * nr_m (1) + nr_m (0),     id + nr_m (0) * nr_m (1) - nr_m (0),         id + nr_m (0) * nr_m (1) + nr_m (0) + 1,
        id + nr_m (0) * nr_m (1) + nr_m (0) - 1, id + nr_m (0) * nr_m (1) - nr_m (0) - 1,     id + nr_m (0) * nr_m (1) - nr_m (0) + 1,
        id - nr_m (0) * nr_m (1),               id - nr_m (0) * nr_m (1) + 1,               id - nr_m (0) * nr_m (1) - 1,
        id - nr_m (0) * nr_m (1) + nr_m (0),     id - nr_m (0) * nr_m (1) - nr_m (0),         id - nr_m (0) * nr_m (1) + nr_m (0) + 1,
        id - nr_m (0) * nr_m (1) + nr_m (0) - 1, id - nr_m (0) * nr_m (1) - nr_m (0) - 1,     id - nr_m (0) * nr_m (1) - nr_m (0) + 1
    };

    /* Test all the 27 cubic boxes to find if the line segment has
       intersection with the triangles in those cubic boxes. */
    for (int k = 0; k < 27; k++) {
        std::map< size_t, std::vector<size_t> >::iterator It = CubicLookupTable_m.find (idc[k]);
        if (It == CubicLookupTable_m.end ())
            continue; // not a boundary box

        // for each triangle in this boundary box
        std::vector<size_t> ::iterator faceIt;
        for (faceIt = (*It).second.begin (); faceIt != (*It).second.end (); faceIt++) {
            if (v != 0 && dot (v, TriNormal_m[*faceIt]) <= 0.0) {
                /* If the particle have none zero momenta and momenta
                   has opposite direction with triangle normal, do
                   the following tests. */
                FindIntersection (
                    temp1,      // IN: particle position in tstep n
                    temp,       // IN: particle position in tstep n+1
                    *faceIt,      // IN: triangle id
                    rI,         // OUT: ratio
                    Isc);       // OUT: intersection points
                if (Isc != out_m) {
                    /* Test if the intersection is between the particle
                       position in tstep n and particle position in
                       tstep n+1 or is in the extension of line segment
                       when particle position in tstep n is already
                       outside the geometry( this may be not accurate
                       and may be the source of problem.) */
                    if ((rI >= -0.00001 && rI <= 1.00001) ||
                        (rI < 0 && dot (temp1 - Tribarycent_m[*faceIt], TriNormal_m[*faceIt]) <= 0.0)) {
                        flg++;
                        intecoords = Isc;
                        triId = (*faceIt);
                        assert (dot (TriNormal_m[*faceIt], v) < 0 || Isc == temp1);
                        Energy = Physics::m_e * (sqrt (1.0 + p_sq) - 1.0) * 1.0e9; //in eV
                        if (Parttype == 0)
                            TriPrPartloss_m[*faceIt] += Qloss;
                        else if (Parttype == 1)
                            TriFEPartloss_m[*faceIt] += Qloss;
                        else
                            TriSePartloss_m[*faceIt] += Qloss;
                        if (TriPrPartloss_m[*faceIt] > 0 ||
                            TriSePartloss_m[*faceIt] > 0 ||
                            TriFEPartloss_m[*faceIt] > 0) {
                            std::cout << "* Loss Data" << *faceIt << " : "
                                      << TriPrPartloss_m[*faceIt] << " "
                                      << TriSePartloss_m[*faceIt] << " "
                                      << TriFEPartloss_m[*faceIt]
                                      << " qloss: " << Qloss << std::endl;
                        }
                        break;
                    }
                }
            }
        } // end for all triangles
        if (flg != 0) {
            ret = 0;
            goto out;
        }
    } // end for all 27 boxes
    if (flg == 0) {
        ret = - 1;
    }

out:
    IpplTimings::stopTimer (TPInside_m);
    return ret;
}

void BoundaryGeometry::updateElement (ElementBase* element) {
}

void BoundaryGeometry::writeGeomToVtk (string fn) {
    if (Ippl::myNode () == 0) {
        std::ofstream of;
        of.open (fn.c_str ());
        assert (of.is_open ());
        of.precision (6);
        of << "# vtk DataFile Version 2.0" << std::endl;
        of << "generated using DataSink::writeGeoToVtk" << std::endl;
        of << "ASCII" << std::endl << std::endl;
        of << "DATASET UNSTRUCTURED_GRID" << std::endl;
        of << "POINTS " << numpoints_global_m << " float" << std::endl;
        for (int i = 0; i < numpoints_global_m; i++)
            of << geo3Dcoords_m[i](0) << " "
               << geo3Dcoords_m[i](1) << " "
               << geo3Dcoords_m[i](2) << std::endl;
        of << std::endl;

        of << "CELLS "
           << numbfaces_global_m << " "
           << 4 * numbfaces_global_m << std::endl;
        for (int i = 0; i < numbfaces_global_m; i++)
            of << "3 "
               << allbfaces_m[4 * i + 1] << " "
               << allbfaces_m[4 * i + 2] << " "
               << allbfaces_m[4 * i + 3] << std::endl;
        of << "CELL_TYPES " << numbfaces_global_m << std::endl;
        for (int i = 0; i < numbfaces_global_m; i++)
            of << "5" << std::endl;
        of << "CELL_DATA " << numbfaces_global_m << std::endl;
        of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
        of << "LOOKUP_TABLE " << "default" << std::endl;
        for (int i = 0; i < numbfaces_global_m; i++)
            of << (float)(i) << std::endl;
        of << std::endl;

    }
}

Inform& BoundaryGeometry::printInfo (Inform& os) const {
    os << "* *************Boundary Geometry Info*********************************************** " << endl;
    os << "* GEOMETRY                   " << getOpalName () << '\n'
       << "* FGEOM                      " << Attributes::getString (itsAttr[FGEOM]) << '\n'
       << "* TOPO                       " << Attributes::getString (itsAttr[TOPO]) << '\n'
       << "* LENGHT                     " << Attributes::getReal (itsAttr[LENGHT]) << '\n'
       << "* S                          " << Attributes::getReal (itsAttr[S]) << '\n'
       << "* A                          " << Attributes::getReal (itsAttr[A]) << '\n'
       << "* B                          " << Attributes::getReal (itsAttr[B]) << '\n';
    if (getTopology () == string ("BOXCORNER")) {
        os << "* C                          " << Attributes::getReal (itsAttr[C]) << '\n'
           << "* L1                         " << Attributes::getReal (itsAttr[L1]) << '\n'
           << "* L1                         " << Attributes::getReal (itsAttr[L2]) << '\n';
    }
    os << "* Total triangle num         " << numbfaces_global_m << '\n'
       << "* Total points num           " << numpoints_global_m << '\n'
       << "* Triangle side(m)   Max=    " << triangle_max_m << '\n'
       << "*                    Min=    " << triangle_min_m << '\n'
       << "* Geometry bounds(m) Max=    " << maxcoords_m << '\n'
       << "*                    Min=    " << mincoords_m << '\n'
       << "* Geometry length(m)         " << len_m << '\n'
       << "* Boundary box grid num      " << nr_m << '\n'
       << "* Boundary box size(m)       " << hr_m << '\n'
       << "* Size of boundary index set " << boundary_ids_m.size () << '\n'
       << "* Number of all boxes        " << nr_m (0) * nr_m (1) * nr_m (2) << '\n'
       << "* Aligned Triangle Number    " << alignedT_m.size () << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}
