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
    Tribarycent_m (NULL),
    TriPrPartloss_m (NULL),
    TriSePartloss_m (NULL),
    TriFEPartloss_m (NULL),
    allbfaces_m (NULL) {
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


void BoundaryGeometry::computeGeometryInterval (void) {

    mincoords_m = get_min_extend (geo3Dcoords_m);
    maxcoords_m = get_max_extend (geo3Dcoords_m);
    len_m = maxcoords_m - mincoords_m;

    /*
      Calculate the maximum dimension of triangles. This value will be used to
      define the cubic box size
    */

    triangle_max_m = 0.0;
    triangle_min_m = 0.01;
    for (int i = 0; i < numbfaces_global_m; i++) {
        // compute length of longest edge
        Vector_t x1 = geo3Dcoords_m[allbfaces_m[4 * i + 1]];
        Vector_t x2 = geo3Dcoords_m[allbfaces_m[4 * i + 2]];
        Vector_t x3 = geo3Dcoords_m[allbfaces_m[4 * i + 3]];
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
        if (triangle_max_m < max) triangle_max_m = max;
        if (triangle_min_m > max) triangle_min_m = max;
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
    nr_m (0) = (int)floor (len_m (0) / triangle_max_m * 5-0);
    nr_m (1) = (int)floor (len_m (1) / triangle_max_m * 5.0);
    nr_m (2) = (int)floor (len_m (2) / triangle_max_m * 10.0);

    hr_m = len_m / nr_m;
    out_m = maxcoords_m + hr_m;
    *gmsg << "*  Geometry interval built done." << endl;
}

void BoundaryGeometry::initialize () {

    h5_int64_t rc;

    *gmsg << "* Iniitializing Boundary Geometry ... ..." << endl;
    IpplTimings::startTimer (TPreProc_m);

    double xyzscale = getXYZScale ();
    *gmsg << "* Scale all points of the geometry by " << xyzscale << endl;

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
        point_coords[i * 3]   = P[0] * xyzscale;
        point_coords[i * 3 + 1] = P[1] * xyzscale;
        point_coords[i * 3 + 2] = P[2] * xyzscale;
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
    delete point_coords;
    *gmsg << "*  Vertex built done." << endl;

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
    *gmsg << "*  Triangle barycent built done." << endl;

    computeGeometryInterval ();
    makeBoundaryIndexSet ();
    makeTriNormal ();

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

    const Vector_t temp1 = r; //particle position in timestep n;
    const Vector_t temp = r + (c * betaP * v * dt); //particle position in tstep n+1;
    double rI = 0.0;
    Vector_t Isc = out_m;

    IpplTimings::startTimer (TPInside_m);

    /* test if particle position in timestep n is inside the cubic bounding box.
       If true, do the following tests */
    int flg = 0;
    int id;
    if (isInGeometry (temp1)) {
        id = map_point2id (temp1);
    } else if (isInGeometry (temp)) {
        id = map_point2id (temp);
    } else {
        // cannot use goto out here (compiler complains)
        IpplTimings::stopTimer (TPInside_m);
        return -1;
    }

    /* Build an array containing the IDs of 27(3*3*3) cubic boxes. The ID
       of the cubic box which contains the particle, is in the center of
       all this boxes. Just for the situation that even the line segment
       cross more than one cubic box.*/
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

    /* Test all the 27 cubic boxes to find if the line segment has
       intersection with the triangles in those cubic boxes. */
    for (int k = 0; k < 27; k++) {
        std::map< size_t, std::set<size_t> >::iterator It;
        It = CubicLookupTable_m.find (idc[k]);
        if (It == CubicLookupTable_m.end ())
            continue; // not a boundary box

        // for each triangle in this boundary box
        std::set<size_t> ::iterator faceIt;
        for (faceIt = (*It).second.begin ();
             faceIt != (*It).second.end ();
             faceIt++) {
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

/*
  P R I V A T E   F U N C T I O N S
 */

static void
write_bbox_mesh (
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

int64_t fcmp (
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

bool is_in_bbox (Vector_t point, Vector_t min, Vector_t max) {
    if (fcmp (point [0], min [0], 10) < 0) return false;
    if (fcmp (point [1], min [1], 10) < 0) return false;
    if (fcmp (point [2], min [2], 10) < 0) return false;
    if (fcmp (point [0], max [0], 10) > 0) return false;
    if (fcmp (point [1], max [1], 10) > 0) return false;
    if (fcmp (point [2], max [2], 10) > 0) return false;
    return true;
}


/*
   Make the boundary set by using
   * triangle vertices
   * bounding box vertices
   * several points in triangle central lines
   * several points in triangle edges

*/
void BoundaryGeometry::makeBoundaryIndexSet () {
    for (int i = 0; i < numbfaces_global_m; i++) {
        std::vector<Vector_t> coords;
        /* Discretize the three central lines and the three triangle edges to
           200 segments to get a more complete boundary index set. */
        Vector_t c1 = 0.5 * (getVertexCoord (i, 2) + getVertexCoord (i, 3)) - getVertexCoord (i, 1);
        Vector_t c2 = 0.5 * (getVertexCoord (i, 3) + getVertexCoord (i, 1)) - getVertexCoord (i, 2);
        Vector_t c3 = 0.5 * (getVertexCoord (i, 1) + getVertexCoord (i, 2)) - getVertexCoord (i, 3);

        Vector_t e1 = getVertexCoord (i, 2) - getVertexCoord (i, 1);
        Vector_t e2 = getVertexCoord (i, 3) - getVertexCoord (i, 2);
        Vector_t e3 = getVertexCoord (i, 1) - getVertexCoord (i, 3);
        for (int j = 1; j <= 200; j++) {
            // discretize the three central lines.
            coords.push_back (getVertexCoord (i, 1) + 0.005 * j * c1);
            coords.push_back (getVertexCoord (i, 2) + 0.005 * j * c2);
            coords.push_back (getVertexCoord (i, 3) + 0.005 * j * c3);
            // discretize the three  triangle edges:
            coords.push_back (getVertexCoord (i, 1) + 0.005 * j * e1);
            coords.push_back (getVertexCoord (i, 2) + 0.005 * j * e2);
            coords.push_back (getVertexCoord (i, 3) + 0.005 * j * e3);
        }

        coords.push_back (Tribarycent_m[i]);

        std::vector<Vector_t> ret = getMinBBoxOfTriangle (i);

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
        min -= hr_m;
        max += hr_m;
        Vector_t P = Vector_t (min[0], min[1], min[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (min[0], min[1], min[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (min[0], min[1], max[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (min[0], max[1], max[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (min[0], max[1], min[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (max[0], min[1], min[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (max[0], min[1], max[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (max[0], max[1], max[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);
        P = Vector_t (max[0], max[1], min[2]);
        if (is_in_bbox (P, mincoords_m, maxcoords_m)) coords.push_back (P);

        std::vector<Vector_t>::iterator point;
        for (point = coords.begin (); point != coords.end (); point++) {
            /*
            if (!is_in_bbox (*point, mincoords_m, maxcoords_m))
                continue;
            */
            int id = map_point2id (*point);
            assert (id > 0);
            // insert ID to std::set! If ID is alread in the set, this is a NOP.
            boundary_ids_m.insert (id);
            
            std::map< size_t, std::set<size_t> >::iterator It;
            It =  CubicLookupTable_m.find (id);
            if (It == CubicLookupTable_m.end ()) {
                std::set<size_t> tmp;
                tmp.insert (i);
                CubicLookupTable_m.insert (std::pair <size_t, std::set<size_t> > (id, tmp));
            } else
                (*It).second.insert (i);
        }
    }
    if(Ippl::myNode() == 0) {
        write_bbox_mesh (boundary_ids_m, hr_m, nr_m, mincoords_m);
    }
    *gmsg << "*  Boundary index set built done." << endl;
}

Vector_t BoundaryGeometry::LineInsTri (
    Vector_t x0,
    Vector_t x1,
    size_t i) {
    Vector_t lseg = x1 - x0;
    Vector_t t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
    Vector_t t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
    Vector_t t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
    Vector_t u = t1 - t0;
    Vector_t v = t2 - t0;
    Vector_t lt = t0 - x0;
    Vector_t n;
    n[0] = u[1] * v[2] - v[1] * u[2];
    n[1] = u[2] * v[0] - v[2] * u[0];
    n[2] = u[0] * v[1] - v[0] * u[1];
    
    if (fabs (dot (n, lseg)) < 1.0e-10) {
        return x1;              // Triangle parallel to line segment
    } else {
        double rI = dot(n, lt) / dot(n, lseg);
        Vector_t ItSec = x0 + rI * lseg;
        Vector_t w = ItSec - t0;
        if((rI < 0) || (rI > 1)) {
            return x1;          // Intesect is on the extended line
        } else {
            double tmp1 = dot (u, v);
            double tmp2 = dot (w, v);
            double tmp3 = dot (u, w);
            double tmp4 = dot (u, u);
            double tmp5 = dot (v, v);
            double tmp6 = tmp1 * tmp1 - tmp4 * tmp5;
            double sI = (tmp1 * tmp2 - tmp5 * tmp3) / tmp6;
            double tI = (tmp1 * tmp3 - tmp4 * tmp2) / tmp6;
            if((sI >= 0) && (tI >= 0) && ((sI + tI) <= 1)) {
                return ItSec;
            } else {
                return x1;      // Intesect is on the extended plane
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
