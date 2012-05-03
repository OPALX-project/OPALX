#ifndef OPAL_BOUNDARY_GEOMETRY_HH
#define OPAL_BOUNDARY_GEOMETRY_HH

class OpalBeamline;
class ElementBase;

#include "Algorithms/PartBunch.h"
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "Distribution/ranlib.h"
#include "Structure/SecondaryEmissionPhysics.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cmath>
#include <assert.h>

namespace BGphysics {
    enum TPHYACTION {
        Nop = 0x01,                 // tringle is transparent to particle like beam window
        Absorption = 0x02,          // triangle has no field emission and secondary emission
        FNEmission = 0x04,          // triangle has field emission
        SecondaryEmission = 0x08    // trangle has secondary emission
    };
}

class BoundaryGeometry : public Definition {

public:
    BoundaryGeometry();
    virtual ~BoundaryGeometry();

    virtual bool canReplaceBy (Object* object);

    virtual BoundaryGeometry* clone (const std::string& name);

    // Check the GEOMETRY data.
    virtual void execute ();

    // Find named GEOMETRY.
    static BoundaryGeometry* find (const std::string& name);

    // Update the GEOMETRY data.
    virtual void update ();

    void updateElement (ElementBase* element);

    void initialize ();
    void createParticlesOnSurface (
        size_t n, double darkinward,
        OpalBeamline& itsOpalBeamline,
        PartBunch& itsBunch);
    void createPriPart (
        size_t n, double darkinward,
        OpalBeamline& itsOpalBeamline,
        PartBunch* itsBunch);

    inline string getFilename () const {
        return (string) Attributes::getString (itsAttr[FGEOM]);
    }

    inline string getTopology () const {
        return (string) Attributes::getString (itsAttr[TOPO]);
    }

    inline string getDistribution () {
        return (string) Attributes::getString (itsAttr[DISTR]);
    }

    inline std::vector<string> getDistributionArray () {
        return Attributes::getStringArray (itsAttr[DISTRS]);
    }

    inline double getZshift () {
        return (double)(Attributes::getReal (itsAttr[ZSHIFT]));
    }

    inline double getXYZScale () {
        return (double)(Attributes::getReal (itsAttr[XYZSCALE]));
    }

    inline size_t getN () {
        return partsr_m.size ();
    }

    inline Vector_t getCooridinate (size_t i) {
        return partsr_m[i];
    }
    inline void clearCooridinateArray () {
        return partsr_m.clear ();
    }
    inline Vector_t getMomenta (size_t i) {
        return partsp_m[i];
    }

    inline void clearMomentaArray () {
        return partsp_m.clear ();
    }

    bool isInside (Vector_t x);
    std::vector<Vector_t>  GridIntersection (Vector_t x0, Vector_t x1);
    Inform& printInfo (Inform& os) const;

    inline double getA() {
        return (double)Attributes::getReal(itsAttr[A]);
    }

    inline double getB() {
        return (double)Attributes::getReal(itsAttr[B]);
    }

    inline double getC() {
        return (double)Attributes::getReal(itsAttr[C]);
    }

    inline double getS() {
        return (double)Attributes::getReal(itsAttr[S]);
    }

    inline double getLenght() {
        return (double)Attributes::getReal(itsAttr[LENGHT]);
    }

    inline double getL1() {
        return (double)Attributes::getReal(itsAttr[L1]);
    }

    inline double getL2() {
        return (double)Attributes::getReal(itsAttr[L2]);
    }

    void makeBoundaryIndexSet ();

    int PartInside (const Vector_t r, const Vector_t v, const double dt,
                    int Parttype, const double Qloss, Vector_t& intecoords,
                    int& triId, double& Energy);

    // non secondary emission version.
    int doBGphysics (const Vector_t& intecoords, const int& triId);

    // call Furman-Pivi's model
    int doBGphysics (const Vector_t& intecoords, const int& triId,
                     const double& incEnergy, const double& incQ,
                     const Vector_t& incMomentum, PartBunch* itsBunch,
                     double& seyNum);

    // call Vaughan's model
    int doBGphysics (const Vector_t& intecoords, const int& triId,
                     const double& incEnergy, const double& incQ,
                     const Vector_t& incMomentum, PartBunch* itsBunch,
                     double& seyNum,
                     const int& para_null); //



    inline void setNEmissionMode (bool nEmissionMode) {
        nEmissionMode_m = nEmissionMode;
    }

    inline void setWorkFunction (double workFunction) {
        workFunction_m = workFunction;
    }

    inline void setFieldEnhancement (double fieldEnhancement) {
        fieldEnhancement_m = fieldEnhancement;
    }

    inline void setMaxFN (size_t maxFNemission) {
        maxFNemission_m = maxFNemission;
    }

    inline void setFNTreshold (double fieldFNthreshold) {
        fieldFNthreshold_m = - 1.0e6 * fieldFNthreshold;
    }

    inline void setFNParameterA (double parameterFNA) {
        parameterFNA_m = parameterFNA;
    }

    inline void setFNParameterB (double parameterFNB) {
        parameterFNB_m = parameterFNB;
    }

    inline void setFNParameterY (double parameterFNY) {
        parameterFNY_m = parameterFNY;
    }

    inline void setFNParameterVYZe (double parameterFNVYZe) {
        parameterFNVYZe_m = parameterFNVYZe;
    }

    inline void setFNParameterVYSe (double parameterFNVYSe) {
        parameterFNVYSe_m = parameterFNVYSe;
    }

    inline void setBoundaryMatType (int BoundaryMatType) {
        seBoundaryMatType_m = BoundaryMatType;
    }

    inline void setEInitThreshold (double einitthreshold) {
        eInitThreshold_m = 1.0e6 * einitthreshold;
    }

    // return sey_0 in Vaughan's model
    inline void setvSeyZero (double vSeyZero) {
        vSeyZero_m = vSeyZero;
    }

    // set the energy related to sey_0 in Vaughan's model
    inline void setvEZero (double vEZero) {
        vEzero_m = vEZero;
    }

    // set sey max in Vaughan's model
    inline void setvSeyMax (double vSeyMax) {
        vSeyMax_m = vSeyMax;
    }

    // return Emax in Vaughan's model
    inline void setvEmax (double vEmax) {
        vEmax_m = vEmax;
    }

    // return fitting parameter denotes the roughness of surface for
    // impact energy in Vaughan's model
    inline void setvKenergy (double vKenergy) {
        vKenergy_m = vKenergy;
    }
    
    // return fitting parameter denotes the roughness of surface for impact
    // angle in Vaughan's model
    inline void setvKtheta (double vKtheta) {
        vKtheta_m = vKtheta;
    }

    // return thermal velocity Maxwellian distribution of secondaries
    // in Vaughan's model
    inline void setvVThermal (double vVThermal) {
        vVThermal_m = vVThermal;
    }

    inline void setVw (double ppVw) {
        ppVw_m = ppVw;
    }

    size_t doFNemission (OpalBeamline& itsOpalBeamline, PartBunch* itsBunch, const double t);

    bool nEmissionMode_m;
    double eInitThreshold_m;

    // energy related to sey_0 in Vaughan's model
    double vSeyZero_m;

    // sey_0 in Vaughan's model
    double vEzero_m;

     // sey max in Vaughan's model
    double vSeyMax_m;

    // Emax in Vaughan's model
    double vEmax_m;

    // fitting parameter denotes the roughness of surface
    // for impact energy in Vaughan's model
    double vKenergy_m;

    // fitting parameter denotes the roughness of surface
    // for impact angle in Vaughan's model
    double vKtheta_m;

    // return thermal velocity of Maxwellian distribution
    // of secondaries in Vaughan's model
    double vVThermal_m;

    // denotes the velocity scalar for Parallel plate benchmark.
    double ppVw_m;

    // stores the user defined material type for secondary emission model.
    int seBoundaryMatType_m;

    // stores the user defined work function Fowler-Nordheim model.
    double workFunction_m;

    // stores the user defined field factor for Fowler-Nordheim model.
    double fieldEnhancement_m;

    // stores the user defined maximum emitted number per triangle for
    // Fowler-model.
    size_t maxFNemission_m;

    // stores the user defined lower threshold electric field for
    // Fowler-Nordheim model.
    double fieldFNthreshold_m;

    /**
       @param parameterFNA_m stores the user defined parameter A for
       Fowler-Nordheim model.Default value:\f$1.54\times 10^{-6}\f$.
    */
    double parameterFNA_m;

    /**
       @param parameterFNB_m stores the user defined parameter B 
       Fowler-Nordheim model. Default value:\f$6.83\times 10^9\f$.
    */
    double parameterFNB_m;

    /**
       @param parameterFNY_m stores the user defined parameter for
       Fowler-Nordheim parameter y.
       Default value:\f$3.795\times 10^{-5}\f$.
    */
    double parameterFNY_m;

    /**
       @param parameterFNVYZe_m stores the user defined parameter for zero
       order fit constant for v(y) of Fowler-Nordheim model.
       Default value:\f$0.9632\f$.
    */
    double parameterFNVYZe_m;
    
    /**
       @param parameterFNVYSe_m stores the user defined parameter for second
       order fit constant for v(y) of Fowler-Nordheim model.
       Default value:\f$1.065\f$.
    */
    double parameterFNVYSe_m;

    /**
       Return number of boundary faces.
    */
    inline int getNumBFaces () {
            return numbfaces_global_m;
    }

    /**
     @param allbfaces_m stores the IDs of triangle vertex, which can be used
     to reference the  coordinates of vertex
    */
    int* allbfaces_m;

    /**
       @param numpoints_global_m stores the number of geometry points to be used
       to represent the geometry
       caution: not only the surface points are included
    */
    int numpoints_global_m;

    /**
       @param  Tribarycent_m store the coordinates of barycentric points of
       triangles, The Id number is the same with triangle Id.
    */
    Vector_t* Tribarycent_m;

    /**
       @param  TriNormal_m store the oriented normal vector of triangles, The Id
       number is the same with triangle Id.
    */
    std::vector<Vector_t>    TriNormal_m;

    /**
       @param TriPrPartloss_m store the number of primary particles hitting the
       Id th triangle. The Id number is the same with triangle Id(not vertex ID).
    */
    double* TriPrPartloss_m;

    /**
       @param TriSePartloss_m store the number of secondary particles hitting the
       Id th triangle. The Id number is the same with triangle Id(not vertex ID).
    */
    double* TriSePartloss_m;

    /**
       @param TriFEPartloss_m store the number of field emission/darkcurrent
       particles hitting the Id th triangle. The Id number is the same with
       triangle Id(not vertex ID).
    */
    double* TriFEPartloss_m;

    /**
       @param TriPrPartlossZ_m is a counter of lost primary particle number in
       each Z intervals.
    */
    double* TriPrPartlossZ_m;

    /**
       @param TriSePartlossZ_m is a counter of lost secondary particle number in
       each Z intervals.
    */
    double* TriSePartlossZ_m;

    /**
       @param TriFEPartlossZ_m is a counter of lost field emission and dark
       current particle number in each Z intervals.
    */
    double* TriFEPartlossZ_m;    //counter along Z for histogram

    /**
       @param geo3Dcoords_m store the geometry coordinates in a STL Vector
       The ID of geo3Dcoords_m is equal to points' ID, and can be referenced by
       triangle to get vertex coordinates
    */
    std::vector<Vector_t> geo3Dcoords_m;

    /**
       @param Triarea_m store area of triangles in a STL Vector
       The ID of Triarea_m is equal to triangles' ID, and can be referenced by
       triangle to get its area
    */
    std::vector<double> Triarea_m;

    /**
       @param TriBGphysicstag_m store the tags of each boundary triangle for
       proper physics action.
    */
    std::vector<short> TriBGphysicstag_m;

    std::string h5FileName_m;

    /**
       Used to determine whether a particle hits boundary.
       @param x stands for the coordinates of the testing particle

       Not sufficient for particle hitting the triangle area.
    */
    bool isInGeometry (Vector_t x) {
        return boundary_ids_m.find (f (x)) != boundary_ids_m.end ();
    }

    /**
       Return the hr_m.
    */
    Vector_t gethr () {
        return hr_m;
    }
    /**
       Return the nr_m.
     */
    Vektor<int, 3> getnr () {
        return nr_m;
    }

    /**
       Return the mincoords_m.
     */
    Vector_t getmincoords () {
        return mincoords_m;
    }
    /**
       Return the maxcoords_m.
    */
    Vector_t getmaxcoords () {
        return maxcoords_m;
    }

    void writeGeomToVtk (string fn);


private:
    double triangle_max_m;
    double triangle_min_m;
    int numbfaces_global_m;             // number of boundary triangles
    double xyzscale_m;                  //  scale of x,y and z coordinates

    std::set<size_t> boundary_ids_m;    // boundary triangle IDs
    Vector_t len_m;                     // length of geometry in 3D Cartesian coordinates.
    Vector_t hr_m;                      // length of cubic box
    Vektor<int, 3> nr_m;                // number of intervals of geometry in X,Y,Z direction
    Vector_t mincoords_m;               // minimum of geometry coordinate.
    Vector_t maxcoords_m;               // maximum of geometry coordinate.
    std::vector<Vector_t> partsp_m;     // particle momenta
    std::vector<Vector_t> partsr_m;     // particle positions

    std::map< size_t, std::vector<size_t> >
            CubicLookupTable_m;         // Maps boundary box ID to included triangles
    std::set<size_t> isOriented_m;      //  IDs of oriented triangles.
    std::map< size_t, std::vector<size_t> >
            triangleLookupTable_m;      // map vertex ID to triangles with this vertex
    std::vector<size_t> NaliT_m;        // IDs of to be oriented triangles
    std::vector<size_t> alignedT_m;     // IDs of oriented triangles
    Vector_t out_m;                     // a point outside the domain


    IpplTimings::TimerRef TPInside_m;   // timers for profiling
    IpplTimings::TimerRef TPreProc_m;
    IpplTimings::TimerRef TRayTrace_m;
    IpplTimings::TimerRef Tinward_m;

    SecondaryEmissionPhysics sec_phys_m;

    // Not implemented.
    BoundaryGeometry(const BoundaryGeometry&);
    void operator= (const BoundaryGeometry&);

    // Clone constructor.
    BoundaryGeometry(const std::string& name, BoundaryGeometry* parent);

    inline Vector_t getVertexCoord (int face_id, int vertex_id) {
        return geo3Dcoords_m[allbfaces_m[4 * face_id + vertex_id]];
    }

    /*
      Recursively get inward normal of triangle with id=idx by switching the order
      of vertex w.r.t the vertex order of triangle with id = caller.
    */
    void orientAllTriangles (size_t caller, size_t idx) {
        orientTriangle (idx, caller);
        isOriented_m.insert (idx);
        std::vector<size_t> neighbours = findNeighbours (idx);
        for (unsigned int i = 0; i < neighbours.size (); i++) {
            if (isOriented_m.find (neighbours[i]) == isOriented_m.end ())
                orientAllTriangles (idx, neighbours[i]);
        }
    }

    std::vector<size_t> findNeighbours (size_t idx) {
        std::vector<size_t> ret, ret1;
        std::set<size_t>temp;
        std::vector<size_t>::iterator retIt;
        std::map< size_t, std::vector<size_t> >::iterator it;
        for (int i = 1; i <= 3; i++) {
            it = triangleLookupTable_m.find (triangleCorner (idx, i));
            ret.insert (ret.end (), (*it).second.begin (), (*it).second.end ());
        }
        for (retIt = ret.begin (); retIt != ret.end (); retIt++) {
            if (*retIt != idx) {
                if (temp.find (*retIt) == temp.end ()) {
                    temp.insert (*retIt);
                } else
                    ret1.push_back (*retIt);
            }
        }
        PAssert (ret1.size () != 0);
        return ret1;
    }

    inline size_t triangleCorner (size_t idx, size_t corner) const {
        return allbfaces_m[4 * idx + corner];
    }

    void orientTriangle (size_t idx, size_t caller) {
        std::vector<size_t> id, ic;
        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                if (triangleCorner (idx, j) == triangleCorner (caller, i)) {
                    id.push_back (j);
                    ic.push_back (i);
                }
            }
        }
        if ((ic[1] - ic[0]) == 1) {
            int idtmp = id[1] - id[0];
            if (idtmp == 1) {
                alignedT_m.push_back (idx);
                allbfaces_m[4 * idx + id[0]] = allbfaces_m[4 * caller + ic[1]];
                allbfaces_m[4 * idx + id[1]] = allbfaces_m[4 * caller + ic[0]];
            }
            if (idtmp == - 1) {
                NaliT_m.push_back (idx);
            }
            if (idtmp == 2) {
                NaliT_m.push_back (idx);
            }
            if (idtmp == - 2) {
                alignedT_m.push_back (idx);
                allbfaces_m[4 * idx + id[0]] = allbfaces_m[4 * caller + ic[1]];
                allbfaces_m[4 * idx + id[1]] = allbfaces_m[4 * caller + ic[0]];
            }
        }
        if ((ic[1] - ic[0]) == 2) {
            int idtmp = id[1] - id[0];
            if (idtmp == - 1) {
                alignedT_m.push_back (idx);
                allbfaces_m[4 * idx + id[0]] = allbfaces_m[4 * caller + ic[1]];
                allbfaces_m[4 * idx + id[1]] = allbfaces_m[4 * caller + ic[0]];
            }
            if (idtmp == 1) {
                NaliT_m.push_back (idx);
            }
            if (idtmp == - 2) {
                NaliT_m.push_back (idx);
            }
            if (idtmp == 2) {
                alignedT_m.push_back (idx);
                allbfaces_m[4 * idx + id[0]] = allbfaces_m[4 * caller + ic[1]];
                allbfaces_m[4 * idx + id[1]] = allbfaces_m[4 * caller + ic[0]];
            }
        }
    }
    /*
      Recursively get inward triangle normal of all surface triangles.

      The basic idea is as follow: get the first triangle's inward normal by
      determine a nearby point is inside or outside boundary geometry (using
      ray-triangle intersection and even/odd intersection number).
      Then use a recursion method to switch the vertex order of adjacent
      triangles. The inward normal is stored in TriNormal_m.
    */
    void makeTriNormal () {
        Vector_t t0, t1, t2;
        t0 = geo3Dcoords_m[allbfaces_m[1]];
        t1 = geo3Dcoords_m[allbfaces_m[2]];
        t2 = geo3Dcoords_m[allbfaces_m[3]];
        Vector_t u = t1 - t0;
        Vector_t v = t2 - t0;
        Vector_t n;
        n[0] = u[1] * v[2] - v[1] * u[2];
        n[1] = u[2] * v[0] - v[2] * u[0];
        n[2] = u[0] * v[1] - v[0] * u[1];
        double nomtmp = sqrt (n (0) * n (0) + n (1) * n (1) + n (2) * n (2));
        n = n / nomtmp;
        TriNormal_m.push_back (n);

        Vector_t mytemp;
        mytemp[0] = t0[0] + TriNormal_m[0](0) * 0.1 * triangle_min_m;
        mytemp[1] = t0[1] + TriNormal_m[0](1) * 0.1 * triangle_min_m;
        mytemp[2] = t0[2] + TriNormal_m[0](2) * 0.1 * triangle_min_m;
        if (!isInside (mytemp)) {
            if (dot (mytemp - t0, TriNormal_m[0]) < 0) {
                // TriNormal_m[i]=TriNormal_m[i];
            } else {
                TriNormal_m[0] = - TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            }
        } else {
            if (dot (mytemp - t0, TriNormal_m[0]) < 0) {
                TriNormal_m[0] = - TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            } else {
                TriNormal_m[0] = TriNormal_m[0];
            }
        }

        // for all triangles find adjacent triangles to each vertex
        for (int i = 0; i < numbfaces_global_m; i++) {
            for (int j = 1; j <= 3; j++) {
                std::map< size_t, std::vector<size_t> >::iterator it;
                it = triangleLookupTable_m.find (allbfaces_m[4 * i + j]);
                if (it == triangleLookupTable_m.end ()) {
                    std::vector <size_t> tmp;
                    tmp.push_back (i);
                    triangleLookupTable_m.insert (
                            std::pair<size_t, std::vector<size_t> > (
                                    allbfaces_m[4 * i + j], tmp));
                } else
                    (*it).second.push_back (i);
            }
        }
        isOriented_m.insert (0);
        std::vector<size_t> neighbours = findNeighbours (0);

        for (unsigned int i = 0; i < neighbours.size (); i++) {
            if (isOriented_m.find (neighbours[i]) == isOriented_m.end ())
                orientAllTriangles (0, neighbours[i]);
        }
        for (int i = 1; i < numbfaces_global_m; i++) {
            Vector_t t0 = geo3Dcoords_m[allbfaces_m[4 * i + 1]];
            Vector_t t1 = geo3Dcoords_m[allbfaces_m[4 * i + 2]];
            Vector_t t2 = geo3Dcoords_m[allbfaces_m[4 * i + 3]];
            Vector_t u = t1 - t0;
            Vector_t v = t2 - t0;
            Vector_t n;
            n[0] = u[1] * v[2] - v[1] * u[2];
            n[1] = u[2] * v[0] - v[2] * u[0];
            n[2] = u[0] * v[1] - v[0] * u[1];
            double nomtmp = sqrt (n (0) * n (0) + n (1) * n (1) + n (2) * n (2));
            if (nomtmp != 0) {
                n = n / nomtmp;
            }
            TriNormal_m.push_back (n);
        }
    }

    Vector_t LineInsTri(
            Vector_t x0,
            Vector_t x1,
            size_t i
        ) {
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

    /*
      Calculate the number of intersects between a line segment and the
      geometry boundary, make sure x1 is outside the geometry.
    */
    std::vector<Vector_t>  PartBoundaryInteNum (Vector_t x0, Vector_t x1) {
        std::vector<Vector_t> Isp;
        for (int i = 0; i < numbfaces_global_m; i++) {
            Vector_t tmp =  LineInsTri(x0, x1, i);
            if (tmp != x1) {
                Isp.push_back (tmp);
            }
        }
        if (Isp.size () == 0)
            Isp.push_back (x1);

        return Isp;
    }

    inline double SQR (double x) {
	    return x * x;
    }

    /*
      Calculate the maximum dimention of triangles. This value will be used to
      define the cubic box size
    */
    void getMaxDimenssion () {
        triangle_max_m = 0.0;
        triangle_min_m = 0.01;
        for (int i = 0; i < numbfaces_global_m; i++) {
            Vector_t x1 = geo3Dcoords_m[allbfaces_m[4 * i + 1]];
            Vector_t x2 = geo3Dcoords_m[allbfaces_m[4 * i + 2]];
            Vector_t x3 = geo3Dcoords_m[allbfaces_m[4 * i + 3]];
            double max = sqrt (
                    SQR (x1[0] - x2[0]) + SQR (x1[1] - x2[1]) + SQR (x1[2] - x2[2]));
            double length_edge2 = sqrt (
                    SQR (x3[0] - x2[0]) + SQR (x3[1] - x2[1]) + SQR (x3[2] - x2[2]));
            double length_edge3 = sqrt (
                    SQR (x3[0] - x1[0]) + SQR (x3[1] - x1[1]) + SQR (x3[2] - x1[2]));
            if (length_edge2 > max) max = length_edge2;
            if (length_edge3 > max) max = length_edge3;

            if (triangle_max_m < max) triangle_max_m = max;
            if (triangle_min_m > max) triangle_min_m = max;
        }
    }

    /*
      Map a 3D coordinate given in x to a cubic box with a unique Id
    */
    int f (Vector_t x) {
        /**
           We know: lenght of the structure len_m
           Mesh size: hr_m
           Number of mesh points nr_m
         */
        size_t ret = 0;
        size_t id_tx, id_ty, id_tz;
        id_tx = floor ((x[0] - mincoords_m[0]) / hr_m[0]);
        id_ty = floor ((x[1] - mincoords_m[1]) / hr_m[1]);
        id_tz = floor ((x[2] - mincoords_m[2]) / hr_m[2]);
        ret = 1 + id_tz * nr_m[0] * nr_m[1] + id_ty * nr_m[0] + id_tx;
        return ret;
    }

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
    void FindIntersection (
            const Vector_t& x0,
            const Vector_t& x1,
            const size_t& i,
            double& rI, Vector_t& Isc
        ) {
        IpplTimings::startTimer (TRayTrace_m);

        const Vector_t lseg = x1 - x0; // length and direction of line segment;
        const Vector_t t0 = geo3Dcoords_m[allbfaces_m[4 * i + 1]];
        const Vector_t t1 = geo3Dcoords_m[allbfaces_m[4 * i + 2]];
        const Vector_t t2 = geo3Dcoords_m[allbfaces_m[4 * i + 3]];
        const Vector_t u = t1 - t0; // side 1 of triangle;
        const Vector_t v = t2 - t0; // side 2 of triangle;
        const Vector_t lt = t0 - x0;
        const Vector_t n = TriNormal_m[i];

        const double dotLT = dot (n, lseg);
        if (dotLT == 0) {
            if ((x0 == Tribarycent_m[i]) && (x0 != x1)) {
                /*
                  Some initialized particles have momenta parallel to its
                  triangle normal, this kind of particles will lose
                  directly
                */
                Isc = Tribarycent_m[i];
            }
        } else {
            // find intersection position w.r.t x0 and the unit is (x1-x0);
            rI = dot (n, lt) / dotLT;

            // find the coordinate of
            // intersection.
            const Vector_t ItSec = x0 + rI * lseg;
            // find if the intersection is inside the triangle.
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
                Isc = ItSec;
            }
        }
        IpplTimings::stopTimer (TRayTrace_m);
    }

    /*
       Get the smallest bounding box of triangle given by ID. The
       smallest bounding box is used to make sure that all part of a triangle
       is known by boundary bounding box.
    */
    std::vector<Vector_t> SetMinMaxBound (size_t id) {
        Vector_t min = geo3Dcoords_m[allbfaces_m[4 * id + 1]];
        Vector_t max = min;
        for (int i = 2; i <= 3; i++) {
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](0) < min[0])
                min[0] = geo3Dcoords_m[allbfaces_m[4 * id + i]](0);
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](1) < min[1])
                min[1] = geo3Dcoords_m[allbfaces_m[4 * id + i]](1);
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](2) < min[2])
                min[2] = geo3Dcoords_m[allbfaces_m[4 * id + i]](2);
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](0) > max[0])
                max[0] = geo3Dcoords_m[allbfaces_m[4 * id + i]](0);
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](1) > max[1])
                max[1] = geo3Dcoords_m[allbfaces_m[4 * id + i]](1);
            if (geo3Dcoords_m[allbfaces_m[4 * id + i]](2) > max[2])
                max[2] = geo3Dcoords_m[allbfaces_m[4 * id + i]](2);
        }
        /*
          add some margins to make sure that boundary bounding box has both
          positive and negtive margins w.r.t boundary. Just make sure that the
          size of bounding box for triangle is smaller than the boundary bounding
          box.
        */
        Vector_t temp = 1.0 * hr_m;
        std::vector<Vector_t> ret;
        Vector_t min_ext = min - temp;
        Vector_t max_ext = max + temp;
        ret.push_back (min);
        ret.push_back (max);
        ret.push_back (min_ext);
        ret.push_back (max_ext);
        PAssert (ret.size () != 0);
        return ret;
    }

    // Calculate the area of triangle given by id.
    double TriangleArea (int id) {
        Vector_t AB = geo3Dcoords_m[allbfaces_m[4 * id + 2]]
                - geo3Dcoords_m[allbfaces_m[4 * id + 1]];
        Vector_t AC = geo3Dcoords_m[allbfaces_m[4 * id + 3]]
                - geo3Dcoords_m[allbfaces_m[4 * id + 1]];

        return(0.5 * sqrt (dot (AB, AB) * dot (AC, AC) - dot (AB, AC) * dot (AB, AC)));
    }

    /*
      We define some tags in namespace BGphysics for each surface triangle to
      identify the physical reactions for each triangle when amplitude of
      electrostatic field exceeds some threshold or particles incident the surface.
    */
    void setBGphysicstag () {
        for (int i = 0; i < numbfaces_global_m; i++) {
            TriBGphysicstag_m.push_back (
                BGphysics::Absorption
                | BGphysics::FNEmission
                | BGphysics::SecondaryEmission);
        }
    }

    enum {
        FGEOM,    // file holding the geometry
        LENGHT,   // length of elliptic tube or boxcorner
        S,        // start of the geometry
        L1,       // in case of BOXCORNER first part of geometry with hight B
        L2,       // in case of BOXCORNER second part of geometry with hight B-C
        A,        // major semi-axis of elliptic tube
        B,        // minor semi-axis of ellitpic tube
        C,        // in case of BOXCORNER hight of corner
        TOPO,     // BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written
        DISTR,    // Add distribution to generate physics model on the surface
        DISTRS,   // Add distribution array to generate physics model on the surface
        ZSHIFT,   // Shift in z direction
        XYZSCALE,  // Multiplicative scaling factor for coordinates
        SIZE
    };
};

inline Inform &operator<< (Inform& os, const BoundaryGeometry& b) {
    return b.printInfo (os);
}


#endif
