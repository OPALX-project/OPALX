#ifndef OPAL_BOUNDARY_GEOMETRY_HH
#define OPAL_BOUNDARY_GEOMETRY_HH

// ------------------------------------------------------------------------
// $RCSfile: Geometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------
class OpalBeamline;
class ElementBase;

#include "Algorithms/PartBunch.h"
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "Ippl.h"

/**
 * @param TPInside_m, @param TPreProc_m, @param TRayTrace_m are timers for out put the time consumptions.
 */
static IpplTimings::TimerRef TPInside_m = IpplTimings::getTimer("Particle Inside");
static IpplTimings::TimerRef TPreProc_m = IpplTimings::getTimer("Pre Processing");
static IpplTimings::TimerRef TRayTrace_m = IpplTimings::getTimer("Ray tracing");

// Class BoundaryGeometry
// ------------------------------------------------------------------------
/// The GEOMETRY definition.
//  A GEOMETRY definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.
//
//  i.e:
//  G1: Geometry, FILE="input.h5"
//  G2: Geometry, L=1.0, A=0.0025, B=0.0001
//

class BoundaryGeometry: public Definition {

public:

    /// Exemplar constructor.
    BoundaryGeometry();

    virtual ~BoundaryGeometry();

    /// Test if replacement is allowed.
    //  Can replace only by another GEOMETRY.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual BoundaryGeometry *clone(const string &name);

    /// Check the GEOMETRY data.
    virtual void execute();

    /// Find named GEOMETRY.
    static BoundaryGeometry *find(const string &name);

    /// Update the GEOMETRY data.
    virtual void update();

    void updateElement(ElementBase *element);

    void initialize();
    void createParticlesOnSurface(size_t n, double darkinward, OpalBeamline &itsOpalBeamline, PartBunch &itsBunch);


    // inline size_t getN() { return geo3Dcoords_m.size(); }
    inline size_t getN() { return partsr_m.size();}
    inline Vector_t getCooridinate(size_t i) { return partsr_m[i]; }
    inline void clearCooridinateArray() { return partsr_m.clear(); }


    //    vector<Vector_t> PartBoundaryInteNum(Vector_t x0,Vector_t x1);
    int isInside(Vector_t x);

    Inform &print(Inform &os) const;
    string getTopology();
    string getFilename();

    double getA();
    double getB();
    double getS();

    string getDistribution();

    void makeBoundaryIndexSet();
    int PartInside(Vector_t r, Vector_t v, double dt, short Parttype, Vector_t &intecoords, int &triId);

    void setBGphysicstag();
    int doBGphysics(Vector_t intecoords, int triId);
    void setWorkFunction(double workFunction);
    void setFieldEnhancement(double fieldEnhancement);
    void setMaxFN(size_t maxFNemission);
    void callFNemission(OpalBeamline &itsOpalBeamline, PartBunch *itsBunch, const double t);
    double workFunction_m;
    double fieldEnhancement_m;
    size_t  maxFNemission_m;
    /**
     * @param numbfaces_global_m stores the number of triangles used to represent the geometry
     */
    int     numbfaces_global_m;
    /**
     * @param allbfaces_m stores the IDs of triangle vertex, which can be used to reference the  coordinates of vertex
     */
    int    *allbfaces_m;
    /**
     * @param numpoints_global_m stores the number of geometry points to be used to represent the geometry
     * caution: not only the surface points are included
     */
    int     numpoints_global_m;
    /**
     * @param  Tribarycent_m store the coordinates of barycentric points of triangles, The Id number is the same with triangle Id.
     */
    Vector_t    *Tribarycent_m;
    /**
     * @param  TriNormal_m store the oriented normal vector of triangles, The Id number is the same with triangle Id.
     */
    vector<Vector_t>    TriNormal_m;
    /**
     * @param TriPrPartloss_m store the number of primary particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    int    *TriPrPartloss_m;
    /**
     * @param TriSePartloss_m store the number of secondary particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    int    *TriSePartloss_m;
    /**
     * @param TriPrPartlossZ_m is a counter of lost particle number in each Z intervals.
     */
    int    *TriPrPartlossZ_m;//counter along Z for histogram
    /**
      * @param TriSePartlossZ_m is a counter of lost particle number in each Z intervals.
      */
    int    *TriSePartlossZ_m;//counter along Z for histogram
    /**
     * @param Triangle_lossmax_m store the maximum number dumped in each triangle.
     */
    int     Triangle_lossmax_m;
    /**
     * @param geo3Dcoords_m store the geometry coordinates in a STL Vector
     * The ID of geo3Dcoords_m is equal to points' ID, and can be referenced by triangle to get vertex coordinates
     */
    vector<Vector_t> geo3Dcoords_m;
    /**
     * @param Triarea_m store area of triangles in a STL Vector
     * The ID of Triarea_m is equal to triangles' ID, and can be referenced by triangle to get its area
     */
    vector<double> Triarea_m;
    /**
     * @param TriBGphysicstag_m store the tags of each boundary triangle for proper physics action.
     */
    vector<short> TriBGphysicstag_m;

    /**
     * Calculate the maximum of coordinates of geometry,i.e the maximum of X,Y,Z
     */
    Vector_t getMaxExtend() {
        const Vector_t x = *max_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessx);
        const Vector_t y = *max_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessy);
        const Vector_t z = *max_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessz);
        return Vector_t(x(0), y(1), z(2));
    }
    /**
     * Calculate the minimum of coordinates of geometry,i.e the minimum of X,Y,Z
     */
    Vector_t getMinExtend() {
        const Vector_t x = *min_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessx);
        const Vector_t y = *min_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessy);
        const Vector_t z = *min_element(geo3Dcoords_m.begin(), geo3Dcoords_m.end(), myLessz);
        return Vector_t(x(0), y(1), z(2));
    }

    /**
     * Used to determine whether a particle hits boundary.
     * @param x stands for the coordinates of the testing particle

     Not sufficient for particle hitting the triangle area.

    */
    bool isInGeometry(Vector_t x) {
        return boundary_ids_m.find(f(x)) != boundary_ids_m.end();
    }
    /**
     * Return the hr_m.
     */
    Vector_t gethr() {
        return hr_m;
    }
    /**
     * Return the nr_m.
     */
    Vektor<int, 3> getnr() {
        return nr_m;
    }
    /**
     * Return the mincoords_m.
     */
    Vector_t getmincoords() {
        return mincoords_m;
    }
    /**
     * Return the maxcoords_m.
     */
    Vector_t getmaxcoords() {
        return maxcoords_m;
    }
private:

    // Not implemented.
    BoundaryGeometry(const BoundaryGeometry &);
    void operator=(const BoundaryGeometry &);

    // Clone constructor.
    BoundaryGeometry(const string &name, BoundaryGeometry *parent);

    /**
     * Recursively get inward normal of triangle with id=idx by switching the order of vertex w.r.t the vertex order of triangle with id = caller;
     * @param caller stands for the id of triangle which has already been aligned.
     * @param idx is the triangle need to be aligned to get the inward normal.
     */
    void orientAllTriangles(size_t caller, size_t idx) {
        // cout << "." << flush;
        orientTriangle(idx, caller);
        isOriented_m.insert(idx);
        vector<size_t> neighbours = findNeighbours(idx);
        for(int i = 0; i < neighbours.size(); i++) {
            if(isOriented_m.find(neighbours[i]) == isOriented_m.end())
                orientAllTriangles(idx, neighbours[i]);
        }
    }


    vector<size_t> findNeighbours(size_t idx) {
        vector<size_t> ret, ret1;
        set<size_t>temp;
        vector<size_t>::iterator retIt;
        map< size_t, vector<size_t> >::iterator it;
        for(int i = 1; i <= 3; i++) {
            it = triangleLookupTable_m.find(triangleCorner(idx, i));
            ret.insert(ret.end(), (*it).second.begin(), (*it).second.end());
        }
        for(retIt = ret.begin(); retIt != ret.end(); retIt++) {
            if(*retIt != idx) {
                if(temp.find(*retIt) == temp.end()) {
                    temp.insert(*retIt);
                } else
                    ret1.push_back(*retIt);
            }
        }
        assert(ret1.size() != 0);
        return ret1;
    }

    size_t triangleCorner(size_t idx, size_t corner) const {
        return allbfaces_m[4*idx + corner];
    }

    void orientTriangle(size_t idx, size_t caller) {
        vector<size_t> id, ic;
        for(int i = 1; i <= 3; i++) {
            for(int j = 1; j <= 3; j++) {
                if(triangleCorner(idx, j) == triangleCorner(caller, i)) {
                    id.push_back(j);
                    ic.push_back(i);
                }
            }
        }
        if((ic[1] - ic[0]) == 1) {
            int idtmp = id[1] - id[0];
            if(idtmp == 1) {
                alignedT_m.push_back(idx);
                allbfaces_m[4*idx+id[0]] = allbfaces_m[4*caller+ic[1]];
                allbfaces_m[4*idx+id[1]] = allbfaces_m[4*caller+ic[0]];
            }
            if(idtmp == -1) {
                NaliT_m.push_back(idx);
            }
            if(idtmp == 2) {
                NaliT_m.push_back(idx);
            }
            if(idtmp == -2) {
                alignedT_m.push_back(idx);
                allbfaces_m[4*idx+id[0]] = allbfaces_m[4*caller+ic[1]];
                allbfaces_m[4*idx+id[1]] = allbfaces_m[4*caller+ic[0]];
            }
        }
        if((ic[1] - ic[0]) == 2) {
            int idtmp = id[1] - id[0];
            if(idtmp == -1) {
                alignedT_m.push_back(idx);
                allbfaces_m[4*idx+id[0]] = allbfaces_m[4*caller+ic[1]];
                allbfaces_m[4*idx+id[1]] = allbfaces_m[4*caller+ic[0]];
            }
            if(idtmp == 1) {
                NaliT_m.push_back(idx);
            }
            if(idtmp == -2) {
                NaliT_m.push_back(idx);
            }
            if(idtmp == 2) {
                alignedT_m.push_back(idx);
                allbfaces_m[4*idx+id[0]] = allbfaces_m[4*caller+ic[1]];
                allbfaces_m[4*idx+id[1]] = allbfaces_m[4*caller+ic[0]];
            }

        }

    }
    /**
      * Recursively get inward triangle normal of all surface triangles.
      *
      * The basic idea is as follow: get the first triangle's inward normal by determine a nearby point
      * is inside or outside boundary geometry( using ray-triangle intersection and even/odd intersection number).
      * Then use a recursion method to switch the vertex order of adjacent triangles. The inward normal is stored
      * in TriNormal_m.
      */
    void makeTriNormal() {

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
        double nomtmp = sqrt(n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
        n = n / nomtmp;
        TriNormal_m.push_back(n);

        Vector_t mytemp;
        mytemp[0] = t0[0] + TriNormal_m[0](0) * 0.00000005;
        mytemp[1] = t0[1] + TriNormal_m[0](1) * 0.00000005;
        mytemp[2] = t0[2] + TriNormal_m[0](2) * 0.00000005;
        if(isInside(mytemp) == 0) {
            if(((mytemp[0] - t0[0])*TriNormal_m[0](0) + (mytemp[1] - t0[1])*TriNormal_m[0](1) + (mytemp[2] - t0[2])*TriNormal_m[0](2)) < 0) {
                // TriNormal_m[i]=TriNormal_m[i];
            } else {
                TriNormal_m[0] = -TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            }
        } else {
            if(((mytemp[0] - t0[0])*TriNormal_m[0](0) + (mytemp[1] - t0[1])*TriNormal_m[0](1) + (mytemp[2] - t0[2])*TriNormal_m[0](2)) < 0) {
                TriNormal_m[0] = -TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            } else {
                TriNormal_m[0] = TriNormal_m[0];
            }
        }
        for(size_t i = 0; i < numbfaces_global_m; i++) { //for every triangle find adjacent triangles for each vertex ;
            for(int j = 1; j <= 3; j++) {
                map< size_t, vector<size_t> >::iterator it;
                it = triangleLookupTable_m.find(allbfaces_m[4*i + j]);
                if(it == triangleLookupTable_m.end()) {
                    vector <size_t> tmp;
                    tmp.push_back(i);
                    triangleLookupTable_m.insert(pair<size_t, vector<size_t> > (allbfaces_m[4*i + j], tmp));
                } else
                    (*it).second.push_back(i);
            }
        }
        /*  int temp = allbfaces_m[5];
        allbfaces_m[5] = allbfaces_m[6];
        allbfaces_m[6] = temp; */
        isOriented_m.insert(0);
        vector<size_t> neighbours = findNeighbours(0);

        for(int i = 0; i < neighbours.size(); i++) {
            if(isOriented_m.find(neighbours[i]) == isOriented_m.end())
                orientAllTriangles(0, neighbours[i]);
        }
        for(int i = 1; i < numbfaces_global_m; i ++) {
            Vector_t t0, t1, t2;
            t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
            t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
            t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
            Vector_t u = t1 - t0;
            Vector_t v = t2 - t0;
            Vector_t n;
            n[0] = u[1] * v[2] - v[1] * u[2];
            n[1] = u[2] * v[0] - v[2] * u[0];
            n[2] = u[0] * v[1] - v[0] * u[1];
            double nomtmp = sqrt(n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
            if(nomtmp != 0) {
                n = n / nomtmp;
            }
            TriNormal_m.push_back(n);
        }
    }






    Vector_t LineInsTri(Vector_t x0, Vector_t x1, size_t i);

    /**
     * Calculate the number of intersects between a line segment and the geometry boundary,make sure x1 is outside the geometry.
     * @param
     */
    vector<Vector_t>  PartBoundaryInteNum(Vector_t x0, Vector_t x1) {

        using namespace std;

        vector<Vector_t> SegDiscrete;
        vector<Vector_t> Isp;
        vector<Vector_t>::iterator TriIscIt;

        vector<int> TriId;
        vector<int>::iterator TriIdIt;


        int Seglen = (((int)floor(sqrt((x0[0] - x1[0]) * (x0[0] - x1[0]) + (x0[1] - x1[1]) * (x0[1] - x1[1]) + (x0[2] - x1[2]) * (x0[2] - x1[2])) / hr_m[0])) + 1);
        int count = 0;
        for(int i = 0; i < Seglen ; i++) {
            SegDiscrete.push_back(x0 + hr_m[0]*i * (x1 - x0) / sqrt((x0[0] - x1[0]) * (x0[0] - x1[0]) + (x0[1] - x1[1]) * (x0[1] - x1[1]) + (x0[2] - x1[2]) * (x0[2] - x1[2])));
        }
        SegDiscrete.push_back(x1);

        vector<Vector_t>::iterator myit;
        for(myit = SegDiscrete.begin() ; myit != SegDiscrete.end() ; myit++) {
            if(!isInGeometry(*myit)) { //segment do not cross the boundary cubic;
                continue;
            } else { //segment cross the boundary cubic;
                int id = f(*myit);
                for(int i = 0; i < numbfaces_global_m; i++) {
                    if(((f(geo3Dcoords_m[allbfaces_m[4*i+1]]) == id) || (f(geo3Dcoords_m[allbfaces_m[4*i+2]]) == id) || (f(geo3Dcoords_m[allbfaces_m[4*i+3]]) == id))) {
                        Vector_t tmp =  LineInsTri(x0, x1, i);
                        TriIdIt = std::find(TriId.begin(), TriId.end(), i);
                        if(((tmp != x1) && ((TriIdIt == TriId.end()) || TriId.size() == 0))) {
                            TriIscIt = std::find(Isp.begin(), Isp.end(), tmp);
                            if((TriIscIt == Isp.end()) || (Isp.size() == 0)) {
                                TriId.push_back(i);
                                Isp.push_back(tmp);
                                ++count;
                            }
                        } else
                            continue;
                    }
                }
            }
        }
        if(count == 0)
            Isp.push_back(x1);
        return Isp;
    }


    /**
     * A struct.
     *
     * Used as a compare function object for STL max_element
     */

    struct myLessclassx {
        bool operator()(Vector_t x1, Vector_t x2) { return x1(0) < x2(0); }
    } myLessx;
    /**
     * A struct.
     *
     * Used as a compare function object for STL max_element
     */
    struct myLessclassy {
        bool operator()(Vector_t x1, Vector_t x2) { return x1(1) < x2(1); }
    } myLessy;
    /**
     * A struct.
     *
     * Used as a compare function object for STL max_element
     */
    struct myLessclassz {
        bool operator()(Vector_t x1, Vector_t x2) { return x1(2) < x2(2); }
    } myLessz;

    /**
     * Calculate the maximum of a triangle length
     * @param i stands for the triangle Id.
     */
    double Trianglelen(int i) {
        Vector_t x1, x2, x3;
        int id1 = allbfaces_m[4*i+1], id2 = allbfaces_m[4*i+2], id3 = allbfaces_m[4*i+3];
        x1 = geo3Dcoords_m[id1];
        x2 = geo3Dcoords_m[id2];
        x3 = geo3Dcoords_m[id3];
        double max =  sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2])) > sqrt((x3[0] - x2[0]) * (x3[0] - x2[0]) + (x3[1] - x2[1]) * (x3[1] - x2[1]) + (x3[2] - x2[2]) * (x3[2] - x2[2])) ? sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2])) : sqrt((x3[0] - x2[0]) * (x3[0] - x2[0]) + (x3[1] - x2[1]) * (x3[1] - x2[1]) + (x3[2] - x2[2]) * (x3[2] - x2[2]));
        max = max > sqrt((x3[0] - x1[0]) * (x3[0] - x1[0]) + (x3[1] - x1[1]) * (x3[1] - x1[1]) + (x3[2] - x1[2]) * (x3[2] - x1[2])) ? max : sqrt((x3[0] - x1[0]) * (x3[0] - x1[0]) + (x3[1] - x1[1]) * (x3[1] - x1[1]) + (x3[2] - x1[2]) * (x3[2] - x1[2]));
        return (max);
    }




    /**
     * Calculate the maximum dimention of triangles. This value will be used to define the cubic box size
     */
    double getMaxDimenssion() {
        double maximum = 0, min = 0.01;
        for(int i = 0; i < numbfaces_global_m; i++) {
            maximum = maximum >= Trianglelen(i) ? maximum : Trianglelen(i);
            min = min < Trianglelen(i) ? min : Trianglelen(i);
        }
        INFOMSG("size of Triangle " << maximum << " " << min << endl);
        return maximum ;
    }


    /**
     * Map a 3D coordinate to a cubic box with a unique Id
     *
     * @param x is an OPAL Vector_t type, stands for the 3D coordinates need to be mapped.
     */
    int f(Vector_t x) {
        /**
         We know: lenght of the structure len_m
         Mesh size: hr_m
         Number of mesh points nr_m
        */
        size_t ret = 0;
        size_t id_tx, id_ty, id_tz;
        id_tx = floor((x[0] - mincoords_m[0]) / hr_m[0]);
        id_ty = floor((x[1] - mincoords_m[1]) / hr_m[1]);
        id_tz = floor((x[2] - mincoords_m[2]) / hr_m[2]);
        ret = 1 + id_tz * nr_m[0] * nr_m[1] + id_ty * nr_m[0] + id_tx;
        return ret;
    }

    /**
     * Return the larger one.
     */

    bool gT(double x, double y) {
        return (x >= y); //Used for initializing the particles
    }
    /**
     * Leading zeros.
     * @param number stands for specified bits number
     */

    string convert2Int(int number) {

        stringstream ss;//create a stringstream
        ss <<   setw(5) << setfill('0') <<  number; //add number to the stream
        return ss.str();//return a string with the contents of the stream
    }
    /**
     * A struct.
     *
     * member rI stands for intersection position w.r.t line segment.
     * Isc stores the coordinates of intersection point.
     */
    struct Intersection {
        double rI;
        Vector_t Isc;
    };

    /**
     * Find a intersection between a line segment and triangle,faster by using the pre-computed oriented normal.
     *
     * @param x0 and @param x1 stands for start point and end point of line segment.
     * @param i stands for id of triangle.
     * Basic algorithms are as follows:
     * 1)find the intersection between line segment \f$\vec{x1-x0}\f$ and plane defined by point t0 and triangle normal;
     *   if the dot product of line segment and plane normal equals to zero and x0 is not the barycentric point of the triangle(in the plane),
     *   then the line segment is parallel to plane return no intersection, else if the particle is really move (\f$ x0 \neq x1 \f$),return
     *   initialized position-triangle barycentric point as intersection.
     *   The intersection position rI w.r.t x0 is obtained from: \f$ rI=\frac{\vec{n} \cdot \vec{(t0-x0)}}{\vec{n} \cdot \vec{(x1-x0)}} \f$,
     *   where t0 is the first vertex of triangle, n is the normal of triangle.
     *   The intersection point Itsec is obtained from: \f$ Itsec=x0+rI(x1-x0) \f$.
     * 2)check if the intersection point is inside the triangle by using parametric coordinates sI and tI of the intersecion point.
     *   First calculate sI and tI. The parametric plane equation is given by:  \f$ t(sI,tI)=t0+sI(t1-t0)+tI(t2-t0)=t0+sI\vec{u}+tI\vec{v} \f$.
     *   \f$\vec{w}=\vec{Itsec-t0}\f$ is also in the plane, solve equation: \f$\vec{w}=t0+sI\vec{u}+tI\vec{v}\f$ , we obtain the sI and tI.
     *   \f$ sI=\frac{(\vec{u} \cdot \vec{v})(\vec{w} \cdot \vec{v})-(\vec{v} \cdot \vec{v})(\vec{w} \cdot \vec{u})}{(\vec{u} \cdot \vec{v})^2-(\vec{u} \cdot \vec{u})(\vec{v} \cdot \vec{v})} \f$,
     *   \f$ tI=\frac{(\vec{u} \cdot \vec{v})(\vec{w} \cdot \vec{u})-(\vec{u} \cdot \vec{u})(\vec{w} \cdot \vec{v})}{(\vec{u} \cdot \vec{v})^2-(\vec{u} \cdot \vec{u})(\vec{v} \cdot \vec{v})} \f$.
     *   If \f$ sI \geq 0 \f$, \f$ tI \geq 0 \f$ and \f$ sI+tI \leq 1 \f$, then the intersection is inside the triangle, and return the intersection
     *   coordinate Itsec.
     */
    Intersection FindIntersection(Vector_t x0, Vector_t x1, size_t i) {
        Vector_t lseg = x1 - x0; //stands for the length and direction of line segment;
        Vector_t t0, t1, t2;
        t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
        t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
        t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
        Vector_t u = t1 - t0; //side 1 of triangle;
        Vector_t v = t2 - t0; //side 2 of triangle;
        Vector_t lt = t0 - x0;
        Vector_t n = TriNormal_m[i];

        Intersection myInter;
        myInter.rI = 0;
        myInter.Isc = maxcoords_m;
        double dotLT = Dotproduct(n, lseg);
        if(dotLT == 0) {
            if((x0 == Tribarycent_m[i]) && (x0 != x1)) { //Some initialized particles have momenta parallel to its triangle normal,this kind of particles will lose directly
                myInter.Isc = Tribarycent_m[i];
            }
            return myInter;
        } else {

            myInter.rI = Dotproduct(n, lt) / dotLT; //find intersection position w.r.t x0 and the unit is (x1-x0);
            Vector_t ItSec = x0 + myInter.rI * lseg; //find the coordinate of intersection.
            //find if the intersection is inside the triangle.
            Vector_t w = ItSec - t0;
            double tmp1, tmp2, tmp3, tmp4, tmp5;
            tmp1 = Dotproduct(u, v);
            tmp2 = Dotproduct(w, v);
            tmp3 = Dotproduct(u, w);
            tmp4 = Dotproduct(u, u);
            tmp5 = Dotproduct(v, v);
            double sI = (tmp1 * tmp2 - tmp5 * tmp3) / (tmp1 * tmp1 - tmp4 * tmp5);
            double tI = (tmp1 * tmp3 - tmp4 * tmp2) / (tmp1 * tmp1 - tmp4 * tmp5);
            if((sI >= 0) && (tI >= 0) && ((sI + tI) <= 1)) {
                myInter.Isc = ItSec;
                return myInter;
            } else {
                //Intesect is on the extended plane, return
                return myInter;
            }

        }
    }



    /**
     * Get the smallest bounding box (with small margin) of id th triangle.The smallest bounding box is used to make sure that all part of a triangle is known by boundary bounding box.
     * @param id  stands for id of triangle.
     */

    vector<Vector_t> SetMinMaxBound(size_t id) {
        Vector_t min = geo3Dcoords_m[allbfaces_m[4*id+1]];
        Vector_t max = min;
        for(int i = 2; i <= 3; i++) {
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](0) < min[0])
                min[0] = geo3Dcoords_m[allbfaces_m[4*id+i]](0);
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](1) < min[1])
                min[1] = geo3Dcoords_m[allbfaces_m[4*id+i]](1);
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](2) < min[2])
                min[2] = geo3Dcoords_m[allbfaces_m[4*id+i]](2);
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](0) > max[0])
                max[0] = geo3Dcoords_m[allbfaces_m[4*id+i]](0);
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](1) > max[1])
                max[1] = geo3Dcoords_m[allbfaces_m[4*id+i]](1);
            if(geo3Dcoords_m[allbfaces_m[4*id+i]](2) > max[2])
                max[2] = geo3Dcoords_m[allbfaces_m[4*id+i]](2);
        }
        Vector_t temp = 0.25 * hr_m; //add some margins to make sure that boundary bounding box has both positive and negtive margins w.r.t boundary.Just make sure that the size of bounding box for triangle is smaller than the boundary bounding box.
        vector<Vector_t> ret;
        min = min - temp;
        max = max + temp;
        ret.push_back(min);
        ret.push_back(max);
        assert(ret.size() != 0);
        return ret;
    }

    vector<Vector_t> SetRayBound(Vector_t x0, Vector_t x1) {
        Vector_t min = x0;
        Vector_t max = min;

        if(x1[0] < min[0])
            min[0] = x1[0];
        else
            max[0] = x1[0];
        if(x1[1] < min[1])
            min[1] = x1[1];
        else
            max[1] = x1[1];
        if(x1[2] < min[2])
            min[2] = x1[2];
        else
            max[2] = x1[2];

        vector<Vector_t> ret;
        ret.push_back(min);
        ret.push_back(max);
        assert(ret.size() != 0);
        return ret;
    }

    double TriangleArea(int id) {
        Vector_t AB = geo3Dcoords_m[allbfaces_m[4*id+2]] - geo3Dcoords_m[allbfaces_m[4*id+1]];
        Vector_t AC = geo3Dcoords_m[allbfaces_m[4*id+3]] - geo3Dcoords_m[allbfaces_m[4*id+1]];

        return(0.5 * sqrt(Dotproduct(AB, AB) * Dotproduct(AC, AC) - Dotproduct(AB, AC) * Dotproduct(AB, AC)));
    }



    inline double Dotproduct(Vector_t x0, Vector_t x1) {
        return (x0[0] * x1[0] + x0[1] * x1[1] + x0[2] * x1[2]);
    }

    // The particle reference data.
    PartData reference;


    /**
     * @param bfaces_idx_m store the Id number of a triangle,not in use at present.
     */
    int    *bfaces_idx_m;
    /**
     * @param  nof_sym_m stores the count number of symmetry planes.
     */
    int     nof_sym_m;

    /**
     * @param boundary_ids_m stores the Ids of triangles which form the boundary in a STL set
     */
    set<size_t> boundary_ids_m;
    /**
     * @param len_m is a OPAL Vector_t type parameter stores the length of geometry in 3D Cartesian coordinates.
     */
    Vector_t      len_m;
    /**
     * @param hr_m is a OPAL Vector_t type parameter stores the length of cubic box.
     */
    Vector_t      hr_m;
    /**
     * @param nr_m is a OPAL Vector_t type parameter stores the number of intervals of geometry in 3D Cartesian coordinates.
     */
    Vektor<int, 3> nr_m;
    /**
     * @param mincoords_m stores minimum of geometry coordinate.
     */
    Vector_t      mincoords_m;
    /**
     * @param maxcoords_m stores maximum of geometry coordinate.
     */
    Vector_t      maxcoords_m;
    /**
     * @param partsr_m stores coordinates of particles using STL vector.
     */
    vector<Vector_t> partsr_m;  // add particle positions

    /**
     * @param lostedPart_m stores coordinates of losted particles
     */
    vector<Vector_t> lostedPart_m;  // store particles hitting the surface
    /**
     * @param CubicLookupTable_m is a stl map type variable.
     * The first member of the map is the id of a boundary cubic.
     * The second member is a vector containing the ids of triangles inside the boundary cubic .
     */
    map< size_t, vector<size_t> >  CubicLookupTable_m;
    /**
     * @param isOriented_m stores the ids of oriented triangles.
     */
    set<size_t> isOriented_m;
    /**
     * @param triangleLookupTable_m stores the ids of vertex as map key and the ids of triangles which contains the vertex as the map value.
     */
    map< size_t, vector<size_t> > triangleLookupTable_m;
    /**
     * @param NaliT_m stores the ids of Triangles need to be oriented.
     */
    vector<size_t> NaliT_m;
    /**
     * @param alignedT_m stores the ids of oriented Triangles.
     */
    vector<size_t> alignedT_m;
    /**
     * @param alignedT_m stores the ids of oriented Triangles.
     */
    map<int, vector<Vector_t> > Boxcoords_m;

};

inline Inform &operator<<(Inform &os, const BoundaryGeometry &b) {
    return b.print(os);
}


#endif // OPAL_BOUNDARY_GEOMETRY_HH
