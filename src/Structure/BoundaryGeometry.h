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
//#include "Structure/SecondaryEmissionPhysics.h"
#include "Distribution/ranlib.h"
#include <vector>
#include <sstream>
#include <set>
#include <cmath>
//#include <algorithm>
//#include "Ippl.h"
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
class SecondaryEmissionPhysics;

class BoundaryGeometry: public Definition {

public:

    /// Exemplar constructor.
    BoundaryGeometry();

    virtual ~BoundaryGeometry();

    /// Test if replacement is allowed.
    //  Can replace only by another GEOMETRY.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual BoundaryGeometry *clone(const std::string &name);

    /// Check the GEOMETRY data.
    virtual void execute();

    /// Find named GEOMETRY.
    static BoundaryGeometry *find(const std::string &name);

    /// Update the GEOMETRY data.
    virtual void update();

    void updateElement(ElementBase *element);

    void initialize();
    void createParticlesOnSurface(size_t n, double darkinward, OpalBeamline &itsOpalBeamline, PartBunch &itsBunch);
    void createPriPart(size_t n, double darkinward,  OpalBeamline &itsOpalBeamline, PartBunch *itsBunch);

    // inline size_t getN() { return geo3Dcoords_m.size(); }
    size_t getN() ;
    inline Vector_t getCooridinate(size_t i) { return partsr_m[i]; }
    inline void clearCooridinateArray() { return partsr_m.clear(); }
    inline Vector_t getMomenta(size_t i) { return partsp_m[i]; }//benchmark code
    inline void clearMomentaArray() { return partsp_m.clear(); }//benchmark code

    //    vector<Vector_t> PartBoundaryInteNum(Vector_t x0,Vector_t x1);
    bool isInside(Vector_t x);
    std::vector<Vector_t>  GridIntersection(Vector_t x0, Vector_t x1);
    Inform &printInfo(Inform &os) const;
    std::string getTopology() const;
    std::string getFilename() const;

    double getA();
    double getB();
    double getC();
    double getS();
    double getLenght();
    double getL1();
    double getL2();

    std::string getDistribution();
    std::vector<std::string> getDistributionArray();
    
    double getZshift();

    double getXYZScale();

    void makeBoundaryIndexSet();

    int PartInside(const Vector_t r, const Vector_t v, const double dt, int Parttype, const double Qloss, Vector_t &intecoords, int &triId, double &Energy);
    
    void setBGphysicstag();

    int doBGphysics(const Vector_t &intecoords, const int &triId); // non secondary emission version.

    int doBGphysics(const Vector_t &intecoords, const int &triId, const double &incEnergy, const double &incQ, const Vector_t &incMomentum, PartBunch *itsBunch, double &seyNum);// call Furman-Pivi's model

    int doBGphysics(const Vector_t &intecoords, const int &triId, const double &incEnergy, const double &incQ, const Vector_t &incMomentum, PartBunch *itsBunch, double &seyNum, const int &para_null);// call Vaughan's model


    void setNEmissionMode(bool nEmissionMode);

    void setEInitThreshold(double einitthreshold);

    void setWorkFunction(double workFunction);
    void setFieldEnhancement(double fieldEnhancement);
    void setMaxFN(size_t maxFNemission);
    void setFNTreshold(double fieldFNthreshold);
    void setFNParameterA(double parameterFNA);
    void setFNParameterB(double parameterFNB);
    void setFNParameterY(double parameterFNY);
    void setFNParameterVYZe(double parameterFNVYZe);
    void setFNParameterVYSe(double parameterFNVYSe);
    size_t doFNemission(OpalBeamline &itsOpalBeamline, PartBunch *itsBunch, const double t);
    void setBoundaryMatType(int BoundaryMatType) ;
    void setvSeyZero(double vSeyZero);// return sey_0 in Vaughan's model
    void setvEZero(double vEZero);// return the energy related to sey_0 in Vaughan's model
    void setvSeyMax(double vSeyMax);// return sey max in Vaughan's model
    void setvEmax(double vEmax);// return Emax in Vaughan's model
    void setvKenergy(double vKenergy);// return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
    void setvKtheta(double vKtheta);// return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
    void setvVThermal(double vVThermal);// return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
    void setVw(double ppVw);

    bool nEmissionMode_m;
    double eInitThreshold_m;

    double vSeyZero_m;// energy related to sey_0 in Vaughan's model

    double vEzero_m;// sey_0 in Vaughan's model

    double vSeyMax_m;// sey max in Vaughan's model

    double vEmax_m;// Emax in Vaughan's model

    double vKenergy_m;// fitting parameter denotes the roughness of surface for impact energy in Vaughan's model

    double vKtheta_m;// fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
    double vVThermal_m;// return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model

    /// @param ppVw_m denotes the velocity scalar for Parallel plate benchmark.
    
    double ppVw_m;
    
    double triangle_max_m;
    double triangle_min_m;

    /**
     * @param seBoundaryMatType_m stores the user defined material type number for secondary emission model.
     */
    int seBoundaryMatType_m;

    /**
     * @param workFunction_m stores the user defined work function for Fowler-Nordheim model.
     */
    double workFunction_m;

    /**
     * @param fieldEnhancement_m stores the user defined field enhancement factor for Fowler-Nordheim model.
     */
    double fieldEnhancement_m;

    /**
     * @param maxFNemission_m stores the user defined maximum emitted particle number per triangle for Fowler-Nordheim model.

     */
    size_t maxFNemission_m;

    /**
     * @param fieldFNthreshold_m stores the user defined lower threshold of electric field for Fowler-Nordheim model.
     */
    double fieldFNthreshold_m;

    /**
     * @param parameterFNA_m stores the user defined parameter A for Fowler-Nordheim model.Default value:\f$1.54\times 10^{-6}\f$.       */
    double parameterFNA_m;

    /**
     * @param parameterFNB_m stores the user defined parameter B for Fowler-Nordheim model.Default value:\f$6.83\times 10^9\f$.
     */
    double parameterFNB_m;

    /**
     * @param parameterFNY_m stores the user defined parameter for Fowler-Nordheim parameter y.Default value:\f$3.795\times 10^{-5}\f$.
     */
    double parameterFNY_m;

    /**
     * @param parameterFNVYZe_m stores the user defined parameter for zero order fit constant for v(y) of Fowler-Nordheim model.Default value:\f$0.9632\f$.
     */
    double parameterFNVYZe_m;

    /**
     * @param parameterFNVYSe_m stores the user defined parameter for second order fit constant for v(y) of Fowler-Nordheim model.Default value:\f$1.065\f$.
     */
    double parameterFNVYSe_m;

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
    std::vector<Vector_t>    TriNormal_m;

    /**
     * @param TriPrPartloss_m store the number of primary particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    double    *TriPrPartloss_m;

    /**
     * @param TriSePartloss_m store the number of secondary particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    double    *TriSePartloss_m;

    /**
     * @param TriFEPartloss_m store the number of field emission/darkcurrent particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    double    *TriFEPartloss_m;

    /**
     * @param TriPrPartlossZ_m is a counter of lost primary particle number in each Z intervals.
     */
    double    *TriPrPartlossZ_m;//counter along Z for histogram

    /**
     * @param TriSePartlossZ_m is a counter of lost secondary particle number in each Z intervals.
     */
    double    *TriSePartlossZ_m;//counter along Z for histogram

    /**
     * @param TriFEPartlossZ_m is a counter of lost field emission and dark current particle number in each Z intervals.
     */
    double    *TriFEPartlossZ_m;//counter along Z for histogram

    /**
     * @param geo3Dcoords_m store the geometry coordinates in a STL Vector
     * The ID of geo3Dcoords_m is equal to points' ID, and can be referenced by triangle to get vertex coordinates
     */
    std::vector<Vector_t> geo3Dcoords_m;

    /**
     * @param Triarea_m store area of triangles in a STL Vector
     * The ID of Triarea_m is equal to triangles' ID, and can be referenced by triangle to get its area
     */
    std::vector<double> Triarea_m;

    /**
     * @param TriBGphysicstag_m store the tags of each boundary triangle for proper physics action.
     */
    std::vector<short> TriBGphysicstag_m;

    std::string h5FileName_m;

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

     * Not sufficient for particle hitting the triangle area.

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

    /**
     * @param xyzscale_m scale of x,y and z coordinates
     */
    double xyzscale_m;

    /**
     * @param TPInside_m, @param TPreProc_m, @param TRayTrace_m are timers for out put the time consumptions.
     */
    IpplTimings::TimerRef TPInside_m;
    IpplTimings::TimerRef TPreProc_m;
    IpplTimings::TimerRef TRayTrace_m;
    IpplTimings::TimerRef Tinward_m;

    SecondaryEmissionPhysics *sec_phys_m;

    // Not implemented.
    BoundaryGeometry(const BoundaryGeometry &);
    void operator=(const BoundaryGeometry &);

    // Clone constructor.
    BoundaryGeometry(const std::string &name, BoundaryGeometry *parent);

    /**
     * Recursively get inward normal of triangle with id=idx by switching the order of vertex w.r.t the vertex order of triangle with id = caller;
     * @param caller stands for the id of triangle which has already been aligned.
     * @param idx is the triangle need to be aligned to get the inward normal.
     */
    void orientAllTriangles(size_t caller, size_t idx) {
        // cout << "." << flush;
        orientTriangle(idx, caller);
        isOriented_m.insert(idx);
        std::vector<size_t> neighbours = findNeighbours(idx);
        for(unsigned int i = 0; i < neighbours.size(); i++) {
            if(isOriented_m.find(neighbours[i]) == isOriented_m.end())
                orientAllTriangles(idx, neighbours[i]);
        }
    }


    std::vector<size_t> findNeighbours(size_t idx) {
        std::vector<size_t> ret, ret1;
        std::set<size_t>temp;
        std::vector<size_t>::iterator retIt;
        std::map< size_t, std::vector<size_t> >::iterator it;
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
        PAssert(ret1.size() != 0);
        return ret1;
    }

    size_t triangleCorner(size_t idx, size_t corner) const {
        return allbfaces_m[4*idx + corner];
    }

    void orientTriangle(size_t idx, size_t caller) {
        std::vector<size_t> id, ic;
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
        mytemp[0] = t0[0] + TriNormal_m[0](0) * 0.1*triangle_min_m;// fixme: if boundary size smaller than 0.0005m? maybe use triangle_min_m
        mytemp[1] = t0[1] + TriNormal_m[0](1) * 0.1* triangle_min_m;
        mytemp[2] = t0[2] + TriNormal_m[0](2) * 0.1* triangle_min_m;
	//*gmsg<<"mytemp: "<<mytemp<<" TriNormal_m[0] "<<TriNormal_m[0]<<endl;
        if(!isInside(mytemp)) {
            if(dot(mytemp - t0, TriNormal_m[0]) < 0) {
                // TriNormal_m[i]=TriNormal_m[i];
            } else {
                TriNormal_m[0] = -TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            }
        } else {
            if(dot(mytemp - t0, TriNormal_m[0]) < 0) {
                TriNormal_m[0] = -TriNormal_m[0];
                int temp = allbfaces_m[2];
                allbfaces_m[2] = allbfaces_m[3];
                allbfaces_m[3] = temp;
            } else {
                TriNormal_m[0] = TriNormal_m[0];
            }
        }
	//*gmsg<<"mytemp: "<<mytemp<<" TriNormal_m[0] after align: "<<TriNormal_m[0]<<" isInside(mytemp) "<<isInside(mytemp)<<endl;
       
        for(int i = 0; i < numbfaces_global_m; i++) { //for every triangle find adjacent triangles for each vertex ;
            for(int j = 1; j <= 3; j++) {
                std::map< size_t, std::vector<size_t> >::iterator it;
                it = triangleLookupTable_m.find(allbfaces_m[4*i + j]);
                if(it == triangleLookupTable_m.end()) {
                    std::vector <size_t> tmp;
                    tmp.push_back(i);
                    triangleLookupTable_m.insert(std::pair<size_t, std::vector<size_t> > (allbfaces_m[4*i + j], tmp));
                } else
                    (*it).second.push_back(i);
            }
        }
	isOriented_m.insert(0);
        std::vector<size_t> neighbours = findNeighbours(0);

        for(unsigned int i = 0; i < neighbours.size(); i++) {
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
     * Calculate the number of intersects between a line segment and the geometry boundary, make sure x1 is outside the geometry.
     * @param
     */
    /*std::vector<Vector_t>  PartBoundaryInteNum(Vector_t x0, Vector_t x1) {

        using namespace std;
	IpplTimings::startTimer(Tinward_m);
        std::vector<Vector_t> SegDiscrete;
        std::vector<Vector_t> Isp;
        std::vector<Vector_t>::iterator TriIscIt;

        std::vector<int> TriId;
        std::vector<int>::iterator TriIdIt;
        double hr_tmp=hr_m[0];
	int dim = 3;
	for (int i=1; i<dim; i++) { // interval should be no larger than the minimum size of box;
	    
	    if(hr_m[i]<hr_tmp)
		hr_tmp = hr_m[i];

	}
	hr_tmp *= 1;
	int Seglen = (int)floor(sqrt(dot(x0 - x1, x0 - x1)) / hr_tmp) + 1;
        int count = 0;
        for(int i = 0; i < Seglen ; i++) {
            SegDiscrete.push_back(x0 + hr_tmp*i * (x1 - x0) / sqrt(dot(x0 - x1, x0 - x1)));
        }
        SegDiscrete.push_back(x1);

        std::vector<Vector_t>::iterator myit;
        for(myit = SegDiscrete.begin() ; myit != SegDiscrete.end() ; myit++) {
            if(!isInGeometry(*myit)) { // segment do not cross the boundary cubic;
                continue;
            } else { // segment cross the boundary cubic;
                int id = f(*myit);
                for(int i = 0; i < numbfaces_global_m; i++) {
		  if(((f(geo3Dcoords_m[allbfaces_m[4*i+1]]) == id) || (f(geo3Dcoords_m[allbfaces_m[4*i+2]]) == id) || (f(geo3Dcoords_m[allbfaces_m[4*i+3]]) == id))) {// Fixme: if the triangle is larger than the box, this will cause problem...
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
	IpplTimings::stopTimer(Tinward_m);
        if(count == 0)
            Isp.push_back(x1);
        return Isp;
	}*/
    std::vector<Vector_t>  PartBoundaryInteNum(Vector_t x0, Vector_t x1) {
        std::vector<Vector_t> Isp;
        for(int i = 0; i < numbfaces_global_m; i++) {
	    
	    Vector_t tmp =  LineInsTri(x0, x1, i);

	    if(tmp != x1) {
	        Isp.push_back(tmp);
	    }

        }
	if (Isp.size()==0)
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
    void getMaxDimenssion() {
        double maximum = 0, min = 0.01;
        for(int i = 0; i < numbfaces_global_m; i++) {
            maximum = maximum >= Trianglelen(i) ? maximum : Trianglelen(i);
            min = min < Trianglelen(i) ? min : Trianglelen(i);
        }
       	triangle_max_m = maximum;
	triangle_min_m = min;
        
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

    std::string convert2Int(int number) {

        std::stringstream ss;//create a stringstream
        ss << std::setw(5) << std::setfill('0') <<  number; //add number to the stream
        return ss.str();// return a string with the contents of the stream
    }


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


    void FindIntersection(const Vector_t &x0, const Vector_t &x1, const size_t &i, double &rI, Vector_t &Isc) {

        const Vector_t lseg = x1 - x0; // stands for the length and direction of line segment;

        const Vector_t t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
        const Vector_t t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
        const Vector_t t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
        const Vector_t u = t1 - t0; // side 1 of triangle;
        const Vector_t v = t2 - t0; // side 2 of triangle;
        const Vector_t lt = t0 - x0;
        const Vector_t n = TriNormal_m[i];


        const double dotLT = dot(n, lseg);
        if(dotLT == 0) {

            if((x0 == Tribarycent_m[i]) && (x0 != x1)) { // Some initialized particles have momenta parallel to its triangle normal,this kind of particles will lose directly
                Isc = Tribarycent_m[i];
            }

        } else {

            rI = dot(n, lt) / dotLT; // find intersection position w.r.t x0 and the unit is (x1-x0);
            const Vector_t ItSec = x0 + rI * lseg; // find the coordinate of intersection.
            // find if the intersection is inside the triangle.
            const Vector_t w = ItSec - t0;
            const double tmp1 = dot(u, v);
            const double tmp2 = dot(w, v);
            const double tmp3 = dot(u, w);
            const double tmp4 = dot(u, u);
            const double tmp5 = dot(v, v);
            const double temp = (tmp1 * tmp1 - tmp4 * tmp5);
            const double sI = (tmp1 * tmp2 - tmp5 * tmp3) / temp;
            const double tI = (tmp1 * tmp3 - tmp4 * tmp2) / temp;
            if((sI >= 0.0) && (tI >= 0.0) && ((sI + tI) <= 1.0)) {
                Isc = ItSec;

            }

        }
    }



    /**
     * Get the smallest bounding box (with small margin) of id th triangle.The smallest bounding box is used to make sure that all part of a triangle is known by boundary bounding box.
     * @param id  stands for id of triangle.
     */

    std::vector<Vector_t> SetMinMaxBound(size_t id) {
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
        Vector_t temp = 1.0*hr_m; // add some margins to make sure that boundary bounding box has both positive and negtive margins w.r.t boundary. Just make sure that the size of bounding box for triangle is smaller than the boundary bounding box.//0.25*hr_m
        std::vector<Vector_t> ret;
        Vector_t min_ext = min - temp;
	Vector_t max_ext = max + temp;
        ret.push_back(min);
        ret.push_back(max);
	ret.push_back(min_ext);
	ret.push_back(max_ext);
        PAssert(ret.size() != 0);
        return ret;
    }
    /**
     * Calculate the area of triangle.
     * @param id  stands for id of triangle.
     */

    double TriangleArea(int id) {
        Vector_t AB = geo3Dcoords_m[allbfaces_m[4*id+2]] - geo3Dcoords_m[allbfaces_m[4*id+1]];
        Vector_t AC = geo3Dcoords_m[allbfaces_m[4*id+3]] - geo3Dcoords_m[allbfaces_m[4*id+1]];

        return(0.5 * sqrt(dot(AB, AB) * dot(AC, AC) - dot(AB, AC) * dot(AB, AC)));
    }




   

    /**
     * @param  nof_sym_m stores the count number of symmetry planes.
     */
    int     nof_sym_m;

    /**
     * @param boundary_ids_m stores the Ids of triangles which form the boundary in a STL set
     */
    std::set<size_t> boundary_ids_m;
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
     * @param partsp_m stores momenta of particles using STL vector.//for benchmark code.
     */
    std::vector<Vector_t> partsp_m;  // add particle momenta

    /**
     * @param partsr_m stores coordinates of particles using STL vector.
     */
    std::vector<Vector_t> partsr_m;  // add particle positions

   
    /**
     * @param CubicLookupTable_m is a stl map type variable.
     * The first member of the map is the id of a boundary cubic.
     * The second member is a vector containing the ids of triangles inside the boundary cubic .
     */
    std::map< size_t, std::vector<size_t> >  CubicLookupTable_m;
    /**
     * @param isOriented_m stores the ids of oriented triangles.
     */
    std::set<size_t> isOriented_m;
    /**
     * @param triangleLookupTable_m stores the ids of vertex as map key and the ids of triangles which contains the vertex as the map value.
     */
    std::map< size_t, std::vector<size_t> > triangleLookupTable_m;
    /**
     * @param NaliT_m stores the ids of Triangles need to be oriented.
     */
    std::vector<size_t> NaliT_m;
    /**
     * @param alignedT_m stores the ids of oriented Triangles.
     */
    std::vector<size_t> alignedT_m;
   
    /**
     * Define an outside domain point
     */
    Vector_t out_m;



};

inline Inform &operator<<(Inform &os, const BoundaryGeometry &b) {
    return b.printInfo(os);
}

namespace BGphysics {
    enum TPHYACTION {
        Nop = 0x01,// tringle is transparent to particle like beam window
        Absorption = 0x02,// triangle has no field emission and secondary emission
        FNEmission = 0x04,// triangle has field emission
        SecondaryEmission = 0x08// trangle has secondary emission 
    };
}
#endif // OPAL_BOUNDARY_GEOMETRY_HH
