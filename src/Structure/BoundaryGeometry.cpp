// ------------------------------------------------------------------------
// $RCSfile: Geometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//   The class for the OPAL GEOMETRY command.
//
// $Date: 2003/08/11 22:09:00 $
// $Author: A. Adelmann $
//
// ------------------------------------------------------------------------

//FIXME: cleanup
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Elements/OpalBeamline.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "hdf5.h"
#include "Distribution/Distribution.h"
#include <sys/stat.h>
#include "Ippl.h"
using namespace Expressions;
using namespace Physics;
using Physics::c;

// Class BoundaryGeometry
// ------------------------------------------------------------------------

// The attributes of class Geometry.
namespace {
    enum {
        FGEOM,     // file holding the geometry
        LENGTH,   // length of elliptic tube
        S,        // start of the geometry
        A,          // major semi-axis of elliptic tube
        B,        // minor semi-axis of ellitpic tube
        TOPO,     // BOX, ELLIPTIC if FGEOM is selected topo is over-written
        DISTR,    // Add distribution to generate on the surface
        SIZE
    };

}

namespace BGphysics {
    enum TPHYACTION {
        Nop = 0x01,
        Absorption = 0x02,
        FNEmission = 0x04,
        SecondaryEmission = 0x08
    };
}


BoundaryGeometry::BoundaryGeometry() : Definition(SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry.") {
    itsAttr[FGEOM] = Attributes::makeString
                     ("FGEOM", "Specifies the geometry file [h5fed]", "");

    itsAttr[TOPO] = Attributes::makeString
                    ("TOPO", "BOX, ELLIPTIC if FGEOM is selected topo is over-written ", "ELLIPTIC");

    itsAttr[LENGTH] = Attributes::makeReal
                      ("LENGTH", "Specifies the length of a tube shaped elliptic beam pipe [m]", 1.0);

    itsAttr[S] = Attributes::makeReal
                 ("S", "Specifies the start of a tube shaped elliptic beam pipe [m]", 0.0);

    itsAttr[A] = Attributes::makeReal
                 ("A", "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);

    itsAttr[B] = Attributes::makeReal
                 ("B", "Specifies the minor semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);


    itsAttr[DISTR] = Attributes::makeString
                     ("DISTR", "Distribution to be generated on the surface", "");


    BoundaryGeometry *defGeometry = clone("UNNAMED_GEOMETRY");
    defGeometry->builtin = true;

    try {
        defGeometry->update();
        OPAL.define(defGeometry);
    } catch(...) {
        delete defGeometry;
    }

}


BoundaryGeometry::BoundaryGeometry(const string &name, BoundaryGeometry *parent):
    Definition(name, parent) {
    //empty so far
}


BoundaryGeometry::~BoundaryGeometry() {
    if(allbfaces_m)
        delete allbfaces_m;
    if(bfaces_idx_m)
        delete bfaces_idx_m;
    if(TriPrPartloss_m)
        delete TriPrPartloss_m;
    if(TriSePartloss_m)
        delete TriSePartloss_m;
    if(Tribarycent_m)
        delete Tribarycent_m;
    if(TriPrPartlossZ_m)
        delete TriPrPartlossZ_m;
    if(TriSePartlossZ_m)
        delete TriSePartlossZ_m;
}


bool BoundaryGeometry::canReplaceBy(Object *object) {
    // Can replace only by another GEOMETRY.
    return dynamic_cast<Geometry *>(object) != 0;
}


BoundaryGeometry *BoundaryGeometry::clone(const string &name) {
    return new BoundaryGeometry(name, this);
}


void BoundaryGeometry::execute() {
    update();
}


BoundaryGeometry *BoundaryGeometry::find(const string &name) {
    BoundaryGeometry *geom = dynamic_cast<BoundaryGeometry *>(OPAL.find(name));

    if(geom == 0)
        throw OpalException("BoundaryGeometry::find()", "Geometry \"" + name + "\" not found.");

    return geom;
}


string BoundaryGeometry::getFilename() {
    return (string)Attributes::getString(itsAttr[FGEOM]);
}

string BoundaryGeometry::getTopology() {
    return (string)Attributes::getString(itsAttr[TOPO]);
}

double BoundaryGeometry::getA() {
    return (double)Attributes::getReal(itsAttr[A]);
}

double BoundaryGeometry::getB() {
    return (double)Attributes::getReal(itsAttr[B]);
}

double BoundaryGeometry::getS() {
    return (double)Attributes::getReal(itsAttr[S]);
}

string BoundaryGeometry::getDistribution() {

    return (string)Attributes::getString(itsAttr[DISTR]);
}


void BoundaryGeometry::update() {
    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_GEOMETRY");
}

void BoundaryGeometry::setBGphysicstag() {
    for(int i = 0; i < numbfaces_global_m; i++) {
        short temp = 0x00;
        if((Tribarycent_m[i](2) > maxcoords_m[2] - 0.00001) && (Tribarycent_m[i](2) < maxcoords_m[2] + 0.00001)) {
            temp = BGphysics::Nop;
        } else {
            temp = (BGphysics::Absorption) | (BGphysics::FNEmission);
        }
        TriBGphysicstag_m.push_back(temp);
    }
}
int BoundaryGeometry::doBGphysics(Vector_t intecoords, int triId) {
    short BGtag = TriBGphysicstag_m[triId];
    int ret = 0;
    if(BGtag & (BGphysics::Nop) == (BGphysics::Nop)) {
        ret = 1;
    } else {
        /*
         * Secondary Emission if there are any;
         */

    }
    return ret;
}
void BoundaryGeometry::callFNemission(OpalBeamline &itsOpalBeamline, PartBunch *itsBunch, const double t) {


    const double fa = 1.54 * pow(10.0, -6 + 4.52 / sqrt(workFunction_m)) / workFunction_m * fieldEnhancement_m * fieldEnhancement_m;

    for(int i = 0; i < numbfaces_global_m; i++) {

        if((TriBGphysicstag_m[i] & (BGphysics::FNEmission)) == (BGphysics::FNEmission)) {

            Vector_t E = (0.0, 0.0, 0.0), B = (0.0, 0.0, 0.0);
            Vector_t centroid = (0.0, 0.0, 0.0);
            unsigned long rtvtmp = itsOpalBeamline.getFieldAt(Tribarycent_m[i], centroid,  t, E, B);
            double Enormal = dot(TriNormal_m[i], E);

            if(Enormal < -3.0e7) {

                double Jtmp = fa * Enormal * Enormal * exp(-6.53 * pow(10.0, 9.0) * pow(workFunction_m, 1.5) / fieldEnhancement_m / (-Enormal));
                double chargeScalFactor = 1.0;
                size_t N = static_cast<size_t>(Jtmp * Triarea_m[i] / -itsBunch->getChargePerParticle() * itsBunch->getdT());
                if(N > maxFNemission_m) {
                    chargeScalFactor = static_cast<double>(N / maxFNemission_m);
                    N = maxFNemission_m;
                    cout << Enormal << " " << "Emitted Particle Number: " << N << " " << chargeScalFactor/*Jtmp *Triarea_m[i]*/ << itsBunch->getChargePerParticle() << endl;
                }

                int pc = 0;
                size_t count = itsBunch->getLocalNum();

                for(size_t k = 0; k < N; k++) {

                    if(pc == Ippl::myNode()) {
                        double r1 = IpplRandom();
                        double r2 = IpplRandom();

                        //double rsqrt = sqrt(r1*r1+r2*r2);
                        itsBunch->create(1);
                        itsBunch->R[count] = geo3Dcoords_m[allbfaces_m[4*i+1]] + 0.5 * r1 * (geo3Dcoords_m[allbfaces_m[4*i+2]] - geo3Dcoords_m[allbfaces_m[4*i+1]]) + 0.5 * r2 * (geo3Dcoords_m[allbfaces_m[4*i+3]] - geo3Dcoords_m[allbfaces_m[4*i+1]]) ; //fix me
                        itsBunch->P[count] = Vector_t(0.0);
                        itsBunch->Bin[count] = 0;
                        itsBunch->PType[count] = 1;
                        itsBunch->Q[count] = chargeScalFactor * itsBunch->getChargePerParticle();
                        itsBunch->LastSection[count] = 0;
                        itsBunch->Ef[count] = Vector_t(0.0); /// Initialize fields to zero.
                        itsBunch->Bf[count] = Vector_t(0.0);
                        itsBunch->dt[count] = itsBunch->getdT();
                        count++;

                    }

                    pc++;

                    if(pc == Ippl::getNodes())
                        pc = 0;

                }

                itsBunch->boundp();
                itsBunch->LastSection = 0;
                itsBunch->update();

            }

        }

    }
    /*    int n = 100;
          double Thre_h = 1.0e7;
          double Thre_l = -1.0e7;
          for(int i = 0; i <  n; i ++) {
          Vector_t E(0.0), B(0.0);

          int k = static_cast<int>(IpplRandom() * numbfaces_global_m) ;
          if(i < 10) {
          cout << Tribarycent_m[k] << " " << TriBGphysicstag_m[k] << " " << ((BGphysics::Absorption) | (BGphysics::FNEmission)) << " " << ((TriBGphysicstag_m[k] & (BGphysics::FNEmission)) == (BGphysics::FNEmission)) << " " << (TriBGphysicstag_m[k] & (BGphysics::FNEmission)) << endl;
          }


          if((TriBGphysicstag_m[k] & (BGphysics::FNEmission)) == (BGphysics::FNEmission)) {
          Vector_t centroid(0.0);
          unsigned long rtvtmp = itsOpalBeamline.getFieldAt(Tribarycent_m[k], centroid,  t, E, B);
          double projE = dot(E, TriNormal_m[k]);
          // if (E!= (0.0,0.0,0.0)) {
          if(projE > 0.000001 || projE < -0.000001) {
          int pc = 0;
          size_t count = itsBunch->getLocalNum();
          size_t N = 1;
          for(int j = 0; j < N; j ++) {
          if(pc == Ippl::myNode()) {
          itsBunch->create(1);
          itsBunch->R[count] = Tribarycent_m[k] ;

          itsBunch->P[count] = Vector_t(0.0);
          itsBunch->Bin[count] = 0;
          itsBunch->PType[count] = 1;
          itsBunch->Q[count] = itsBunch->getChargePerParticle();
          itsBunch->LastSection[count] = 0;
          itsBunch->Ef[count] = Vector_t(0.0); /// Initialize fields to zero.
          itsBunch->Bf[count] = Vector_t(0.0);
          itsBunch->dt[count] = itsBunch->getdT();
          count++;
          }
          pc++;
          if(pc == Ippl::getNodes())
          pc = 0;

          }
          itsBunch->boundp();
          itsBunch->LastSection = 0;
          //itsBunch->update();
          }
          }
          }*/
}

void BoundaryGeometry::setWorkFunction(double workFunction) { workFunction_m = workFunction; }
void BoundaryGeometry::setFieldEnhancement(double fieldEnhancement) { fieldEnhancement_m = fieldEnhancement; }
void BoundaryGeometry::setMaxFN(size_t maxFNemission) { maxFNemission_m = maxFNemission; }
/**
 * Determine if a line segment  has intersects with a triangle.
 * Algorithm :http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
 * @param x0 and @param x1 define the line segment.
 * @param i is the id of the triangle.
 */

Vector_t BoundaryGeometry::LineInsTri(Vector_t x0, Vector_t x1, size_t i) {
    //    INFOMSG(" LineInsTri called "<< endl);
    Vector_t lseg = x1 - x0;
    Vector_t t0, t1, t2;
    t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
    t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
    t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
    Vector_t u = t1 - t0;
    Vector_t v = t2 - t0;
    Vector_t lt = t0 - x0;
    Vector_t n;
    n[0] = u[1] * v[2] - v[1] * u[2];
    n[1] = u[2] * v[0] - v[2] * u[0];
    n[2] = u[0] * v[1] - v[0] * u[1];

    if((n[0]*lseg[0] + n[1]*lseg[1] + n[2]*lseg[2]) == 0) {
        //     INFOMSG("Triangle parallell to line segment, return x1 "<< endl);
        return x1;
    } else {

        double rI = (n[0] * lt[0] + n[1] * lt[1] + n[2] * lt[2]) / (n[0] * lseg[0] + n[1] * lseg[1] + n[2] * lseg[2]);
        Vector_t ItSec = x0 + rI * lseg;
        Vector_t w = ItSec - t0;
        if((rI < 0) || (rI > 1)) {
            //  INFOMSG("Intesect is on the extended line, return x1 "<< endl);
            return x1;
        } else {
            double tmp1, tmp2, tmp3, tmp4, tmp5;
            tmp1 = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            tmp2 = w[0] * v[0] + w[1] * v[1] + w[2] * v[2];
            tmp3 = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
            tmp4 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
            tmp5 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            double sI = (tmp1 * tmp2 - tmp5 * tmp3) / (tmp1 * tmp1 - tmp4 * tmp5);
            double tI = (tmp1 * tmp3 - tmp4 * tmp2) / (tmp1 * tmp1 - tmp4 * tmp5);
            if((sI >= 0) && (tI >= 0) && ((sI + tI) <= 1)) {
                return ItSec;
            } else {
                //    INFOMSG("Intesect is on the extended plane, return x1 "<< endl);
                return x1;
            }
        }
    }
}



/**
 * Determine if a point is outside (return 0), inside (return 1) or just on (return 2) the boundary.
 * @param x stands for coordinates of the test point.
 * The basic idea is if a line segment starting from the test point has odd intersects with a closed boundary,
 * then the test point is inside the geometry; if the intersects have even number, then the test points is outside the geometry;
 * if the test point is amoung the intersects, then the test point is just on the boundary. Makesure the end point of the line
 * segment is outside the geometry boundary.
 */

int  BoundaryGeometry::isInside(Vector_t x) {
    Vector_t x0 = x;
    Vector_t x1;
    x1[0] = x0[0];
    x1[1] = x0[1];
    x1[2] = maxcoords_m[2] * (1 + IpplRandom()); //Random number could avoid some specific situation,like line parallel to boundary......
    // x1 could be any point outside the boundary like :x1[0] = hr_m[0]+maxcoords_m[0]; x1[1] = hr_m[1]+maxcoords_m[1]  x1[2] = x0[2] ;
    if(isInGeometry(x0)) {

        vector<Vector_t> IntesecNum = PartBoundaryInteNum(x0, x1);
        vector<Vector_t>::iterator InIt;
        for(InIt = IntesecNum.begin(); InIt != IntesecNum.end(); InIt++) {
            INFOMSG("IntesecNum " << *InIt << endl);
        }
        if(IntesecNum[0] == x0) {
            INFOMSG("IntesecNum[0] " << IntesecNum[0] << endl);
            return 2;//x0 is just on the boundary;
        } else {
            if(((IntesecNum.size() % 2) == 0) || (*IntesecNum.begin() == x1)) {
                return 0;//x0 is  outside the boundary;
            } else
                return 1;//x0 is inside the boundary;
        }
    } else
        return 1;//if x0 is in boundary cubic start to use the accurate algorithm, or the x0 is just inside the geometry.
    //(we have to make sure that x0 will not jump outside in one time step and the initial particle will not outside the geometry)
}




void  BoundaryGeometry::initialize() {
    IpplTimings::startTimer(TPreProc_m);
    string dir("vtk");
    mkdir(dir.c_str(), O_CREAT);
    chmod(dir.c_str(), S_IRWXU);

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS); //Property list identifier
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    // Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fopen(getFilename().c_str(), H5F_ACC_RDONLY, plist_id);
    assert(file_id >= 0);
    H5Pclose(plist_id);

    /////////////////////////////////////////////
    //   Read dataset "surface" from .h5 file
    ////////////////////////////////////////////

    hsize_t dimsf[2];//dataset dimensions
    herr_t  status;

    hid_t dset_id = H5Dopen(file_id, "/surface");
    assert(dset_id >= 0);
    // Create the dataspace for the dataset.
    hid_t x = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(x, dimsf, NULL);
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);

    numbfaces_global_m = dimsf[0];
    INFOMSG("total_num= " << numbfaces_global_m << endl);

    allbfaces_m = new int [numbfaces_global_m * 4];

    TriPrPartloss_m = new int [numbfaces_global_m];
    TriSePartloss_m = new int [numbfaces_global_m];

    bfaces_idx_m = new int [numbfaces_global_m];

    Tribarycent_m = new Vector_t[numbfaces_global_m];
    // Create property list for collective dataset write.
    hid_t plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id2, H5FD_MPIO_COLLECTIVE);
    status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace, plist_id2, allbfaces_m);
    assert(status >= 0);
    // Store local boundary faces, discard others
    nof_sym_m = 0; //Count number of symmetry planes.
    int const nogeosym_flag = 0; //heronion outputs a index, case 0 indicates no symmetry plane, 1 for (x,y) sym, 2 for (y,z), 3 for (z,x),none 0 means not a surface.
    for(int i = 0; i < numbfaces_global_m; i ++) {
        bfaces_idx_m[i] = i;
        if(allbfaces_m[4 * i] > nogeosym_flag) {
            nof_sym_m += 1;
            if(i < numbfaces_global_m - 1) {
                for(int j = 0; j < numbfaces_global_m - i; j ++) {
                    allbfaces_m[4*(i+j)] = allbfaces_m[4*(i+j+1)];
                    allbfaces_m[4*(i+j)+1] = allbfaces_m[4*(i+j+1)+1];
                    allbfaces_m[4*(i+j)+2] = allbfaces_m[4*(i+j+1)+2];
                    allbfaces_m[4*(i+j)+3] = allbfaces_m[4*(i+j+1)+3];
                }
            } else
                numbfaces_global_m = numbfaces_global_m - 1;
        }
    }
    H5Dclose(dset_id);
    H5Sclose(filespace);
    ////////////////////////////////
    //  Also read dataset "coords"
    ////////////////////////////////
    hsize_t dimsf_c[2];
    herr_t status_c;
    hid_t dset_id_c = H5Dopen(file_id, "/coords");
    assert(dset_id_c >= 0);

    // Create the dataspace for the dataset.
    hid_t cp_space = H5Dget_space(dset_id_c);
    H5Sget_simple_extent_dims(cp_space, dimsf_c, NULL);
    hid_t filespace_c = H5Screate_simple(2, dimsf_c, NULL);

    numpoints_global_m = dimsf_c[0];
    double *point_coords = new double[3 * numpoints_global_m];
    hid_t plist_id3 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id3, H5FD_MPIO_COLLECTIVE);
    status_c = H5Dread(dset_id_c, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_c, plist_id3, point_coords);
    assert(status_c >= 0);
    Vector_t geo3d_tmp;
    for(int i = 0; i < numpoints_global_m; i++) {
        geo3d_tmp[0] = point_coords[3*i];
        geo3d_tmp[1] = point_coords[3*i+1];
        geo3d_tmp[2] = point_coords[3*i+2];

        geo3Dcoords_m.push_back(geo3d_tmp);
    }
    INFOMSG("Points_tot_num= " << numpoints_global_m << " " << geo3Dcoords_m.size() << endl);
    for(int i = 0; i < numbfaces_global_m; i ++) {
        Tribarycent_m[i] = (geo3Dcoords_m[allbfaces_m[4*i+1]] + geo3Dcoords_m[allbfaces_m[4*i+2]] + geo3Dcoords_m[allbfaces_m[4*i+3]]) / 3.0 ;
        Triarea_m.push_back(TriangleArea(i));
    }
    if(point_coords)
        delete point_coords;
    // Close HDF5 stuff
    H5Dclose(dset_id_c);
    H5Sclose(filespace_c);

    mincoords_m = getMinExtend();

    maxcoords_m = getMaxExtend();

    // fixme

    len_m = maxcoords_m - mincoords_m;  // have to fix !

    nr_m(0) = (int)floor(len_m(0) / getMaxDimenssion() / 2.0);
    nr_m(1) = (int)floor(len_m(1) / getMaxDimenssion() / 2.0);
    nr_m(2) = (int)floor(len_m(2) / getMaxDimenssion() / 2.0);

    hr_m = len_m / nr_m;

    TriPrPartlossZ_m = new int [nr_m(2)];
    TriSePartlossZ_m = new int [nr_m(2)];

    INFOMSG("maxExtend= " << getMaxExtend() << " minExtend= " << getMinExtend()  << endl);
    INFOMSG("len= " << len_m << endl);
    INFOMSG("nr= " << nr_m << endl);
    INFOMSG("hr= " << hr_m << endl);
    makeTriNormal();
    setBGphysicstag();
    IpplTimings::stopTimer(TPreProc_m);
}


/**
 * Initialize some particles near the surface with inward momenta.

 FixMe: have to eliminate nodes which are on several triangles
*/
void  BoundaryGeometry::createParticlesOnSurface(size_t n, double darkinward,  OpalBeamline &itsOpalBeamline, PartBunch &itsBunch) {

    int pc = 0;

    const myNode = Ippl::myNode();
    const nNodes = Ippl::getNodes();
    for(int i = 0; i < n; i ++) {

        if(pc == myNode) {
            Vector_t E = (0.0, 0.0), B = (0.0, 0.0);
            while(E == (0.0, 0.0, 0.0)) {
                int k = (int)(IpplRandom() * numbfaces_global_m) ;
                Vector_t centroid = (0.0, 0.0);
                const unsigned long sindex = 0;
                unsigned long temp = itsOpalBeamline.getFieldAt(Tribarycent_m[k], centroid, itsBunch.getT(), E, B);
                if((fabs(E(0)) > 0.00001) || (fabs(E(1)) > 0.00001) || (fabs(E(2)) > 0.00001))
                    partsr_m.push_back(Tribarycent_m[k] + darkinward * TriNormal_m[k]);
            }

        }
        pc++;
        if(pc == nNodes)
            pc = 0;
    }
}

/**
 * Make the boundary set by using triangle vertex, bounding box vertex ,as well as points in triangle central lines.
 */


void BoundaryGeometry::makeBoundaryIndexSet() {
    set<size_t>::iterator bbIt;
    vector<Vector_t> BBox;
    for(int i = 0; i < numbfaces_global_m; i++) {
        vector<Vector_t> discreteTri;
        vector<Vector_t>::iterator disIt;
        for(int j = 0; j < 20; j++) { //Discrete three central line to get a more complete boundary index set.
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]] + 0.05 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+2]] + geo3Dcoords_m[allbfaces_m[4*i+3]]) - geo3Dcoords_m[allbfaces_m[4*i+1]]));
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]] + 0.05 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+3]] + geo3Dcoords_m[allbfaces_m[4*i+1]]) - geo3Dcoords_m[allbfaces_m[4*i+2]]));
            discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]] + 0.05 * (j + 1) * (0.5 * (geo3Dcoords_m[allbfaces_m[4*i+1]] + geo3Dcoords_m[allbfaces_m[4*i+2]]) - geo3Dcoords_m[allbfaces_m[4*i+3]]));
        }
        discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]]);
        discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]]);
        discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]]);
        discreteTri.push_back(Tribarycent_m[i]);
        for(disIt = discreteTri.begin(); disIt != discreteTri.end(); disIt++) {
            int id = f(*disIt);

            //  assert(id > 0);
            bbIt = boundary_ids_m.find(id);
            // bbIt = boundary_ids_m.rbegin();
            if((bbIt == boundary_ids_m.end())) {
                Vector_t temp;
                temp[0] = floor(((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];
                temp[1] = floor(((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor(((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back(temp);
                boundary_ids_m.insert(id);
            }


            vector<size_t> tmp;
            map< size_t, vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find(id);
            if(It == CubicLookupTable_m.end()) {
                tmp.push_back(i);
                CubicLookupTable_m.insert(pair <size_t, vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back(i);
        }
        vector<Vector_t> ret, cubic_coords;
        ret = SetMinMaxBound(i);
        cubic_coords.push_back(ret[0]);
        cubic_coords.push_back((ret[0](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back((ret[0](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back((ret[0](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back((ret[1](0), ret[1](1), ret[0](2)));
        cubic_coords.push_back((ret[1](0), ret[1](1), ret[1](2)));
        cubic_coords.push_back((ret[1](0), ret[0](1), ret[1](2)));
        cubic_coords.push_back((ret[1](0), ret[0](1), ret[0](2)));
        // vector<Vector_t>::iterator disIt;
        for(disIt = cubic_coords.begin(); disIt != cubic_coords.end(); disIt++) {
            int id = f(*disIt);
            bbIt = boundary_ids_m.find(id);

            //   assert(id > 0);


            if((bbIt == boundary_ids_m.end())) {
                Vector_t temp;
                temp[0] = floor(((*disIt)[0] - mincoords_m[0]) / hr_m[0]) * hr_m[0] + mincoords_m[0];
                temp[1] = floor(((*disIt)[1] - mincoords_m[1]) / hr_m[1]) * hr_m[1] + mincoords_m[1];
                temp[2] = floor(((*disIt)[2] - mincoords_m[2]) / hr_m[2]) * hr_m[2] + mincoords_m[2];
                BBox.push_back(temp);
                boundary_ids_m.insert(id);
            }

            vector<size_t> tmp;
            map< size_t, vector<size_t> >::iterator It;
            It =  CubicLookupTable_m.find(id);
            if(It == CubicLookupTable_m.end()) {
                tmp.push_back(i);
                CubicLookupTable_m.insert(pair <size_t, vector<size_t> > (id, tmp));
            } else
                (*It).second.push_back(i);
        }

    }
    INFOMSG("size of boundary index set " << boundary_ids_m.size() << BBox.size() << endl);
    INFOMSG("size of all index " << nr_m(0)*nr_m(1)*nr_m(2) << endl);
    INFOMSG("size of geo3Dcoords " << geo3Dcoords_m.size() << endl);

    /*-------------------------------------------------------------------------*/
    size_t numpoints = 8 * BBox.size();
    vector<Vector_t>::iterator bIt;
    std::ofstream of;
    of.open(string("vtk/testBBox.vtk").c_str());
    assert(of.is_open());
    of.precision(6);
    of << "# vtk DataFile Version 2.0" << endl;                      // first line of a vtk file
    of << "generated using BoundaryGeometry::makeBoundaryIndexSet" << endl;   // header
    of << "ASCII" << endl << endl;                                   // file format
    of << "DATASET UNSTRUCTURED_GRID" << endl;                 // dataset structrue
    of << "POINTS " << numpoints << " float" << endl;  // data num and type
    for(bIt = BBox.begin(); bIt != BBox.end() ; bIt ++) {
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] << " " << (*bIt)[1] + hr_m[1]  << " " << (*bIt)[2] + hr_m[2] << endl;
        of << (*bIt)[0] + hr_m[0] << " " << (*bIt)[1] + hr_m[1] << " " << (*bIt)[2] + hr_m[2] << endl;
    }
    of << endl;

    of << "CELLS " << BBox.size() << " " << 9 * BBox.size() << endl;
    for(int i = 0; i < BBox.size() ; i ++)
        of << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3 << " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << endl;
    of << "CELL_TYPES " << BBox.size() << endl;
    for(int i = 0; i <  BBox.size() ; i ++)
        of << "11" << endl;
    of << "CELL_DATA " << BBox.size() << endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << endl;
    of << "LOOKUP_TABLE " << "default" << endl ;
    for(int i = 0; i <  BBox.size() ; i ++)
        of << (float)(i) << endl;
    of << endl;
    of << "COLOR_SCALARS " << "BBoxColor " << 4 << endl ;
    for(int i = 0; i < BBox.size() ; i ++) {

        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << endl;

    }

    of << endl;
}

/**
 * Determine whether a particle with position @param r, momenta @param v, and time step @param dt will hit the boundary.
 *
 * Basic algorithms are as follows:
 * 1) Determine if the particle is near the boundary by checking r is in boundary bounding box index set.
 *    if r is in, then do the following checking step, else return 1 to indicate that particle is far from boundary and to be integrated.
 * 2) Traversal all the triangles in the bounding cubic box which the particle is in, as well as triangles in the adjacent 26 bounding cubic boxes
 *    if the momenta has oppsite direction with those triangles' normals, then check if the particle has intersection with those triangles. If
 *    intersection exsists, then return 0.
 */

int BoundaryGeometry::PartInside(Vector_t r, Vector_t v, double dt, short Parttype, Vector_t &intecoords, int &triId) {

    int ret;
    double betaP = 1.0 / sqrt(1.0 + Dotproduct(v, v));
    Vector_t temp1 = r;//particle position in tstep n;
    Vector_t temp = r + (c * betaP * v * dt); //particle position in tstep n+1;

    IpplTimings::startTimer(TPInside_m);

    //test if particle position in tstep n is inside the cubic bounding box. If true, do the following tests
    if(isInGeometry(temp1)) {
        int flg = 0;
        int id = f(temp1);

        // Build an array contains the ids of 27(3*3*3) cubic boxes. The id of the cubic box which contains the particle,
        // lies in the center of the 27 cubic boxes. Just for the situation that even the line segment cross more than one cubic box.
        int idc[27] = {id, id + 1, id - 1, id + nr_m(0), id - nr_m(0), id + nr_m(0) + 1, id + nr_m(0) - 1, id - nr_m(0) - 1, id - nr_m(0) + 1,
                       id + nr_m(0) *nr_m(1), id + nr_m(0) *nr_m(1) + 1, id + nr_m(0) *nr_m(1) - 1, id + nr_m(0) *nr_m(1) + nr_m(0),
                       id + nr_m(0) *nr_m(1) - nr_m(0), id + nr_m(0) *nr_m(1) + nr_m(0) + 1, id + nr_m(0) *nr_m(1) + nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1), id - nr_m(0) *nr_m(1) + 1, id - nr_m(0) *nr_m(1) - 1, id - nr_m(0) *nr_m(1) + nr_m(0), id - nr_m(0) *nr_m(1) - nr_m(0), id - nr_m(0) *nr_m(1) + nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1) + nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) + 1
                      };

        // Test all the 27 cubic boxes to find if the line segment has intersection with the triangles in those cubic boxes.
        for(int k = 0; k < 27; k++) {
            // If cubic box among the 27 cubic boxes is not the boundary bounding box, do not need to be tested at all. If cubic box is the boundary bounding box then do the following tests.
            if(boundary_ids_m.find(idc[k]) != boundary_ids_m.end()) {
                map< size_t, vector<size_t> >::iterator It;
                vector<size_t> ::iterator myIt;
                It = CubicLookupTable_m.find(idc[k]);
                // Lookup table (stl map) contains the id of cubic box as key and triangle ids inside the id th cubic box as value
                for(myIt = (*It).second.begin(); myIt != (*It).second.end(); myIt++) {
                    // For each particle inside the cubic box, test if this particle have intersect with each triangle inside the same cubic box
                    if((v != 0) && (Dotproduct(v, TriNormal_m[*myIt]) <= 0.0)) {
                        // If the particle have none zero momenta and momenta has opposite direction with triangle normal, do the following tests.
                        IpplTimings::startTimer(TRayTrace_m);
                        Intersection tmp = FindIntersection(temp1, temp, *myIt); //Ray tiangle test, parameters are particle position in tstep n, particle position in tstep n+1 and tiangle id.
                        IpplTimings::stopTimer(TRayTrace_m);
                        if((tmp.Isc != maxcoords_m)) {
                            // Test if the intersection is between the particle position in tstep n and particle position in tstep n+1 or is in the extension of line
                            // segment when particle position in tstep n is already outside the geometry( this may be not accurate and may be the source of problem.)
                            if(((tmp.rI >= -0.00001) && (tmp.rI <= 1.00001))  || ((tmp.rI < 0) && (Dotproduct(temp1 - Tribarycent_m[*myIt], TriNormal_m[*myIt]) <= 0))) {
                                flg++;
                                intecoords = tmp.Isc;
                                triId = *myIt;
                                if(Parttype == 1)
                                    TriSePartloss_m[*myIt]++;
                                else
                                    TriPrPartloss_m[*myIt]++;

                                break;
                            }
                        }
                    }
                }
                if(flg != 0) {
                    ret = 0;
                    break;
                }
            }
        }
        if(flg == 0) {
            ret = 1;
        }
    } else if(isInGeometry(temp)) {
        int flg = 0;
        int id = f(temp);

        // Build an array contains the ids of 27(3*3*3) cubic boxes. The id of the cubic box which contains the particle,
        // lies in the center of the 27 cubic boxes. Just for the situation that even the line segment cross more than one cubic box.
        int idc[27] = {id, id + 1, id - 1, id + nr_m(0), id - nr_m(0), id + nr_m(0) + 1, id + nr_m(0) - 1, id - nr_m(0) - 1, id - nr_m(0) + 1,
                       id + nr_m(0) *nr_m(1), id + nr_m(0) *nr_m(1) + 1, id + nr_m(0) *nr_m(1) - 1, id + nr_m(0) *nr_m(1) + nr_m(0),
                       id + nr_m(0) *nr_m(1) - nr_m(0), id + nr_m(0) *nr_m(1) + nr_m(0) + 1, id + nr_m(0) *nr_m(1) + nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) - 1, id + nr_m(0) *nr_m(1) - nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1), id - nr_m(0) *nr_m(1) + 1, id - nr_m(0) *nr_m(1) - 1, id - nr_m(0) *nr_m(1) + nr_m(0), id - nr_m(0) *nr_m(1) - nr_m(0), id - nr_m(0) *nr_m(1) + nr_m(0) + 1,
                       id - nr_m(0) *nr_m(1) + nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) - 1, id - nr_m(0) *nr_m(1) - nr_m(0) + 1
                      };

        // Test all the 27 cubic boxes to find if the line segment has intersection with the triangles in those cubic boxes.
        for(int k = 0; k < 27; k++) {
            // If cubic box among the 27 cubic boxes is not the boundary bounding box, do not need to be tested at all. If cubic box is the boundary bounding box then do the following tests.
            if(boundary_ids_m.find(idc[k]) != boundary_ids_m.end()) {
                map< size_t, vector<size_t> >::iterator It;
                vector<size_t> ::iterator myIt;
                It = CubicLookupTable_m.find(idc[k]);
                // Lookup table (stl map) contains the id of cubic box as key and triangle ids inside the id th cubic box as value
                for(myIt = (*It).second.begin(); myIt != (*It).second.end(); myIt++) {
                    // For each particle inside the cubic box, test if this particle have intersect with each triangle inside the same cubic box
                    if((v != 0) && (Dotproduct(v, TriNormal_m[*myIt]) <= 0.0)) {
                        // If the particle have none zero momenta and momenta has opposite direction with triangle normal, do the following tests.
                        IpplTimings::startTimer(TRayTrace_m);
                        Intersection tmp = FindIntersection(temp1, temp, *myIt); //Ray tiangle test, parameters are particle position in tstep n, particle position in tstep n+1 and tiangle id.
                        IpplTimings::stopTimer(TRayTrace_m);
                        if((tmp.Isc != maxcoords_m)) {
                            // Test if the intersection is between the particle position in tstep n and particle position in tstep n+1 or is in the extension of line
                            // segment when particle position in tstep n is already outside the geometry( this may be not accurate and may be the source of problem.)
                            if(((tmp.rI >= -0.00001) && (tmp.rI <= 1.00001))  || ((tmp.rI < 0) && (Dotproduct(temp1 - Tribarycent_m[*myIt], TriNormal_m[*myIt]) <= 0))) {
                                flg++;
                                intecoords = tmp.Isc;
                                triId = *myIt;
                                if(Parttype == 1)
                                    TriSePartloss_m[*myIt]++;
                                else
                                    TriPrPartloss_m[*myIt]++;

                                break;
                            }
                        }
                    }
                }
                if(flg != 0) {
                    ret = 0;
                    break;
                }
            }
        }
        if(flg == 0) {
            ret = 1;
        }
    } else {
        ret = 1;
    }
    IpplTimings::stopTimer(TPInside_m);
    return ret;
}



void BoundaryGeometry::updateElement(ElementBase *element) {

}



Inform &BoundaryGeometry::print(Inform &os) const {
    os << "* ************* W A K E ************************************************************ " << endl;
    os << "* GEOMETRY     " << getOpalName() << '\n'
       << "* FGEOM        " << Attributes::getString(itsAttr[FGEOM]) << '\n'
       << "* TOPO         " << Attributes::getString(itsAttr[TOPO]) << '\n'
       << "* LENGTH       " << Attributes::getReal(itsAttr[LENGTH]) << '\n'
       << "* S            " << Attributes::getReal(itsAttr[S]) << '\n'
       << "* A            " << Attributes::getReal(itsAttr[A]) << '\n'
       << "* B            " << Attributes::getReal(itsAttr[B]) << endl ;
    os << "* ********************************************************************************** " << endl;








}
