#include "CycMagneticField.h"

CycMagneticField::CycMagneticField(const std::string fmapfn,
                                   const double& symmetry) :
    Cyclotron(),
    nullGeom_m(NullGeometry()),
    nullField_m(NullField())
{
    this->setFieldMapFN(fmapfn);
    this->setSymmetry(symmetry);
}

void CycMagneticField::read(const std::string& type,
                            const double& scaleFactor) {
    
    if ( type == "CARBONCYCL" )
        this->getFieldFromFile_Carbon(scaleFactor);
    else if ( type == "CYCIAE" )
        this->getFieldFromFile_CYCIAE(scaleFactor);
    else if ( type == "AVFEQ" )
        this->getFieldFromFile_AVFEQ(scaleFactor);
    else if ( type == "FFAG" )
        this->getFieldFromFile_FFAG(scaleFactor);
    else if ( type == "BANDRF" )
        this->getFieldFromFile_BandRF(scaleFactor);
    else if ( type == "SYNCHROCYCLOTRON" )
        this->getFieldFromFile_Synchrocyclotron(scaleFactor);
    else
        this->getFieldFromFile(scaleFactor);
    
    if ( BP.rarr.empty() ) {
        // calculate the radii of initial grid.
        this->initR(BP.rmin, BP.delr, Bfield.nrad);
    }
    
    if ( Bfield.dbr.empty() || Bfield.dbt.empty() ) {
        // calculate the remaining derivatives
        this->getdiffs();
    }
}


void CycMagneticField::interpolate(double& bint,
                                   double& brint,
                                   double& btint,
                                   const double& rad,
                                   const double& theta
                                  )
{
    // x horizontal
    // y longitudinal
    // z is vertical
    const double xir = (rad * 1000 - BP.rmin) / (BP.delr);

    // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
    const int    ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;
    

    double tet_rad = theta;

    // the actual angle of particle
//     tet_rad = theta / Physics::pi * 180.0;
    
    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    double tet_map = fmod(tet_rad, 360.0 / this->getSymmetry());

    double xit = tet_map / BP.dtet;

    int it = (int) xit;

    const double wt1 = xit - (double)it;
    const double wt2 = 1.0 - wt1;

    // it : the number of point on the inner path whoes angle is less then the particle' corresponding angle.
    // include zero degree point
    it = it + 1;

    int r1t1, r2t1, r1t2, r2t2;
//     int ntetS = Bfield.ntet + 1;

    // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
    // considering  the array start with index of zero, minus 1.

//     if(myBFieldType_m != FFAGBF) {
        /*
          For FFAG this does not work
        */
//         r1t1 = it + ntetS * ir - 1;
//         r1t2 = r1t1 + 1;
//         r2t1 = r1t1 + ntetS;
//         r2t2 = r2t1 + 1 ;

//     } else {
//         /*
//           With this we have B-field AND this is far more
//           intuitive for me ....
//         */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
//     }

    bint = 0.0;
    brint = 0.0;
    btint = 0.0;

    if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)) {
        
        // dB_{z}/dr
        brint = (Bfield.dbr[r1t1] * wr2 * wt2 +
                Bfield.dbr[r2t1] * wr1 * wt2 +
                Bfield.dbr[r1t2] * wr2 * wt1 +
                Bfield.dbr[r2t2] * wr1 * wt1) * 1000.0;
        
        // dB_{z}/dtheta
        btint = Bfield.dbt[r1t1] * wr2 * wt2 +
                Bfield.dbt[r2t1] * wr1 * wt2 +
                Bfield.dbt[r1t2] * wr2 * wt1 +
                Bfield.dbt[r2t2] * wr1 * wt1;
        
        // B_{z}
        bint = Bfield.bfld[r1t1] * wr2 * wt2 +
               Bfield.bfld[r2t1] * wr1 * wt2 +
               Bfield.bfld[r1t2] * wr2 * wt1 +
               Bfield.bfld[r2t2] * wr1 * wt1;
    }
}

double CycMagneticField::getSlices() const {
    return -1.0;
}

double CycMagneticField::getStepsize() const {
    return -1.0;
}

const EMField &CycMagneticField::getField() const {
    return nullField_m;
}

EMField &CycMagneticField::getField() {
    return nullField_m;
}

ElementBase* CycMagneticField::clone() const {
    return nullptr;
}

const BGeometryBase &CycMagneticField::getGeometry() const {
    return nullGeom_m;
}

BGeometryBase &CycMagneticField::getGeometry() {
    return nullGeom_m;
}