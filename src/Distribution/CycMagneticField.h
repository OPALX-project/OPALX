#ifndef CYCMAGNETICFIELD_H
#define CYCMAGNETICFIELD_H

#include <cmath>

#include "AbsBeamline/Cyclotron.h"


class CycMagneticField : public Cyclotron {
    
public:
    
    void interpolate(double& bint,
                     double& brint,
                     double& btint,
                     const double& r,
                     const double& theta);
    
private:
    bool interpolate_m(const double& r,
                       const double& theta,
                       const double& z,
                       double& bint
                      );
    
};


void CycMagneticField::interpolate(double& bint,
                                   double& brint,
                                   double& btint,
                                   const double& r,
                                   const double& theta,
                                   const double& z
                                  )
{
    
    interpolate_m(r, theta, z, bint);
    
    
    // derivative w.r.t. the radius
    double dr = 0.0001
    
    
    double b1, b2;
    
    double tmp = r + dr;
    interpolate_m(tmp, theta, z, b1);
    
    tmp = r - dr;
    interpolate_m(tmp, theta, z, b2);
    
    brint = (b2 - b1) / (2.0 * dr);
    
    tmp = theta + dr;
    
    interpolate_m(r, tmp, z, b1);
    
    tmp = theta - dr;
    
    interpolate_m(r, tmp, z, b2);
    
    btint = (b2 - b1) / (2.0 * dr);
}

bool CycMagneticField::interpolate_m(const double& rad,
                                     const double& theta,
                                     const double& z,
                                     double& bint
                                    )
{
    // x horizontal
    // y longitudinal
    // z is vertical
    
    const double xir = (rad - BP.rmin) / BP.delr;

    // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
    const int    ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;
    

//     double tet_rad = theta;

    // the actual angle of particle
    theta = theta / pi * 180.0;
    
    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    tet_map = fmod(theta, 360.0 / symmetry_m);

    xit = tet_map / BP.dtet;

    int it = (int) xit;

    //    *gmsg << R << " tet_map= " << tet_map << " ir= " << ir << " it= " << it << " bf= " << Bfield.bfld[idx(ir,it)] << endl;

    const double wt1 = xit - (double)it;
    const double wt2 = 1.0 - wt1;

    // it : the number of point on the inner path whoes angle is less then the particle' corresponding angle.
    // include zero degree point
    it = it + 1;

    int r1t1, r2t1, r1t2, r2t2;
    int ntetS = Bfield.ntet + 1;

    // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
    // considering  the array start with index of zero, minus 1.

    if(myBFieldType_m != FFAGBF) {
        /*
          For FFAG this does not work
        */
        r1t1 = it + ntetS * ir - 1;
        r1t2 = r1t1 + 1;
        r2t1 = r1t1 + ntetS;
        r2t2 = r2t1 + 1 ;

    } else {
        /*
          With this we have B-field AND this is far more
          intuitive for me ....
        */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
    }

    double brf = 0.0
    bint = 0.0;

    if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)) {

        /* Br */
        brf = (Bfield.dbr[r1t1] * wr2 * wt2 +
               Bfield.dbr[r2t1] * wr1 * wt2 +
               Bfield.dbr[r1t2] * wr2 * wt1 +
               Bfield.dbr[r2t2] * wr1 * wt1) * z;

        bint = - brf;
        
    } else {
      return true;
    }
    
    return false;
}

#endif