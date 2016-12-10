#ifndef CYCMAGNETICFIELD_H
#define CYCMAGNETICFIELD_H


#include <iostream>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "Utilities/GeneralClassicException.h"
#include "AbsBeamline/Cyclotron.h"

#include "Physics/Physics.h"


class CycMagneticField {
    
public:
    /*!
     * @param fmapfn specifies the fieldmap file. We support only
     * the CARBONCYCL format at the moment.
     */
    CycMagneticField(const std::string fmapfn, const double& symmetry);
    
    /*!
     * @param scaleFactor for the magnetic field.
     */
    void read(const double& scaleFactor);
    
    /*!
     * Get the value of B_{z}, dB_{z}/dr and dB_{z}/dtheta
     * where theta is the azimuth angle and r the radius.
     * @param bint is the magnetic field B_{z} [kG]
     * @param brint is dB_{z}/dr [kG/m]
     * @param btint is dB_{z}/dtheta [kG/deg]
     * @param r is the radial interpolation point [m]
     * @param theta is the azimuthal interpolation point [deg]
     */
    void interpolate(double& bint,
                     double& brint,
                     double& btint,
                     const double& r,
                     const double& theta);
private:
    /*!
     * @param rmin is the minimum radius [mm]
     * @param dr is the radial spacing [mm]
     * @param nrad is the number of radial points
     */
    void initR_m(double rmin, double dr, int nrad);
    
    /*!
     * Interpolation of the derivatives of the magnetic field using
     * the 5-point Lagrange's formula.
     * 
     * @param f is the startaddress for the 6 support points
     * @param dx is the stepwidth
     * @param kor is the order of the derivative (1, 2 or 3)
     * @param krl is the number of support points, where the derivative is to be
     *        calculated
     * @param lpr is the distance of the 5 storage posiitions (=1 if they are
     *        neighbors or length of columnlength of a matrix, if the support
     *        points are on a line)
     * @return the interpolated value of the derivative of the magnetic field.
     */
    double gutdf5d_m(double *f, double dx, const int kor, const int krl, const int lpr);
    
    /*!
     * Initialize the derivatives of the magnetic field on
     * the grid
     */
    void getdiffs_m();
    
private:
    
    BfieldData Bfield;      ///< object of Matrics including magnetic field map and its derivates
    
    BPositions BP;          ///< object of parameters about the map grid
    
    inline int idx(int irad, int ktet) {return (ktet + Bfield.ntetS * irad);}
    
    std::string fmapfn_m;   ///< Fieldmap file
    
    double symmetry_m;      ///< Number of sectors
};

#endif