/**
 * @file MagneticField.h
 * This class provides an interpolation algorithm for the computation of \f$ B(r,\theta),\ \frac{\partial B(r,\theta)}{\partial r},\
 * \frac{\partial B(r,\theta)}{\partial\theta} \f$
 *
 * @author Dr. Christian Baumgarten
 * @version 1.0
 */
#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H

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

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "Utilities/GeneralClassicException.h"

extern Inform *gmsg;

#define CHECK_CYC_FSCANF_EOF(arg) if(arg == EOF)\
throw GeneralClassicException("Cyclotron::getFieldFromFile",\
                              "fscanf returned EOF at " #arg);


template<typename Value_type>
class MagneticField
{
public:
    struct BPositions {
        // this 4 parameters are need to be read from field file.
        double  rmin, delr;
        double  tetmin, dtet;
    
        // Radii and step width of initial Grid
        std::vector<double> rarr;
    
        //  int     ThetaPeriodicity; // Periodicity of Magnetic field
        double  Bfact;      // MULTIPLICATION FACTOR FOR MAGNETIC FIELD
    };
    
    struct BfieldData {
    std::string filename;
    // known from file: field and three theta derivatives
    //~ double *bfld;   //Bz
    //~ double *dbt;    //dBz/dtheta
    //~ double *dbtt;   //d2Bz/dtheta2
    //~ double *dbttt;  //d3Bz/dtheta3
//~
    //~ // to be calculated in getdiffs: all other derivatives:
    //~ double *dbr;    // dBz/dr
    //~ double *dbrr;   // ...
    //~ double *dbrrr;
//~
    //~ double *dbrt;
    //~ double *dbrrt;
    //~ double *dbrtt;
//~
    //~ // used to get (Br,Btheta,Bz) at any off-plane point
    //~ double *f2;  // for Bz
    //~ double *f3;  // for Br
    //~ double *g3;  // for Btheta
//~
    std::vector<double> bfld;   //Bz
    std::vector<double> dbt;    //dBz/dtheta
    std::vector<double> dbtt;   //d2Bz/dtheta2
    std::vector<double> dbttt;  //d3Bz/dtheta3

    // to be calculated in getdiffs: all other derivatives:
    std::vector<double> dbr;    // dBz/dr
    std::vector<double> dbrr;   // ...
    std::vector<double> dbrrr;

    std::vector<double> dbrt;
    std::vector<double> dbrrt;
    std::vector<double> dbrtt;

    // used to get (Br,Btheta,Bz) at any off-plane point
    std::vector<double> f2;  // for Bz
    std::vector<double> f3;  // for Br
    std::vector<double> g3;  // for Btheta

    // Grid-Size
    //need to be read from inputfile.
    int nrad, ntet;

    // one more grid line is stored in azimuthal direction:
    int ntetS;

    // total grid points number.
    int ntot;

    // Mean and Maximas
    double bacc, dbtmx, dbttmx, dbtttmx;
};
    
    MagneticField(const std::string fmapfn);
    
    /*!
     * Read in a CARBONCYCL fieldmap.
     */
    void read(const Value_type &scaleFactor);
    
    /*!
     * @param bint stores the interpolated magnetic field
     * @param brint stores the interpolated radial derivative of the magnetic field
     * @param btint stores the interpolated "theta" derivative of the magnetic field
     * @param r is radial position where to interpolate
     * @param theta is the azimuth angle where to interpolate
     */
    void interpolate(Value_type& bint, Value_type& brint, Value_type& btint,
                     const Value_type& r, Value_type theta);
  
  
    ~MagneticField() {
        gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
    }
    
    inline int getNradSteps();
    
    inline int getNtetSteps();
  
private:
    
  // object of Matrics including magnetic field map and its derivates
  BfieldData Bfield;
  
  // object of parameters about the map grid
  BPositions BP;
  
  inline int idx(int irad, int ktet) {return (ktet + Bfield.ntetS * irad);}
  
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  
  std::string fmapfn_m;     ///< Fieldmap file
  
};


template<typename Value_type>
MagneticField<Value_type>::MagneticField(const std::string fmapfn) :
    fmapfn_m(fmapfn)
{ }

template<typename Value_type>
void MagneticField<Value_type>::read(const Value_type& scaleFactor) {
    
    FILE *f = NULL;
    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*      READ IN CARBON CYCLOTRON FIELD MAP       " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Cyclotron::getFieldFromFile_Carbon",
                                      "failed to open file '" + fmapfn_m + "', please check if it exists");
    }

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.delr));
    //if the value is negative, the actual value is its reciprocal.
    if(BP.delr < 0.0) BP.delr = 1.0 / (-BP.delr);
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.dtet));
    //if the value is negative, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuthal direction: " << BP.dtet << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Grid points along azimuth (ntet): " << Bfield.ntet << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Grid points along radius (nrad): " << Bfield.nrad << endl;

    //Bfield.ntetS = Bfield.ntet;
    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    //Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;
    Bfield.ntot = Bfield.nrad * Bfield.ntetS;
    
    *gmsg << "* Adding a guard cell along azimuth" << endl;
    *gmsg << "* Total stored grid point number ((ntet+1) * nrad) : " << Bfield.ntot << endl;
    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    std::vector<Value_type> r(Bfield.nrad), theta(Bfield.ntet), bval(Bfield.ntot);
    for(int i = 0; i < Bfield.nrad; i++) {
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] *= BP.Bfact;
        }
    }
    
    std::fstream fp1;
    fp1.open("gnu.out", std::ios::out);
    for (int i = 0; i < Bfield.nrad; ++i) {
        for (int k = 0; k < Bfield.ntet; ++k) {
            fp1 << BP.rmin + (i * BP.delr) << " \t " << k * (BP.tetmin + BP.dtet)
                << " \t " << Bfield.bfld[idx(i, k)] << std::endl;
        }
    }
    fp1.close();
    
    const gsl_interp2d_type *T = gsl_interp2d_bicubic; //gsl_interp2d_bilinear;
    spline = gsl_spline2d_alloc(T, Bfield.nrad, Bfield.ntet);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
    
    for (int i = 0; i < Bfield.ntet; ++i)
        theta[i] = i * (BP.tetmin + BP.dtet);
    
    for (int i = 0; i < Bfield.nrad; ++i) {
        r[i] = BP.rmin * 0.001 + (i * BP.delr * 0.001);
        for (int k = 0; k < Bfield.ntet; ++k) {
            bval[Bfield.nrad * k + i] = Bfield.bfld[idx(i, k)];
        }
    }
    
    gsl_spline2d_init(spline, &r[0], &theta[0], &bval[0], r.size(), theta.size());
    
    fclose(f);

    *gmsg << "* Field Maps read successfully!" << endl << endl;
}

template<typename Value_type>
void MagneticField<Value_type>::interpolate(Value_type& bint,
                                            Value_type& brint,
                                            Value_type& btint,
                                            const Value_type &r,
                                            Value_type theta) {

    bint = 0;
    brint = 0;
    btint = 0;
    
    // return if we are out of the radius range
    if ( r < BP.rmin * 0.001 || r > (BP.rmin + (Bfield.nrad - 1)* BP.delr) * 0.001 )
        return;
    
    // scale the angle to the correct range
    if ( theta < 0 )
        theta = 360.0 + theta;
    
    if ( theta > BP.tetmin + (Bfield.ntet - 1) * BP.dtet ) {
        Value_type max = BP.tetmin + (Bfield.ntet - 1) * BP.dtet;
        theta -= int(theta / max) * max;
    }
    
    bint = gsl_spline2d_eval(spline, r, theta, xacc, yacc);
    
    // derivative w.r.t. the radius, i.e. dB/dr
    brint = gsl_spline2d_eval_deriv_x(spline, r, theta, xacc, yacc);
    
    // derivative w.r.t. the azimuth angle, i.e. dB/dtheta
    btint = gsl_spline2d_eval_deriv_y(spline, r, theta, xacc, yacc);
}

template<typename Value_type>
int MagneticField<Value_type>::getNradSteps() {
    return Bfield.nrad;
}

template<typename Value_type>
int MagneticField<Value_type>::getNtetSteps() {
    return Bfield.ntet;
}

#endif
