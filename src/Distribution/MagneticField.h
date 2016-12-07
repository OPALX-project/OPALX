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

// #include "Utilities/LogicalError.h"

// #include <Ippl.h>

/// Memory allocator
#define MALLOC(n,a) ((a*)malloc((n)*sizeof(a)))

#define CHECK_CYC_FSCANF_EOF(arg) if(arg == EOF)\
throw GeneralClassicException("Cyclotron::getFieldFromFile",\
                              "fscanf returned EOF at " #arg);


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
                     const Value_type& r, const Value_type& theta);
  
  
    ~MagneticField() {
        gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
    }
  
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




#endif
