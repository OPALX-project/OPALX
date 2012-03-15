#ifndef CLASSIC_FIELDMAP2DMAGNETOSTATIC_CSPLINE_HH
#define CLASSIC_FIELDMAP2DMAGNETOSTATIC_CSPLINE_HH

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "Fields/Fieldmap.hh"

class FM2DMagnetoStatic_cspline: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FM2DMagnetoStatic_cspline(std::string aFilename);
    ~FM2DMagnetoStatic_cspline();

    virtual void readMap();
    virtual void freeMap();

    gsl_spline **Bz_interpolants_m;
    gsl_spline **Br_interpolants_m;

    gsl_interp_accel **Bz_accel_m;
    gsl_interp_accel **Br_accel_m;

    double *Bz_values_m;
    double *Br_values_m;

    double *zvals_m;
    double *rvals_m;

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    double hr_m;                   /**< length between points in grid, r-direction, m*/
    int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

#endif
