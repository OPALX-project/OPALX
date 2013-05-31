#ifndef CLASSIC_FIELDMAP2DDYNAMIC_CSPLINE_HH
#define CLASSIC_FIELDMAP2DDYNAMIC_CSPLINE_HH

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "Fields/Fieldmap.hh"

class FM2DDynamic_cspline: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FM2DDynamic_cspline(std::string aFilename);
    ~FM2DDynamic_cspline();

    virtual void readMap();
    virtual void freeMap();

    gsl_spline **Ez_interpolants_m;
    gsl_spline **Er_interpolants_m;
    gsl_spline **Ht_interpolants_m;

    gsl_interp_accel **Ez_accel_m;
    gsl_interp_accel **Er_accel_m;
    gsl_interp_accel **Ht_accel_m;

    double *Ez_values_m;
    double *Er_values_m;
    double *Ht_values_m;

    double *zvals_m;
    double *rvals_m;

    double frequency_m;

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
