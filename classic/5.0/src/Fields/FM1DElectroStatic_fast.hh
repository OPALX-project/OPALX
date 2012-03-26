#ifndef CLASSIC_FIELDMAP1DELECTROSTATICFAST_HH
#define CLASSIC_FIELDMAP1DELECTROSTATICFAST_HH

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "Fields/Fieldmap.hh"

class FM1DElectroStatic_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FM1DElectroStatic_fast(std::string aFilename);
    ~FM1DElectroStatic_fast();

    virtual void readMap();
    virtual void freeMap();

    double *onAxisField_m;
    gsl_spline *onAxisInterpolants_m[4];
    gsl_interp_accel *onAxisAccel_m[4];

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double length_m;
    double hz_m;

    int num_gridpz_m;

    friend class Fieldmap;
};

#endif
