#ifndef CLASSIC_FIELDMAP1DDYNAMICFAST_HH
#define CLASSIC_FIELDMAP1DDYNAMICFAST_HH

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DDynamic_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(vector<pair<double, double> > & F);

private:
    FM1DDynamic_fast(string aFilename);
    ~FM1DDynamic_fast();

    virtual void readMap();
    virtual void freeMap();

    double* onAxisField_m;
    gsl_spline *onAxisInterpolants_m[4];
    gsl_interp_accel *onAxisAccel_m[4];

    double hz_m;

    double frequency_m;
    double xlrep_m;

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double length_m;

    int num_gridpz_m;

    friend class Fieldmap;
};

#endif
