#ifndef CLASSIC_FIELDMAP1DDYNAMICFAST_HH
#define CLASSIC_FIELDMAP1DDYNAMICFAST_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DDynamic_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FM1DDynamic_fast(string aFilename);
    ~FM1DDynamic_fast();

    virtual void readMap();
    virtual void freeMap();

    double *FieldstrengthEz_m;    /**< 2D array with Ez, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
    double *FieldstrengthEr_m;    /**< 2D array with Er, read in like Ez*/
    double *FieldstrengthBt_m;    /**< 2D array with Er, read in like Ez*/

    double frequency_m;

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    double hr_m;                   /**< length between points in grid, r-direction, m*/
    int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    friend class Fieldmap;
};

#endif
