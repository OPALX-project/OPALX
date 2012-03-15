#ifndef CLASSIC_FIELDMAP1DELECTROSTATICFAST_HH
#define CLASSIC_FIELDMAP1DELECTROSTATICFAST_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DElectroStatic_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FM1DElectroStatic_fast(string aFilename);
    ~FM1DElectroStatic_fast();

    virtual void readMap();
    virtual void freeMap();

    double *FieldstrengthEz_m;
    double *FieldstrengthEr_m;

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
