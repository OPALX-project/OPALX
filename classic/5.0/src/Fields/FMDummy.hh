#ifndef CLASSIC_FIELDMAPDUMMY
#define CLASSIC_FIELDMAPDUMMY

#include "Fields/Fieldmap.hh"

class FMDummy: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FMDummy(string aFilename);
    ~FMDummy();

    virtual void readMap();
    virtual void freeMap();

    double zbegin_m;
    double zend_m;

    friend class Fieldmap;
};

#endif
