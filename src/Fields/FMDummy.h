#ifndef CLASSIC_FIELDMAPDUMMY
#define CLASSIC_FIELDMAPDUMMY

#include "Fields/Fieldmap.h"

class FMDummy: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t<double, 3> &R, Vector_t<double, 3> &E, Vector_t<double, 3> &B) const;
    virtual bool getFieldDerivative(const Vector_t<double, 3> &R, Vector_t<double, 3> &E, Vector_t<double, 3> &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    FMDummy(std::string aFilename);
    ~FMDummy();

    virtual void readMap();
    virtual void freeMap();

    double zbegin_m;
    double zend_m;

    friend class Fieldmap;
};

#endif
