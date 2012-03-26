#ifndef CLASSIC_FIELDMAP1DDYNAMIC_HH
#define CLASSIC_FIELDMAP1DDYNAMIC_HH

#include "Fields/Fieldmap.hh"

class FM1DDynamic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

private:
    FM1DDynamic(std::string aFilename);
    ~FM1DDynamic();

    virtual void readMap();
    virtual void freeMap();

    double *restrict FourCoefs_m;

    double frequency_m;
    double xlrep_m;

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double length_m;
    double trafo_phase_m;

    int accuracy_m;
    int num_gridpz_m;

    friend class Fieldmap;
};

#endif
