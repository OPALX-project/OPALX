#ifndef CSRWAKEFUNCTION_HH
#define CSRWAKEFUNCTION_HH

#include "Filters/Filter.h"
#include "Utilities/OpalFilter.h"
#include "Solvers/WakeFunction.hh"
#include <vector>

class ElementBase;

class CSRWakeFunction:public WakeFunction
{
public:
    CSRWakeFunction(const std::string &name, ElementBase* element, vector<Filter*> filters, const unsigned int &N);
  
    void apply(PartBunch &bunch);

    virtual const string getType() const;

private:
    double calcPsi(const double& x, const double& Ds) const;

    vector<Filter*> filters_m;
    LineDensity lineDensity_m;
    LineDensity dlineDensitydz_m;
    LineDensity d2lineDensitydz2_m;

    vector<double> Ez_m;
    double Begin_m;
    double Length_m;
    double R_m;
    unsigned int NBin_m;

};

#endif //CSRWAKEFUNCTION_HH
