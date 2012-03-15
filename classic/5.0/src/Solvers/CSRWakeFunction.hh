#ifndef CSRWAKEFUNCTION_HH
#define CSRWAKEFUNCTION_HH

#include "Filters/Filter.h"
#include "Utilities/OpalFilter.h"
#include "Solvers/WakeFunction.hh"
#include <vector>

class ElementBase;

class CSRWakeFunction: public WakeFunction {
public:
    CSRWakeFunction(const std::string &name, ElementBase *element, std::vector<Filter *> filters, const unsigned int &N);

    void apply(PartBunch &bunch);

    virtual const string getType() const;

private:
    double calcPsi(const double &psiInitial, const double &x, const double &Ds) const;

    std::vector<Filter *> filters_m;
    LineDensity lineDensity_m;
    LineDensity dlineDensitydz_m;
    LineDensity d2lineDensitydz2_m;

    // Longitudinal CSR field.
    std::vector<double> Ez_m;

    // Retarded angle Psi;
    std::vector<double> Psi_m;

    // Start position of CSR wake.
    double Begin_m;

    // Start position of equivalent hard edge dipole that approximates actual
    // dipole.
    double FieldBegin_m;

    // Effective length of dipole.
    double Length_m;

    // Radius of curvature of effective dipole.
    double R_m;

    // Bend angle.
    double Phi_m;

    unsigned int NBin_m;

};

#endif //CSRWAKEFUNCTION_HH
