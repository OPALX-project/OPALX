#ifndef CSRWAKEFUNCTION_HH
#define CSRWAKEFUNCTION_HH

#include "Algorithms/Filters.h"
#include "Solvers/WakeFunction.hh"
#include <vector>

class ElementBase;

class CSRWakeFunction:public WakeFunction
{
 public:
  CSRWakeFunction(PartData &reference, ElementBase* element, const S_G_FilterOptions &options1 = S_G_FilterOptions(33,16,16,4), 
                  const S_G_FilterOptions &options2 = S_G_FilterOptions(33,16,16,4));
  
  void apply(PartBunch &bunch);

  virtual const string getType() const;

 private:
  double calcPsi(double x, double Ds);

  SavitzkyGolayFilter smoother1_m;
  SavitzkyGolayFilter smoother2_m;
  LineDensity lineDensity_m;
  LineDensity dlineDensitydz_m;
  LineDensity d2lineDensitydz2_m;

  vector<double> Ez_m;
  double Begin_m;
  double Length_m;
  double R_m;
};

#endif //CSRWAKEFUNCTION_HH
