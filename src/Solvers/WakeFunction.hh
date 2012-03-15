#ifndef WAKEFUNCTION_HH
#define WAKEFUNCTION_HH

#include "Algorithms/PartData.h"
#include "Algorithms/PartBunch.h"

class WakeFunction
{
 public:
  WakeFunction(PartData &ref);

  virtual void apply(PartBunch &bunch) = 0;

  virtual const string getType() const = 0;

 private:
  PartData &reference_m;
};

inline WakeFunction::WakeFunction(PartData &ref):
  reference_m(ref)
{ }

class LineDensity: public vector<double>
{
public:
  LineDensity(int size = 0, double defaultValue = 0.0):
    vector<double>(size,defaultValue)
  { }
  
  void getFirstDerivative(vector<double> &firstDerivative, const double &hz);
};

inline void LineDensity::getFirstDerivative(vector<double> &firstDerivative, const double &hz)
{
  const int size = this->size();
  if (firstDerivative.size() != size)
    firstDerivative.resize(size,0.0);

  firstDerivative[0] = ((*this)[1] - (*this)[0])/hz;
  for (int i = 1; i < size - 1; ++i)
    firstDerivative[i] = ((*this)[i + 1] - (*this)[i - 1])/hz;
  firstDerivative[size - 1] = ((*this)[size - 1] - (*this)[size - 2])/hz;
}

#endif // WAKEFUNCTION_HH
