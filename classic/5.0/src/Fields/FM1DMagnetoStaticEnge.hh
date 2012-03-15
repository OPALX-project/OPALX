#ifndef CLASSIC_FIELDMAP1DMAGNETOSTATICENGE_HH
#define CLASSIC_FIELDMAP1DMAGNETOSTATICENGE_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DMagnetoStaticEnge:public Fieldmap
{

 public:
  virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
  virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
  virtual void swap();
  virtual void rescale(double factor);
  virtual void getInfo(Inform *);
  virtual double getFrequency() const;
  virtual void setFrequency(double freq);

 private:
  FM1DMagnetoStaticEnge(string aFilename);
  ~FM1DMagnetoStaticEnge();

  virtual void readMap();
  virtual void freeMap();

  string Filename_m;
  double *EngeCoefs_entry_m;
  double *EngeCoefs_exit_m;
  
  double scaleFactor_m;

  double zbegin_m;
  double zend_m;
  double length_m;
  double field_sections_m[6];
  // field_sections_m[0] == start of entry fringe field (where field > 0)
  // field_sections_m[1] == origin of the polynomial of the entry fringe field
  // field_sections_m[2] == end of the entry fringe field (where normalised field == 1)
  // field_sections_m[3] == start of exit fringe field (where normalised field < 1)
  // field_sections_m[4] == origin of the polynomial of the exit fringe field
  // field_sections_m[5] == end of the exit fringe field (where field == 0)
  double fringe_length_entry_m;
  double fringe_length_exit_m;

  int polynomialOrder_m;

  friend class Fieldmap;
};

namespace QRDecomposition
{
  void solve(double* Matrix, double* Solution, double* rightHandSide, const int &M, const int &N);
}

#endif
