#ifndef CLASSIC_FIELDMAP1DPROFILE1_HH
#define CLASSIC_FIELDMAP1DPROFILE1_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DProfile1:public Fieldmap
{

 public:
  virtual bool getFieldstrength(const Vector_t &X, Vector_t &strength, Vector_t &info) const;
  virtual bool getFieldstrength_fdiff(const Vector_t &X, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
  virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
  virtual void swap();
  virtual void rescale(double factor);
  virtual void getInfo(Inform *);
  virtual double getFrequency() const;
  virtual void setFrequency(double freq);

 private:
  FM1DProfile1(string aFilename);
  ~FM1DProfile1();

  virtual void readMap();
  virtual void freeMap();

  string Filename_m;
  double *EngeCoefs_entry_m;
  double *EngeCoefs_exit_m;
  
  double scaleFactor_m;

  double zbegin_entry_m;
  double zend_entry_m;
  double polynomialOrigin_entry_m;
  int polynomialOrder_entry_m;

  double zbegin_exit_m;
  double zend_exit_m;
  double polynomialOrigin_exit_m;
  int polynomialOrder_exit_m;

  double length_m;
  double gapHeight_m;
  ofstream *testout_m;

  friend class Fieldmap;
};

#endif
