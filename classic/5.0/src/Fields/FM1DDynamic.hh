#ifndef CLASSIC_FIELDMAP1DDYNAMIC_HH
#define CLASSIC_FIELDMAP1DDYNAMIC_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DDynamic:public Fieldmap
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
  FM1DDynamic(string aFilename);
  ~FM1DDynamic();

  virtual void readMap();
  virtual void freeMap();

  string Filename_m;
  double* restrict realFourCoefs_m;
  double* restrict imagFourCoefs_m;
  
  double frequency_m;
  double xlrep_m;

  double scaleFactor_m;

  double rbegin_m;
  double rend_m;
  double zbegin_m;
  double zend_m;
  double length_m;
  
  int accuracy_m;

  friend class Fieldmap;
};

#endif
