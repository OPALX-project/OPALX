#ifndef CLASSIC_FIELDMAP1DDYNAMIC_HH
#define CLASSIC_FIELDMAP1DDYNAMIC_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM1DDynamic:public Fieldmap
{

 public:
  bool getFieldstrength(Vector_t R, Vector_t &E, Vector_t &B);
  void readMap();
  void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
  void swap();
  void rescale(double factor);
  void getInfo(Inform *);
  double getFrequency() const;
  void setFrequency(double freq);

 private:
  FM1DDynamic(string aFilename);
  ~FM1DDynamic();

  string Filename_m;
  double *realFourCoefs_m;
  double *imagFourCoefs_m;
  
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
