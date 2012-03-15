#ifndef CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH
#define CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM2DMagnetoStatic:public Fieldmap
{

 public:
  bool getFieldstrength(Vector_t R, Vector_t &E, Vector_t &B);
  void readMap();
  double GetFrequency();
  void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
  void swap();
  void rescale(double factor);
  void getInfo(Inform *msg);
  double getFrequency() const;
  void setFrequency(double freq);

 private:
  FM2DMagnetoStatic(string aFilename);
  ~FM2DMagnetoStatic();

  string Filename_m;
  double *FieldstrengthBz_m;    /**< 2D array with Ez, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
  double *FieldstrengthBr_m;    /**< 2D array with Er, read in like Ez*/
  
  double scaleFactor_m;

  double rbegin_m;
  double rend_m;
  double zbegin_m;
  double zend_m;
  double hz_m;                   /**< length between points in grid, z-direction, m*/
  double hr_m;                   /**< length between points in grid, r-direction, m*/
  int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
  int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/
  
  bool swap_m;
  friend class Fieldmap;
};

#endif
