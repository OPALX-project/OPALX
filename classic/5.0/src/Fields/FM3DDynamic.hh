#ifndef CLASSIC_FIELDMAP3DDYNAMIC_HH
#define CLASSIC_FIELDMAP3DDYNAMIC_HH

#include "Fields/Fieldmap.hh"

using namespace std;

class FM3DDynamic:public Fieldmap
{

 public:
  virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
  virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
  virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
  virtual void swap();
  virtual void getInfo(Inform *msg);
  virtual double getFrequency() const;
  virtual void setFrequency(double freq);

 private:
  FM3DDynamic(string aFilename);
  ~FM3DDynamic();

  virtual void readMap();
  virtual void freeMap();

  double *FieldstrengthEz_m;    /**< 3D array with Ez */
  double *FieldstrengthEx_m;    /**< 3D array with Ex */
  double *FieldstrengthEy_m;    /**< 3D array with Ey */
  double *FieldstrengthHz_m;    /**< 3D array with Hz */
  double *FieldstrengthHx_m;    /**< 3D array with Hx */
  double *FieldstrengthHy_m;    /**< 3D array with Hy */

  double frequency_m;

  double xbegin_m;
  double xend_m;

  double ybegin_m;
  double yend_m;

  double zbegin_m;
  double zend_m;

  double hx_m;                   /**< length between points in grid, x-direction */
  double hy_m;                   /**< length between points in grid, y-direction */
  double hz_m;                   /**< length between points in grid, z-direction */
  int num_gridpx_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
  int num_gridpy_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
  int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/
  
  bool swap_m;
  friend class Fieldmap;
};

#endif
