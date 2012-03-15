#ifndef CLASSIC_Solenoid_HH
#define CLASSIC_Solenoid_HH

// ------------------------------------------------------------------------
// $RCSfile: Solenoid.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Fields/Fieldmap.hh"

// Class Solenoid
// ------------------------------------------------------------------------
/// Interface for solenoids.
//  Class Solenoid defines the abstract interface for solenoid magnets.


class Solenoid: public Component {

public:

  /// Constructor with given name.
  explicit Solenoid(const string &name);

  Solenoid();
  Solenoid(const Solenoid &);
  virtual ~Solenoid();

  /// Apply visitor to Solenoid.
  virtual void accept(BeamlineVisitor &) const;

  /// Get solenoid field Bz in Teslas.
  virtual double getBz() const = 0;

  void setKS(double ks);

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply( const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);

  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;

  virtual void goOnline();

  virtual void goOffline();

  //  Assign the field filename.
  void setFieldMapFN(string fn);

  void setFast(bool fast);

  bool getFast() const;

  void setMisalignment(double x, double y, double z);

  void getMisalignment(double &x, double &y, double &z) const;

  virtual string getType() { return "Solenoid";}

private:

  //  string name;                   /**< The name of the object*/
  string filename_m;               /**< The name of the inputfile*/
  Fieldmap *myFieldmap;
  double scale_m;                /**< scale multiplier*/
  double startField_m;           /**< startingpoint of field, m*/
  double lengthUnit_m;

  double dx_m;
  double dy_m;
  double dz_m;

  bool fast_m;
  // Not implemented.
  void operator=(const Solenoid &);
};

#endif // CLASSIC_Solenoid_HH
