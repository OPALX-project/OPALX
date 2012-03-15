#ifndef CLASSIC_RBend_HH
#define CLASSIC_RBend_HH

// ------------------------------------------------------------------------
// $RCSfile: RBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Structure/Wake.h"
#include <vector>


// Class RBend
// ------------------------------------------------------------------------
/// Interface for rectangular bend.
//  Class RBend defines the abstract interface for rectangular bend magnets.
//  A rectangular bend magnet has a rectilinear geometry about which its
//  multipole components are specified.

class RBend: public Component {

public:

  /// Constructor with given name.
  explicit RBend(const string &name);

  RBend();
  RBend(const RBend &);
  virtual ~RBend();

  /// Apply visitor to RBend.
  virtual void accept(BeamlineVisitor &) const;

  /// Get dipole field of RBend.
  virtual double getB() const = 0;

  /// Get RBend geometry.
  //  Version for non-constant object.
  virtual RBendGeometry &getGeometry() = 0;

  /// Get RBend geometry
  //  Version for constant object.
  virtual const RBendGeometry &getGeometry() const = 0;

  /// Get multipole expansion of field.
  //  Version for non-constant object.
  virtual BMultipoleField &getField() = 0;

  /// Get multipole expansion of field.
  //  Version for constant object.
  virtual const BMultipoleField &getField() const = 0;

  /// Get normal component.
  //  Return the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getNormalComponent(int) const;

  /// Get skew component.
  //  Return the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the return value is zero.
  double getSkewComponent(int) const;

  /// Set normal component.
  //  Set the normal component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setNormalComponent(int, double);

  /// Set skew component.
  //  Set the skew component of order [b]n[/b] in T/m**(n-1).
  //  If [b]n[/b] is larger than the maximum order, the component is created.
  void setSkewComponent(int, double);

  /// Get pole entry face rotation.
  //  Return the rotation of the entry pole face with respect to the x-axis.
  //  A positive angle rotates the pole face normal away from the centre
  //  of the machine.
  virtual double getEntryFaceRotation() const = 0;

  /// Get exit pole face rotation.
  //  Return the rotation of the exit pole face with respect to the x-axis.
  //  A positive angle rotates the pole face normal away from the centre
  //  of the machine.
  virtual double getExitFaceRotation() const = 0;

  /// Get entry pole face curvature.
  //  Return the curvature of the entry pole face.
  //  A positive curvature creates a convex pole face.
  virtual double getEntryFaceCurvature() const = 0;
  
  /// Get exit pole face curvature.
  //  Return the curvature of the exit pole face.
  //  A positive curvature creates a convex pole face.
  virtual double getExitFaceCurvature() const = 0;

  /// Get number of slices.
  //  Slices and stepsize used to determine integration step.
  virtual double getSlices() const = 0;

  /// Get stepsize.
  //  Slices and stepsize used to determine integration step.
  virtual double getStepsize() const = 0;

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);
  
  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;

  void setAmplitudem(double vPeak);

  void setGapWidth(double gapwidth);

  void setFaceAngleEntry(double angle);

  void setFaceAngleExit(double angle);

  /// Set the name of the "field map"; 
  /// For a bend this file contains the coefficients for the Enge function 
  void setFieldMapFN(string fmapfn);

  string getFieldMapFN() const;

  void setEngeCoefs(const vector<double> EngeCoefs);

  virtual string getType() { return "RBend";}

  void setWakeFunction(Wake *wfstr);

private:

  // Not implemented.
  void operator=(const RBend &);
  string filename_m;             /**< The name of the inputfile*/
  double amplitude_m;
  double cos_B_xy_angle_m;
  double sin_B_xy_angle_m;
  double s_entry_m;
  double s_exit_m;
  double sin_face_angle_entry_m;
  double sin_face_angle_exit_m;
  double cos_face_angle_entry_m;
  double cos_face_angle_exit_m;
  double tan_face_angle_entry_m;
  double tan_face_angle_exit_m;
  double fringe_length_entry_m;
  double fringe_length_exit_m;
  double gap_width_m;
  double* enge_coefficients_m;
  unsigned int enge_polynomial_order_m;
  double enge_shift_entry_m;
  double enge_shift_exit_m;

  /// the name of the wake field to be used
  Wake  *wf_func_m;
};

#endif // CLASSIC_RBend_HH
