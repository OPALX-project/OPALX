#ifndef CLASSIC_PatchRep_HH
#define CLASSIC_PatchRep_HH

// ------------------------------------------------------------------------
// $RCSfile: PatchRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PatchRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Patch.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Fields/NullField.h"


// Class PatchRep
// ------------------------------------------------------------------------
/// Representation for a geometry patch.

class PatchRep: public Patch {

public:

  /// Constructor with given name.
  explicit PatchRep(const string &name);

  PatchRep();
  PatchRep(const PatchRep &);
  virtual ~PatchRep();

  /// Return clone.
  //  Return an identical deep copy of the element.
  virtual ElementBase *clone() const;

  /// Construct a read/write channel.
  //  This method constructs a Channel permitting read/write access to
  //  the attribute [b]aKey[/b] and returns it.
  //  If the attribute does not exist, it returns NULL.
  virtual Channel *getChannel(const string &aKey);

  /// Get field.
  //  Version for non-constant object.
  virtual NullField &getField();

  /// Get field.
  //  Version for constant object.
  virtual const NullField &getField() const;

  /// Get geometry.
  //  Version for non-constant object.
  virtual NullGeometry &getGeometry();

  /// Get geometry.
  //  Version for constant object.
  virtual const NullGeometry &getGeometry() const;

  /// Get patch.
  //  Return the geometric patch to transform from input to output.
  virtual const Euclid3D &getPatch() const;

  /// Set patch.
  //  Assign the geometric patch to transform from input to output.
  void setPatch(const Euclid3D &euclid);

  /// Set patch.
  //  Assign the geometric patch to transform from input to output,
  //  using the displacement (x,y,z) and the rotation vector (vx,vy,vz).
  void setPatch(double x, double y, double z,
		double vx, double vy, double vz);

  /// Construct an image.
  //  Return the image of the element, containing the name and type string 
  //  of the element, and a copy of the user-defined attributes.
  virtual ElementImage *getImage() const;

  /// Get displacement.
  //  Return the x-displacement.
  double getX() const;

  /// Get displacement.
  //  Return the y-displacement.
  double getY() const;

  /// Get displacement.
  //  Return the z-displacement.
  double getZ() const;

  /// Get rotation.
  //  Return the rotation around the x-axis
  double getVX() const;

  /// Get rotation.
  //  Return the rotation around the y-axis
  double getVY() const;

  /// Get rotation.
  //  Return the rotation around the z-axis
  double getVZ() const;

  /// Set displacement.
  //  Assign the x-displacement.
  void setX(double);

  /// Set displacement.
  //  Assign the y-displacement.
  void setY(double);

  /// Set displacement.
  //  Assign the z-displacement.
  void setZ(double);

  /// Set rotation.
  //  Assign the rotation around the x-axis
  void setVX(double);

  /// Set rotation.
  //  Assign the rotation around the y-axis
  void setVY(double);

  /// Set rotation.
  //  Assign the rotation around the z-axis
  void setVZ(double);

private:

  // Not implemented.
  void operator=(const PatchRep &);

  /// The zero magnetic field.
  NullField field;

  /// The geometry.
  NullGeometry geometry;

  /// The geometry patch.
  Euclid3D patch;
};

#endif // CLASSIC_PatchRep_HH
