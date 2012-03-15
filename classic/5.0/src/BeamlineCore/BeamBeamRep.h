#ifndef CLASSIC_BeamBeamRep_HH
#define CLASSIC_BeamBeamRep_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamBeamRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeamRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamBeam.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Fields/NullField.h"


// Class BeamBeamRep
// ------------------------------------------------------------------------
/// Representation of beam-beam interaction.
//  A concrete representation for a beam-beam interaction.
//  The opposite bunch has a Gaussian distribution of the form
//  [P]
//  N(x,y,s) = c * exp(- transpose(z - delta) sigma**(-1) (z - delta)),
//  [P]
//  with the definitions
//  [DL]
//  [DT]z[DD]     position vector in space,
//  [DT]delta[DD] centroid of opposite bunch in global reference system,
//  [DT]sigma[DD] beam matrix in the TRANSPORT sense,
//  [DT]c[DD]     a normalising factor, such that the total charge is Q,
//  [DT]sigma(i,i)[DD] standard deviation sigma(i) in direction i,
//  [DT]r(i,j)[DD] = sigma(i,j)/sqrt(sigma(i,i)*sigma(i,j)),
//    correlations between phase space coordinates i and j.
//  [/DL]

class BeamBeamRep: public BeamBeam {

public:

  /// Constructor with given name.
  explicit BeamBeamRep(const string &name);

  BeamBeamRep();
  BeamBeamRep(const BeamBeamRep &);
  virtual ~BeamBeamRep();

  /// Return clone.
  //  Return an identical deep copy of the element.
  virtual ElementBase *clone() const;

  /// Construct a read/write channel.
  //  This method constructs a Channel permitting read/write access to
  //  the attribute[b]aKey[/b] and returns it.
  //  If the attribute does not exist, it returns NULL.
  virtual Channel *getChannel(const string &aKey);

  /// Get bunch charge.
  //  Return the number of particles times the particle charge in the 
  //  opposite bunch.  Units are proton charges.
  virtual double getBunchCharge() const;

  /// Get moments.
  //  Return the moment matrix for the opposite bunch (the matrix of
  //  second momenta).  Units are square metres.
  virtual const Matrix3D &getBunchMoment() const;

  /// Get displacement.
  //  Return the displacement vector for position of opposite bunch.
  //  Units are metres.
  virtual const Vector3D &getBunchDisplacement() const;

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

  /// Construct an image.
  //  Return the image of the element, containing the name and type string 
  //  of the element, and a copy of the user-defined attributes.
  virtual ElementImage *getImage() const;

  /// Return type name string
  virtual const string &getType() const;

  /// Set the bunch charge.
  //  Units are proton charges.
  void setBunchCharge(double);

  /// Set moments.
  //  Assign the moment matrix for the opposite bunch
  //  (the matrix of second momenta).  Units are square metres.
  void setBunchMoment(const Matrix3D &);

  /// Set displacement.
  //  Assign the displacement vector for the opposite bunch.
  //  Units are metres.
  void setBunchDisplacement(const Vector3D &);

private:

  // Not implemented.
  void operator=(const BeamBeamRep &);

  // Get a moment for the opposite bunch.
  double getMoment(int index) const;

  // Get the displacement vector for the opposite bunch.
  double getDisplacement(int index) const;

  // Set a moment for the opposite bunch.
  void setMoment(int index, double);

  // Set a displacement for the opposite bunch.
  void setDisplacement(int index, double);

  // The zero magnetic field.
  NullField field;

  // The geometry (a straight geometry with zero length).
  NullGeometry geometry;

  // Total charge of the opposite bunch in proton charges.
  double Q;

  // Displacement of the centre of gravity of the opposite bunch
  // with respect to the local frame for the affected bunch.
  Vector3D delta;

  // Second moments of the opposite bunch.
  Matrix3D sigma;
};

#endif // CLASSIC_BeamBeamRep_HH
