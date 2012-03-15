#ifndef CLASSIC_SeptumRep_HH
#define CLASSIC_SeptumRep_HH

// ------------------------------------------------------------------------
// $RCSfile: SeptumRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SeptumRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2009/09/21 10:21:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Septum.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

// Class SeptumRep
// ------------------------------------------------------------------------
/// Representation for Septum.

class SeptumRep: public Septum {

public:

  /// Constructor with given name.
  explicit SeptumRep(const string &name);

  SeptumRep();
  SeptumRep(const SeptumRep &);
  virtual ~SeptumRep();

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
  //  Return the element geometry.
  //  Version for non-constant object.
  /// Get geometry.
  //  Return the element geometry.
  //  Version for non-constant object.
   virtual StraightGeometry &getGeometry();

  /// Get geometry.
  //  Return the element geometry
  //  Version for constant object.
  virtual const StraightGeometry &getGeometry() const;

  /// Construct an image.
  //  Return the image of the element, containing the name and type string 
  //  of the element, and a copy of the user-defined attributes.
  virtual ElementImage *getImage() const;

  /// Set active flag.
  virtual void setActive(bool = true);

protected:

  /// The zero magnetic field.
  NullField field;

  /// The septa's geometry.
  StraightGeometry geometry;

  /// The active/inactive flag.
  bool active;

private:

  // Not implemented.
  void operator=(const SeptumRep &);


 
};

#endif // CLASSIC_SeptumRep_HH
