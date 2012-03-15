#ifndef CLASSIC_MonitorRep_HH
#define CLASSIC_MonitorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: MonitorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MonitorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Monitor.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

// Class MonitorRep
// ------------------------------------------------------------------------
/// Representation for an orbit position monitor.
//  The base class observes both planes.

class MonitorRep: public Monitor {

public:

  /// Constructor with given name.
  explicit MonitorRep(const string &name);

  MonitorRep();
  MonitorRep(const MonitorRep &);
  virtual ~MonitorRep();

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
  virtual StraightGeometry &getGeometry();

  /// Get geometry.
  //  Return the element geometry.
  //  Version for constant object.
  virtual const StraightGeometry &getGeometry() const;

  /// Construct an image.
  //  Return the image of the element, containing the name and type string 
  //  of the element, and a copy of the user-defined attributes.
  virtual ElementImage *getImage() const;

  /// Get planes.
  //  Return the plane(s) observed by this monitor.
  virtual Plane getPlane() const;

  /// Get element type string.
  virtual const string &getType() const;

  /// Set active flag.
  //  If [b]flag[/b] is true, the monitor is activated, otherwise it is
  //  deactivated.
  virtual void setActive(bool = true);

protected:

  /// The zero magnetic field.
  NullField field;

  /// The monitor geometry.
  StraightGeometry geometry;

  /// The active/inactive flag.
  bool active;

private:

  // Not implemented.
  void operator=(const MonitorRep &);
};

#endif // CLASSIC_MonitorRep_HH
