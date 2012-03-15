#ifndef CLASSIC_Monitor_HH
#define CLASSIC_Monitor_HH

// ------------------------------------------------------------------------
// $RCSfile: Monitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/StraightGeometry.h"


// Class Monitor
// ------------------------------------------------------------------------
/// Interface for beam position monitors.
//  Class Monitor defines the abstract interface for general beam position
//  monitors.

class Monitor: public Component {

public:

  /// Plane selection.
  enum Plane {
    /// Monitor is off (inactive).
    OFF,
    /// Monitor acts on x-plane.
    X,
    /// Monitor acts on y-plane.
    Y,
    /// Monitor acts on both planes.
    XY
  };

  /// Constructor with given name.
  explicit Monitor(const string &name);

  Monitor();
  Monitor(const Monitor &);
  virtual ~Monitor();

  /// Apply visitor to Monitor.
  virtual void accept(BeamlineVisitor &) const;

  /// Get geometry.
  virtual StraightGeometry &getGeometry() = 0;

  /// Get geometry. Version for const object.
  virtual const StraightGeometry &getGeometry() const = 0;

  /// Get plane on which monitor observes.
  virtual Plane getPlane() const = 0;

private:

  // Not implemented.
  void operator=(const Monitor &);
};

#endif // CLASSIC_Monitor_HH
