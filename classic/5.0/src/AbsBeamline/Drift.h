#ifndef CLASSIC_Drift_HH
#define CLASSIC_Drift_HH

// ------------------------------------------------------------------------
// $RCSfile: Drift.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
//   Defines the abstract interface for a drift space.
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


// Class Drift
// ------------------------------------------------------------------------
/// Interface for drift space.
//  Class Drift defines the abstract interface for a drift space.

class Drift: public Component {

public:

  /// Constructor with given name.
  explicit Drift(const string &name);

  Drift();
  Drift(const Drift &right);
  virtual ~Drift();

  /// Apply visitor to Drift.
  virtual void accept(BeamlineVisitor &) const;

private:

  // Not implemented.
  void operator=(const Drift &);
};

#endif // CLASSIC_Drift_HH
