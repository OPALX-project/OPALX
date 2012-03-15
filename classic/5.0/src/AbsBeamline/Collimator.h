#ifndef CLASSIC_Collimator_HH
#define CLASSIC_Collimator_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Collimator
//   Defines the abstract interface for a beam Collimator.
//   *** MISSING *** Collimator interface is still incomplete.
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


// Class Collimator
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class Collimator defines the abstract interface for a collimator.

class Collimator: public Component {

public:

  /// Constructor with given name.
  explicit Collimator(const string &name);

  Collimator();
  Collimator(const Collimator &rhs);
  virtual ~Collimator();

  /// Apply visitor to Collimator.
  virtual void accept(BeamlineVisitor &) const;

  /// Return the horizontal half-aperture.
  virtual double getXsize() const = 0;

  /// Return the vertical half-aperture.
  virtual double getYsize() const = 0;

private:

  // Not implemented.
  void operator=(const Collimator &);
};

#endif // CLASSIC_Collimator_HH
