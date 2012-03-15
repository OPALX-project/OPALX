#ifndef CLASSIC_RFQuadrupole_HH
#define CLASSIC_RFQuadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: RFQuadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Class: RFQuadrupole
//   *** MISSING *** RFQuadrupole interface is incomplete.
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class RFQuadrupole
// ------------------------------------------------------------------------
/// Interface for RF Quadrupole.
//  Class RFQuadrupole defines the abstract interface for a RF Quadrupole.

class RFQuadrupole: public Component {

public:

  /// Constructor with given name.
  explicit RFQuadrupole(const string &name);

  RFQuadrupole();
  RFQuadrupole(const RFQuadrupole &);
  virtual ~RFQuadrupole();

  /// Apply visitor to RFQuadrupole.
  virtual void accept(BeamlineVisitor &) const;

private:

  // Not implemented.
  void operator=(const RFQuadrupole &);
};

#endif // CLASSIC_RFQuadrupole_HH
