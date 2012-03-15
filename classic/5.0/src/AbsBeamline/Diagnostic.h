#ifndef CLASSIC_Diagnostic_HH
#define CLASSIC_Diagnostic_HH

// ------------------------------------------------------------------------
// $RCSfile: Diagnostic.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Diagnostic
//   *** MISSING *** Diagnostic interface is still incomplete.
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


// Class Diagnostic
// ------------------------------------------------------------------------
/// Interface for beam diagnostics.
//  Class Diagnostic defines the abstract interface for a beam diagnostic.

class Diagnostic: public Component {

public:

  /// Constructor with given name.
  explicit Diagnostic(const string &name);

  Diagnostic();
  Diagnostic(const Diagnostic &rhs);
  virtual ~Diagnostic();

  /// Apply visitor to Diagnostic.
  virtual void accept(BeamlineVisitor &) const;

private:

  // Not implemented.
  void operator=(const Diagnostic &);
};

#endif // CLASSIC_Diagnostic_HH
