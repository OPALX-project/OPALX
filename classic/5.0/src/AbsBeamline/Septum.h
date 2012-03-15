#ifndef CLASSIC_Septum_HH
#define CLASSIC_Septum_HH

// ------------------------------------------------------------------------
// $RCSfile: Septum.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Septum
//   *** MISSING *** Septum interface is still incomplete.
// 
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Septum
// ------------------------------------------------------------------------
/// Interface for septum magnet.
//  Class Septum defines the abstract interface for a septum magnet.

class Septum: public Component {

public:

  /// Constructor with given name.
  explicit Septum(const string &name);

  Septum();
  Septum(const Septum &);
  virtual ~Septum();

  /// Apply visitor to Septum.
  virtual void accept(BeamlineVisitor &) const;

private:

  // Not implemented.
  void operator=(const Septum &);
};

#endif // CLASSIC_Septum_HH
