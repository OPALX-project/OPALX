#ifndef CLASSIC_Lambertson_HH
#define CLASSIC_Lambertson_HH

// ------------------------------------------------------------------------
// $RCSfile: Lambertson.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   *** MISSING *** Lambertson interface is still incomplete.
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


// Class Lambertson
// ------------------------------------------------------------------------
/// Interface for a Lambertson septum.
//  Class Lambertson defines the abstract interface for a Lambertson
//  septum magnet.

class Lambertson: public Component {

public:

  /// Constructor with given name.
  explicit Lambertson(const string &name);

  Lambertson();
  Lambertson(const Lambertson &);
  virtual ~Lambertson();

  /// Apply visitor to Lambertson.
  virtual void accept(BeamlineVisitor &) const;

private:

  // Not implemented.
  void operator=(const Lambertson &);
};

#endif // CLASSIC_Lambertson_HH
