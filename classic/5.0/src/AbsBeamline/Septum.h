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

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);
  
  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;

  virtual string getType() {return "Septum";}

private:

  // Not implemented.
  void operator=(const Septum &);
};

#endif // CLASSIC_Septum_HH
