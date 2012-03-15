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

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);
  
  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;
  
  virtual string getType() {return "RFQuadrupole";}

  virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

  // Not implemented.
  void operator=(const RFQuadrupole &);
};

#endif // CLASSIC_RFQuadrupole_HH
