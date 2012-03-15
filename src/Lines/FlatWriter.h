#ifndef OPAL_FlatWriter_HH
#define OPAL_FlatWriter_HH 1

// ------------------------------------------------------------------------
// $RCSfile: FlatWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FlatWriter
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "Lines/Sequence.h"
#include <iosfwd>

class Beamline;
class Drift;
class ElementBase;
class FlaggedElmPtr;


// Class FlatWriter
// ------------------------------------------------------------------------
/// Write a line or sequence as a flat sequence in OPAL-8 format.

class FlatWriter: public DefaultVisitor {

public:

  FlatWriter(const Beamline &, Sequence::TLine &seq, std::ostream &os);
  ~FlatWriter();

  /// Return the temporary sequence.
  const Sequence::TLine &getSequence() const;
  
  /// Override drift operation.
  //  Accumulate length, but do not save anything.
  virtual void visitDrift(const Drift &);
  
  /// Override beamline exit.
  virtual void visitMapIntegrator(const MapIntegrator &);

protected:

  /// Apply default.
  //  Add the element to the sequence list.
  virtual void applyDefault(const ElementBase &);
  
private:
  
  // The current longitudinal position.
  double itsPosition;

  // The temporary flat sequence.
  Sequence::TLine &itsSequence;

  // The output stream.
  std::ostream &itsStream;
};

#endif // OPAL_FlatWriter_HH
