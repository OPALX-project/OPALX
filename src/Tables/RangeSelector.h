#ifndef OPAL_RangeSelector_HH
#define OPAL_RangeSelector_HH

// ------------------------------------------------------------------------
// $RCSfile: RangeSelector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeSelector
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "AbstractObjects/RangeRep.h"

class Beamline;


// Class RangeSelector
// ------------------------------------------------------------------------
/// An abstract visitor which calls the pure virtual method
//  [tt]RangeSelector::handleXXX()[/tt] for each element or beamline in range.

class RangeSelector: public DefaultVisitor {

public:

  /// Constructor.
  //  Attach visitor to a beamline, remember the range.
  RangeSelector(const Beamline &, const RangeRep &range);

  virtual ~RangeSelector();

  /// Execute the algorithm.
  virtual void execute();

  /// Apply the visitor to an FlaggedElmPtr.
  virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

protected:

  /// The operation to be done for beamlines.
  //  When overriding, make sure the beamline members are handled.
  virtual void handleBeamline(const FlaggedElmPtr &);

  /// The operation to be done for elements.
  virtual void handleElement(const FlaggedElmPtr &);

  /// Working data for range.
  RangeRep itsRange;

private:

  // Not implemented.
  RangeSelector();
  RangeSelector(const RangeSelector &);
  void operator=(const RangeSelector &);
};

#endif // OPAL_RangeSelector_HH
