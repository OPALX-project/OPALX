#ifndef OPAL_TrackRun_HH
#define OPAL_TrackRun_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackRun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackRun
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:12 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Algorithms/Tracker.h"
#include "Structure/FieldSolver.h"
#include "Structure/DataSink.h"
#include "Algorithms/bet/SLDataSink.h"
#include "Distribution/Distribution.h"

// Class TrackRun
// ------------------------------------------------------------------------
/// The RUN command.

class TrackRun: public Action {

public:

  /// Exemplar constructor.
  TrackRun();

  virtual ~TrackRun();

  /// Make clone.
  virtual TrackRun *clone(const string &name);
  
  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  TrackRun(const TrackRun &);
  void operator=(const TrackRun &);

  // Clone constructor.
  TrackRun(const string &name, TrackRun *parent);

  // Pointer to tracking algorithm.
  Tracker *itsTracker;

  Distribution *dist;

  FieldSolver  *fs;
  
  DataSink *ds;

#ifdef HAVE_ENVELOPE_SOLVER
  SLDataSink *slds;
#endif


};

#endif // OPAL_TrackRun_HH
