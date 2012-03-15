#ifndef OPAL_ThreadBpm_HH
#define OPAL_ThreadBpm_HH

// ------------------------------------------------------------------------
// $RCSfile: ThreadBpm.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThreadBpm
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/CorrectionBase.h"
#include "Algorithms/OrbitTracker.h"
#include "MemoryManagement/OwnPtr.h"

/// Class ThreadBpm
// ------------------------------------------------------------------------
/// OPAL command THREADBPM.

class ThreadBpm: public CorrectionBase {

public:

  /// Exemplar constructor.
  ThreadBpm();

  virtual ~ThreadBpm();

  /// Make clone.
  virtual ThreadBpm *clone(const string &name);

  /// Check validity of the table definition.
  virtual void execute();

private:

  /// The additional attributes for the "THREADBPM" command.
  enum {
    TOL = CorrectionBase::SIZE, // Tolerance for the closed orbit positions.
    LISTC,      // List the correctors after correction.
    LISTM,      // List the monitors after correction.
    SIZE
  };

  // Not implemented.
  ThreadBpm(const ThreadBpm &);
  void operator=(const ThreadBpm &);

  /// Clone constructor.
  ThreadBpm(const string &name, ThreadBpm *parent);

  // Backtrack to correct in one plane.
  void correct(int plane, TLine::iterator &);

  /// Find the closed orbit by threading.
  bool thread(double tol);

  // The tracker.
  OwnPtr<OrbitTracker> itsTracker;
};

#endif // OPAL_ThreadBpm_HH
