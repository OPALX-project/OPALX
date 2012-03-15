#ifndef OPAL_ThreadAll_HH
#define OPAL_ThreadAll_HH

// ------------------------------------------------------------------------
// $RCSfile: ThreadAll.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThreadAll
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


/// Class ThreadAll
// ------------------------------------------------------------------------
/// OPAL command THREADALL.

class ThreadAll: public CorrectionBase {

public:

  /// Exemplar constructor.
  ThreadAll();

  virtual ~ThreadAll();

  /// Make clone.
  virtual ThreadAll *clone(const string &name);

  /// Check validity of the table definition.
  virtual void execute();

private:

  /// The additional attributes for the "THREADALL" command.
  enum {
    TOLQ = CorrectionBase::SIZE, // Tolerance for the closed orbit positions.
    TOLP,        // Tolerance for the closed orbit angles.
    LISTC,       // List the correctors after correction.
    LISTM,       // List the monitors after correction.
    SIZE
  };

  // Not implemented.
  ThreadAll(const ThreadAll &);
  void operator=(const ThreadAll &);

  /// Clone constructor.
  ThreadAll(const string &name, ThreadAll *parent);

  // Backtrack to correct in one plane.
  virtual void correct(int plane, TLine::iterator &);

  /// Find the closed orbit by threading.
  bool thread(double tolq, double tolp);

  // The tracker.
  OwnPtr<OrbitTracker> itsTracker;
};

#endif // OPAL_ThreadAll_HH
