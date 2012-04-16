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

class OpalData;
class DataSink;
class Distribution;
class ParallelTTracker;

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

    ParallelTTracker *setupForAutophase();

    // Pointer to tracking algorithm.
    Tracker *itsTracker;

    Distribution *dist;

    std::vector<Distribution *> distrs_m;

    FieldSolver  *fs;

    DataSink *ds;

    OpalData *OPAL;

};

#endif // OPAL_TrackRun_HH
