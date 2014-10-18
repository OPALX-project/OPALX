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
#ifdef HAVE_AMR_SOLVER
	#include <Amr.H>
	#include <ParallelDescriptor.H>
	#include <fstream>
	#include <iomanip>
	#include <iostream>
	#include <string>
	#include <sstream>
	#include <algorithm>
	#include <iterator>
	#include <utility>
#endif


class Beam;
class OpalData;
class DataSink;
class Distribution;
class Tracker;
class ParallelTTracker;
class FieldSolver;

// Class TrackRun
// ------------------------------------------------------------------------
/// The RUN command.

class TrackRun: public Action {

public:

    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /// Make clone.
    virtual TrackRun *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    TrackRun(const TrackRun &);
    void operator=(const TrackRun &);

    // Clone constructor.
    TrackRun(const std::string &name, TrackRun *parent);

    double SetDistributionParallelT(Beam *beam);
    ParallelTTracker *setupForAutophase();

    // Pointer to tracking algorithm.
    Tracker *itsTracker;

    Distribution *dist;

    std::vector<Distribution *> distrs_m;

    FieldSolver  *fs;

    DataSink *ds;

    OpalData *OPAL;
#ifdef HAVE_AMR_SOLVER
    std::vector<std::string>  filterString(std::string str);
    std::pair<Box,unsigned int> getBlGrids(std::string str);
#endif
};

#endif // OPAL_TrackRun_HH
