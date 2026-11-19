// ------------------------------------------------------------------------
// $RCSfile: TrackEnd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackEnd
//   The class for the OPAL ENDTRACK command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackEnd.h"
#include "Track/TrackCmd.h"

// Class TrackEnd
// ------------------------------------------------------------------------

TrackEnd::TrackEnd() : Action(0, "ENDTRACK", "The \"ENDTRACK\" command closes the open TRACK.") {}

TrackEnd::~TrackEnd() {}

void TrackEnd::execute() { TrackCmd::closeTrack(); }
