// ------------------------------------------------------------------------
// $RCSfile: TrackParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackParser
//   The parser class for the OPAL tracking module.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackParser.h"
#include "Track/TrackEnd.h"
#include "Track/TrackRun.h"
#include "Track/TrackSave.h"
#include "Track/TrackStart.h"


// Class TrackParser
// ------------------------------------------------------------------------


TrackParser::TrackParser():
    trackDirectory() {
    trackDirectory.insert("ENDTRACK", new TrackEnd());
    //trackDirectory.insert("NOISE",    new TrackNoise());
    trackDirectory.insert("RUN",      new TrackRun());
    trackDirectory.insert("TSAVE",    new TrackSave());
    trackDirectory.insert("START",    new TrackStart());
}


TrackParser::~TrackParser()
{}


Object *TrackParser::find(const string &name) const {
    return trackDirectory.find(name);
}
