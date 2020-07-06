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


// Class TrackParser
// ------------------------------------------------------------------------


TrackParser::TrackParser():
    trackDirectory() {
    trackDirectory.insert("ENDTRACK", new TrackEnd());
    trackDirectory.insert("RUN",      new TrackRun());
}


TrackParser::~TrackParser()
{}


Object *TrackParser::find(const std::string &name) const {
    return trackDirectory.find(name);
}
