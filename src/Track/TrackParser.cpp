//
// Class TrackParser
//   The parser class used by the OPAL tracking module.
//   As long as control remains in this class, OPAL recognizes only the
//   commands allowed in tracking mode.  Thus this parser has its own
//   command directory with a find() method which is used to find commands.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Track/TrackParser.h"
#include "Track/TrackEnd.h"
#include "Track/TrackRun.h"

TrackParser::TrackParser() {}

TrackParser::~TrackParser() {}

Object* TrackParser::makeFrom(const std::string& name, const std::string& newName) const {
    Object* command = nullptr;
    if (name == "RUN") {
        command = new TrackRun();
    } else if (name == "ENDTRACK") {
        command = new TrackEnd();
    } else {
        return nullptr;
    }
    command->setOpalName(newName);
    return command;
}
