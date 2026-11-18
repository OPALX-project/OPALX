//
// Class TrackParser
//   The parser class used by the OPAL tracking module.
//   As long as control remains in this class, OPAL recognizes only the
//   commands allowed in tracking mode: RUN and ENDTRACK.
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
#ifndef OPAL_TrackParser_HH
#define OPAL_TrackParser_HH

#include "OpalParser/OpalParser.h"

class TrackParser : public OpalParser {
public:
    TrackParser();
    virtual ~TrackParser();

protected:
    /// A new RUN or ENDTRACK command named @p newName; nullptr for any other @p name.
    virtual Object* makeFrom(const std::string& name, const std::string& newName) const;

private:
    // Not implemented.
    TrackParser(const TrackParser&);
    void operator=(const TrackParser&);
};

#endif  // OPAL_TrackParser_HH
