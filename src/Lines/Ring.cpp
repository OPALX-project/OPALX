//
// Class Ring
//   The RING definition.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
#include "Lines/Ring.h"

Ring::Ring()
    : Line("RING",
           "The \"RING\" statement defines a closed-topology beamline list.\n"
           "\t<name> : ring = (<list>)") {}

void Ring::parse(Statement& stat) { parseDefinition(stat, false); }

const char* Ring::getSequenceKeyword() const { return "RING"; }
