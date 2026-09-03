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

#include <algorithm>

namespace {
    void cloneAndAssignMembership(FlaggedBeamline& line, const std::string& ringName) {
        for (auto& member : line) {
            ElementBase* element = member.getElement();
            if (!element || (!element->getName().empty() && element->getName().front() == '#')) {
                continue;
            }

            member.setElement(element->clone());
            element = member.getElement();
            element->setBeamlineMembership(BeamlineTopology::RING, ringName);

            if (auto* nestedLine = dynamic_cast<FlaggedBeamline*>(element)) {
                cloneAndAssignMembership(*nestedLine, ringName);
            }
        }
    }
}  // namespace

Ring::Ring()
    : Line("RING",
           "The \"RING\" statement defines a closed-topology beamline list.\n"
           "\t<name> : ring = (<list>)") {
    fetchLine()->setBeamlineMembership(BeamlineTopology::RING, "RING");
}

Ring::Ring(const std::string& name, Ring* parent) : Line(name, parent) {
    fetchLine()->setBeamlineMembership(BeamlineTopology::RING, name);
}

Ring* Ring::clone(const std::string& name) { return new Ring(name, this); }

Ring* Ring::copy(const std::string& name) {
    Ring* clone              = new Ring(name, this);
    FlaggedBeamline* oldLine = fetchLine();
    FlaggedBeamline* newLine = clone->fetchLine();
    std::copy(oldLine->begin(), oldLine->end(), std::back_inserter(*newLine));
    return clone;
}

void Ring::parse(Statement& stat) { parseDefinition(stat, false); }

const char* Ring::getSequenceKeyword() const { return "RING"; }

void Ring::prepareForTracking() {
    FlaggedBeamline* line = fetchLine();
    line->setBeamlineMembership(BeamlineTopology::RING, getOpalName());
    cloneAndAssignMembership(*line, getOpalName());
}
