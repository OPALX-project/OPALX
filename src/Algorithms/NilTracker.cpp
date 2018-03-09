//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Algorithms/NilTracker.h"

NilTracker::NilTracker(const Beamline &beamline,
                       const PartData &reference,
                       bool revBeam,
                       bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack)
{ }


NilTracker::~NilTracker() {

}

void NilTracker::execute() {

}