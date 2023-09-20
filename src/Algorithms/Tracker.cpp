//
// Class Tracker
//  Track particles or bunches.
//  An abstract base class for all visitors capable of tracking particles
//  through a beam element.
//  [P]
//  Phase space coordinates (in this order):
//  [DL]
//  [DT]x:[DD]
//    horizontal displacement (metres).
//  [DT]p_x/p_r:[DD]
//     horizontal canonical momentum (no dimension).
//  [DT]y:[DD]
//    vertical displacement (metres).
//  [DT]p_y/p_r:[DD]
//    vertical canonical momentum (no dimension).
//  [DT]delta_p/p_r:[DD]
//    relative momentum error (no dimension).
//  [DT]v*delta_t:[DD]
//    time difference delta_t w.r.t. the reference frame which moves with
//    uniform velocity
//  [P]
//    v_r = c*beta_r = p_r/m
//  [P]
//    along the design orbit, multiplied by the instantaneous velocity v of
//    the particle (metres).
//  [/DL]
//  Where
//  [DL]
//  [DT]p_r:[DD]
//    is the constant reference momentum defining the reference frame velocity.
//  [DT]m:[DD]
//    is the rest mass of the particles.
//  [/DL]
//  Other units used:
//  [DL]
//  [DT]reference momentum:[DD]
//    electron-volts.
//  [DT]accelerating voltage:[DD]
//    volts.
//  [DT]separator voltage:[DD]
//    volts.
//  [DT]frequencies:[DD]
//    hertz.
//  [DT]phase lags:[DD]
//    multiples of (2*pi).
//  [/DL]
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
#include "Algorithms/Tracker.h"
#include "Fields/BMultipoleField.h"

// FIXME Remove headers and dynamic_cast in readOneBunchFromFile
#include "Algorithms/PartBunch.hpp"

#include <cfloat>
#include <cmath>
#include <limits>

// Class Tracker
// ------------------------------------------------------------------------

Tracker::Tracker(const Beamline& beamline, const PartData& reference, bool backBeam, bool backTrack)
    : Tracker(beamline, nullptr, reference, backBeam, backTrack) {
}

Tracker::Tracker(
    const Beamline& beamline, PartBunch_t* bunch, const PartData& reference, bool backBeam,
    bool backTrack)
    : AbstractTracker(beamline, reference, backBeam, backTrack),
      itsBeamline_m(beamline),
      itsBunch_m(bunch) {
}

Tracker::~Tracker() {
}

const PartBunch_t* Tracker::getBunch() const {
    return itsBunch_m;
}

void Tracker::addToBunch(const OpalParticle& part) {
}

//~ void Tracker::setBunch(const PartBunch &bunch) {
//~ itsBunch_m = &bunch;
//~ }

void Tracker::visitComponent(const Component& comp) {
    comp.trackBunch(itsBunch_m, itsReference, back_beam, back_track);
}
