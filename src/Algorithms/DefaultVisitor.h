//
// Class DefaultVisitor
//   The default interface for a BeamlineVisitor.
//   A default implementation for all visitors that can iterate over a
//   beam line representation.
//   This abstract base class implements the default behaviour for the
//   structural classes Beamline and FlaggedElmPtr.
//   It also holds the data required for all visitors in a protected area.
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
#ifndef OPALX_DefaultVisitor_HH
#define OPALX_DefaultVisitor_HH

#include "AbsBeamline/BeamlineVisitor.h"

class ElementBase;

class DefaultVisitor : public BeamlineVisitor {
public:
    /// Constructor.
    //  Arguments:
    //  [ol]
    //  [li]The beamline to be used.
    //  [li]If true, the beam runs backwards through the line.
    //  [li]If true, we track against the beam.
    //  [/ol]
    DefaultVisitor(const Beamline& beamline, bool backBeam, bool backTrack);

    ~DefaultVisitor() override = default;

    /// Apply the algorithm to the top-level beamline.
    void execute() override;

    void visitComponent(const Component&) override;

    /// Apply the algorithm to a beam line.
    void visitBeamline(const Beamline&) override;

    /// Apply the algorithm to a constant E-field cavity element.
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override;

    /// Apply the algorithm to a drift space.
    void visitDrift(const Drift&) override;

    /// Apply the algorithm to a FlaggedElmPtr.
    void visitFlaggedElmPtr(const FlaggedElmPtr&) override;

    /// Apply the algorithm to a marker.
    void visitMarker(const Marker&) override;

    /// Apply the algorithm to a beam position monitor.
    void visitMonitor(const Monitor&) override;

    /// Apply the algorithm to a multipole.
    void visitMultipole(const Multipole&) override;

    /// Apply the algorithm to an arbitrary multipole.
    void visitMultipoleT(const MultipoleT&) override;

    /// Apply the algorithm to a Ring.
    void visitRing(const Ring&) override;

    /// Apply the algorithm to a RF cavity.
    void visitRFCavity(const RFCavity&) override;

    /// Apply the algorithm to a Solenoid.
    void visitSolenoid(const Solenoid&) override;

    /// Apply the algorithm to a traveling wave.
    void visitTravelingWave(const TravelingWave&) override;

    /// Apply the algorithm to a scaling FFA magnet.
    void visitScalingFFAMagnet(const ScalingFFAMagnet& spiral) override;

    /// Apply the algorithm to a vertical FFA magnet.
    void visitVerticalFFAMagnet(const VerticalFFAMagnet&) override;

    /// Apply the algorithm to a variable RF cavity.
    void visitVariableRFCavity(const VariableRFCavity&) override;

    /// Apply the algorithm to a Probe
    void visitProbe(const Probe& prob) override;

protected:
    // The top level beamline.
    const Beamline& itsLine;

    // The direction flags and corresponding factors.
    bool back_beam;   // true, if beam runs from right (s=C) to left (s=0).
    bool back_track;  // true, if tracking opposite to the beam direction.
    bool back_path;   // true, if tracking from right (s=C) to left (s=0).
    // back_path = back_beam && ! back_track || back_track && ! back_beam.

    double flip_B;  // set to -1.0 to flip B fields, when back_beam is true.
    double flip_s;  // set to -1.0 to flip direction of s,
    // when back_path is true.

private:
    // Not implemented.
    DefaultVisitor() = delete;
    DefaultVisitor(const DefaultVisitor&) = delete;
    void operator=(const DefaultVisitor&) = delete;

    // Default do-nothing routine.
    virtual void applyDefault(const ElementBase&);

    // The element order flag. Initially set to back_path.
    // This flag is reversed locally for reflected beam lines.
    bool local_flip;
};

#endif  // OPALX_DefaultVisitor_HH
