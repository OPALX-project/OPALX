//
// Class BeamlineVisitor
//   The abstract class BeamlineVisitor is the base class for all visitors
//   (algorithms) that can iterator over a beam line representation.
//   A BeamlineVisitor applies itself to the representation via the
//   ``Visitor'' pattern, see
//   [p]
//   E. Gamma, R. Helm, R. Johnson, and J. Vlissides,
//   [BR]
//   Design Patterns, Elements of Reusable Object-Oriented Software.
//   [p]
//   By using only pure abstract classes as an interface between the
//   BeamlineVisitor and the beam line representation,
//   we decouple the former from the implementation details of the latter.
//   [p]
//   The interface is defined in such a way that a visitor cannot modify the
//   structure of a beam line, but it can assign special data like misalignments
//   or integrators without problems.
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
#ifndef CLASSIC_BeamlineVisitor_HH
#define CLASSIC_BeamlineVisitor_HH

// Generic element classes interacting with a BeamlineVisitor.
class Component;

// Beam line structure classes.
class Beamline;
class FlaggedElmPtr;

// Specific element classes interacting with a BeamlineVisitor
class Drift;
class Marker;
class Monitor;
class Multipole;
class MultipoleT;
class MultipoleTStraight;
class MultipoleTCurvedConstRadius;
class MultipoleTCurvedVarRadius;
class RFCavity;
class TravelingWave;

class BeamlineVisitor {

public:

    BeamlineVisitor();
    virtual ~BeamlineVisitor();

    /// Execute the algorithm on the attached beam line.
    virtual void execute() = 0;

    /// Apply the algorithm to a beam line.
    virtual void visitBeamline(const Beamline &) = 0;
    
    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &) = 0;

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift &) = 0;

    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &) = 0;

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &) = 0;

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor &) = 0;

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &) = 0;

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT &) = 0;

    /// Apply the algorithm to an arbitrary straight multipole.
    virtual void visitMultipoleTStraight(const MultipoleTStraight &) = 0;

    /// Apply the algorithm to an arbitrary curved multipole of constant radius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &) = 0;

    /// Apply the algorithm to an arbitrary curved multipole of variable radius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &) = 0;

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &) = 0;

    /// Apply the algorithm to a traveling wave.
    virtual void visitTravelingWave(const TravelingWave &) = 0;

private:

    // Not implemented.
    BeamlineVisitor(const BeamlineVisitor &);
    void operator=(const BeamlineVisitor &);
};


#endif // CLASSIC_BeamlineVisitor_HH
