// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/TrackReferenceStep.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

track_reference::State track_reference::advanceInBeamline(
        OpalBeamline& beamline, const PartData& reference, State state, double dt, double endTime,
        bool diagnostics) {
    return advance(
            state, dt, endTime, reference.getM(), reference.getQ(),
            [&](const auto& r, const auto& p, double t, auto& electric, auto& magnetic) {
                const ElementSet elements = beamline.getElements(r);
                const bool hit = beamline.applyOnHost(elements, r, t, electric, magnetic);
                if (diagnostics) beamline.recordReference(elements, r, p, t);
                return hit;
            });
}
