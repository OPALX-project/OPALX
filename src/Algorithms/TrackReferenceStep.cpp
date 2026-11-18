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
                bool hit = false;
                for (const auto& element : beamline.getElements(r)) {
                    if (!diagnostics && element->getType() == ElementType::MONITOR) continue;
                    const auto& transform            = beamline.getCSTrafoLab2Local(element);
                    const Vector_t<double, 3> localR = transform.transformTo(r);
                    const Vector_t<double, 3> localP = transform.rotateTo(p);
                    Vector_t<double, 3> localE(0), localB(0);
                    hit = element->applyToReferenceParticle(localR, localP, t, localE, localB)
                          || hit;
                    electric += transform.rotateFrom(localE);
                    magnetic += transform.rotateFrom(localB);
                }
                return hit;
            });
}
