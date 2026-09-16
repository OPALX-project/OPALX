// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/TrackReferenceStep.h"
#include "Algorithms/PartData.h"
#include "Algorithms/ExternalFieldRayTracker.h"
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

track_reference::State track_reference::advanceResolvedInBeamline(
        OpalBeamline& beamline, const PartData& reference, State state, double dt, double endTime,
        bool diagnostics) {
    if (!std::isfinite(dt) || dt == 0 || !std::isfinite(endTime))
        throw std::invalid_argument("Invalid resolved TRACK step duration/time");
    ExternalFieldRayTracker tracker(beamline, reference);
    ExternalFieldRayTracker::State start;
    start.position = state.position;
    start.momentum = state.momentum;
    start.time = endTime - dt;
    std::vector<ExternalFieldRayTracker::Step> accepted;
    const auto end = tracker.advance(start, dt, diagnostics ? &accepted : nullptr);
    if (diagnostics) {
        for (const auto& step : accepted) {
            for (const auto& element : beamline.getElements(step.midpoint.position)) {
                if (element->getType() != ElementType::MONITOR) continue;
                Vector_t<double, 3> electric(0), magnetic(0);
                state.hitMaterial = element->applyToReferenceParticle(
                        beamline.transformToLocalCS(element, step.midpoint.position),
                        beamline.rotateToLocalCS(element, step.midpoint.momentum),
                        step.midpoint.time, electric, magnetic) || state.hitMaterial;
            }
        }
    }
    state.position = end.position;
    state.momentum = end.momentum;
    return state;
}

double track_reference::collectiveTerminalStep(
        MPI_Comm communicator, const std::function<double()>& search) {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    double duration = 0;
    std::string error;
    if (rank == 0) {
        try {
            duration = search();
            if (!(duration > 0) || !std::isfinite(duration))
                throw std::runtime_error("Invalid terminal-turn duration");
        } catch (const OpalException& exception) {
            error = exception.where() + ": " + exception.what();
        } catch (const std::exception& exception) {
            error = exception.what();
            if (error.empty()) error = "Terminal-turn reference trial failed";
        } catch (...) {
            error = "Unknown failure during terminal-turn reference trial";
        }
    }
    int length = static_cast<int>(error.size());
    MPI_Bcast(&length, 1, MPI_INT, 0, communicator);
    if (length > 0) {
        error.resize(length);
        MPI_Bcast(error.data(), length, MPI_CHAR, 0, communicator);
        throw OpalException("ParallelTracker::terminalTurn", error);
    }
    MPI_Bcast(&duration, 1, MPI_DOUBLE, 0, communicator);
    return duration;
}
