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
