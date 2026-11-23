//
// Class BeamStops
//   The stop rules of a run: DT segment switch, ZSTOP planes and the end of the run.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "Algorithms/BeamStop.h"

#include "Algorithms/StepSizeConfig.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Utilities/Util.h"

#include <functional>

extern Inform* gmsg;

void BeamStops::announce(std::size_t nContainers) {
    stops_m.resize(nContainers);
    for (size_t i = 0; i < stops_m.size(); ++i) {
        const BeamStop& stop = stops_m[i];
        if (stop.zStop) {
            *gmsg << level1 << "* Beam " << i << " stops when its last particle passes ZSTOP = "
                  << Util::getLengthString(*stop.zStop) << "." << endl;
        } else if (stop.untilMaxSteps) {
            *gmsg << level1 << "* Beam " << i << " is tracked until MAXSTEPS." << endl;
        }
    }
}

void BeamStops::keepDoneBeamsStopped(PartBunch_t& bunch) const {
    for (size_t i = 0; i < stops_m.size(); ++i) {
        if (stops_m[i].done) bunch.setPcAtSStop(i);
    }
}

BeamStops::StopRules BeamStops::applyStopRules(
        PartBunch_t& bunch, const StepSizeConfig& stepSizes, bool referenceHit) {
    Inform m("BeamStops::applyStopRules");
    const auto& particleContainers = bunch.getParticleContainers();

    // The primary reference switches to the next DT segment at each ZSTOP entry before the
    // last, and a beam with a ZSTOP stops once its last particle has passed its stop plane.
    StopRules stop;
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc || !bunch.isPcActive(i)) {
            continue;
        }
        ippl::Vector<double, 3> pdivg = pc->getRefPartP() / Util::getGamma(pc->getRefPartP());
        double beta                   = euclidean_norm(pdivg);
        double driftPerTimeStep       = std::abs(bunch.getdT()) * Physics::c * beta;
        m << level4 << "Calculated drift per time step (container " << i
          << "): " << Util::getLengthString(driftPerTimeStep) << "." << endl;

        if (i == 0 && !stepSizes.isLastSegment()
            && beam_stop::reached(pc->get_sPos(), stepSizes.getSStop(), driftPerTimeStep)) {
            m << level2
              << "Primary reference reached sStop = " << Util::getLengthString(stepSizes.getSStop())
              << "; switching to the next DT segment." << endl;
            stop.switchSegment = true;
        }
        if (i < stops_m.size() && !bunch.isFollower(i)) {
            BeamStop& beamStop = stops_m[i];
            if (beamStop.zStop && !beamStop.planeSet
                && beam_stop::reached(pc->get_sPos(), *beamStop.zStop, driftPerTimeStep)) {
                beam_stop::setPlane(beamStop, pc->getRefPartR(), pc->getRefPartP(), pc->get_sPos());
                *gmsg << level1 << "* Beam " << i
                      << ": reference reached ZSTOP = " << Util::getLengthString(*beamStop.zStop)
                      << " at t = " << Util::getTimeString(bunch.getT())
                      << "; tracking until its last particle passes." << endl;
            }
        }
    }
    updateBeamStops(bunch);

    ippl::Comm->allreduce(referenceHit, 1, std::logical_or<bool>());
    if (referenceHit) {
        *gmsg << level1 << "* A reference particle hit an element; the run ends." << endl;
    }
    const bool beamsDone = beam_stop::runFinished(stops_m);
    if (beamsDone) {
        *gmsg << level1 << "* Every beam with a ZSTOP is done; the run ends." << endl;
    }
    stop.runFinished = referenceHit || beamsDone;
    return stop;
}

void BeamStops::updateBeamStops(PartBunch_t& bunch) {
    for (size_t i = 0; i < stops_m.size(); ++i) {
        BeamStop& stop = stops_m[i];
        if (!stop.planeSet || stop.done) continue;

        // Particles are in the bunch frame: move the lab plane into it.
        double minimum = std::numeric_limits<double>::max();
        if (const auto pc = bunch.getParticleContainer(i)) {
            const CoordinateSystemTrafo& toLab = pc->getToLabTrafo();
            minimum                            = beam_stop::minDistanceToPlane(
                    pc->R.getView(), pc->getLocalNum(), toLab.transformFrom(stop.planeOrigin),
                    toLab.rotateFrom(stop.planeNormal));
        }
        ippl::Comm->allreduce(minimum, 1, std::less<double>());
        if (minimum >= 0.0) {
            stop.done = true;
            bunch.setPcAtSStop(i);
            *gmsg << level1 << "* Beam " << i
                  << ": its last particle passed ZSTOP at t = " << Util::getTimeString(bunch.getT())
                  << "; it is no longer tracked." << endl;
        }
    }
}
