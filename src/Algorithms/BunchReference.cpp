//
// Class BunchReference
//   The reference particles of the containers of a bunch: their start, the reference pass of the
//   primary beam and their advance in each step.
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
#include "Algorithms/BunchReference.h"

#include "AbsBeamline/RFCavity.h"
#include "Algorithms/ReferenceTracer.h"
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/TrackReferenceStep.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"
#include "Moments/BeamMoments.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <utility>

extern Inform* gmsg;

using opalx::moments::Moment;

BunchReference::BunchReference(
        PartBunch_t& bunch, OpalBeamline& beamline,
        std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta, int autoPhase,
        CavityPhases cavityPhases)
    : bunch_m(bunch),
      beamline_m(beamline),
      beamReferenceMomenta_m(std::move(beamReferenceMomenta)),
      autoPhase_m(autoPhase),
      cavityPhases_m(std::move(cavityPhases)),
      referencePassTimer_m(IpplTimings::getTimer("ReferencePass")) {}

void BunchReference::setInitialPoses(
        const std::vector<std::vector<std::shared_ptr<SamplingBase>>>& emittingSamplers) {
    Inform m("BunchReference::setInitialPoses");
    const auto& particleContainers      = bunch_m.getParticleContainers();
    CoordinateSystemTrafo beamlineToLab = beamline_m.getCSTrafoLab2Local().inverted();
    for (size_t ci = 0; ci < particleContainers.size(); ++ci) {
        const auto& pc = particleContainers[ci];
        if (!pc || bunch_m.isFollower(ci)) {
            continue;  // a decay daughter container takes its root's reference below
        }
        if (!pc->getReference()) {
            throw OpalException(
                    "BunchReference::setInitialPoses",
                    "Particle container has null PartData reference during lab-frame init.");
        }
        pc->setToLabTrafo(beamlineToLab);

        // Resolve the reference particle's pose in the beamline frame. Position and
        // momentum follow the same 3-tier rule: the bunch mean when particles already
        // exist; otherwise the input-specified emission offsets (R0, P0) reported by
        // the sampler; otherwise the design pose (lattice origin, beta*gamma along +z).
        Vector_t<double, 3> refR = 0.0;
        Vector_t<double, 3> refP = 0.0;
        if (pc->getTotalNum() > 0) {
            const auto moments =
                    opalx::moments::compute(*pc, Moment::MeanPosition | Moment::MeanMomentum);
            refR = moments.meanR();
            refP = moments.meanP();
        } else {
            // Empty container (e.g. an emitted distribution before its first emission):
            // the bunch mean is undefined, so take the emission offsets from the sampler.
            bool useSamplerPosition         = false;
            bool useSamplerMomentum         = false;
            Vector_t<double, 3> samplerRefR = 0.0;
            Vector_t<double, 3> samplerRefP = 0.0;
            if (ci < emittingSamplers.size()) {
                for (const auto& sampler : emittingSamplers[ci]) {
                    if (!sampler) {
                        continue;
                    }
                    if (!useSamplerMomentum && sampler->hasInitialReferenceMomentum()) {
                        samplerRefP        = sampler->getInitialReferenceMomentum();
                        useSamplerMomentum = true;
                    }
                    if (!useSamplerPosition && sampler->hasInitialReferencePosition()) {
                        samplerRefR        = sampler->getInitialReferencePosition();
                        useSamplerPosition = true;
                    }
                }
            }

            // Position: emission offset R0, else the lattice origin.
            refR = useSamplerPosition ? samplerRefR : Vector_t<double, 3>(0.0);

            // Momentum: emission offset P0, else the design beta*gamma along +z.
            if (useSamplerMomentum) {
                if (dot(samplerRefP, samplerRefP) <= 0.0) {
                    throw OpalException(
                            "BunchReference::setInitialPoses",
                            "Sampler-provided initial reference momentum is zero.");
                }
                refP = samplerRefP;
            } else {
                const PartData& pref = *pc->getReference();
                const double P0      = pref.getP() / pref.getM();  // beta*gamma from BEAM pc
                refP                 = Vector_t<double, 3>(0.0, 0.0, P0);
            }
        }

        // A file distribution with PC, ENERGY or GAMMA on BEAM takes that momentum along
        // REFDIRECTION. The mean of a wide angular spread is much shorter than the particle
        // momenta.
        if (ci < beamReferenceMomenta_m.size() && beamReferenceMomenta_m[ci]) {
            refP = *beamReferenceMomenta_m[ci];
        }

        pc->getRefPartR() = beamlineToLab.transformTo(refR);
        pc->getRefPartP() = beamlineToLab.rotateTo(refP);
    }

    m << level4
      << "Transformed reference particle position and momentum to lab frame (all containers)."
      << endl;
}

bool BunchReference::findStartPositions(
        double sStart, const StepSizeConfig& stepSizes,
        const std::function<void(double)>& selectTimeStep) {
    // Primary container (index 0) drives segment advances.
    auto primary = bunch_m.getParticleContainer(0);

    if (sStart <= primary->get_sPos()) {
        return false;
    }

    StepSizeConfig stepSizesCopy(stepSizes);

    double t = 0.0;
    bunch_m.setT(t);

    selectTimeStep(stepSizesCopy.getdT());

    const auto& containers = bunch_m.getParticleContainers();

    // The last, fractional step takes every reference to sTarget; the bunch time advances by the
    // primary's share.
    const auto applyLastStep = [&](double sTarget, double speed) {
        const double tauPrimary = (sTarget - primary->get_sPos()) / speed;
        bunch_m.setT(bunch_m.getT() + tauPrimary);
        for (size_t i = 0; i < containers.size(); ++i) {
            const auto& pc = containers[i];
            if (!pc || bunch_m.isFollower(i)) {
                continue;
            }
            Vector_t<double, 3> pv =
                    pc->getRefPartP() * Physics::c / Util::getGamma(pc->getRefPartP());
            double speed_i = euclidean_norm(pv);
            double tau_i   = (sTarget - pc->get_sPos()) / speed_i;
            pc->applyFractionalStep(tau_i, sStart);
        }
    };

    while (true) {
        autophaseCavities();

        t += bunch_m.getdT();
        bunch_m.setT(t);

        std::vector<Vector_t<double, 3>> oldRs(containers.size());
        for (size_t i = 0; i < containers.size(); ++i) {
            if (containers[i]) {
                oldRs[i] = containers[i]->getRefPartR();
            }
        }

        updateReferenceParticles();

        for (size_t i = 0; i < containers.size(); ++i) {
            if (!containers[i] || bunch_m.isFollower(i)) {
                continue;
            }
            Vector_t<double, 3> dR = containers[i]->getRefPartR() - oldRs[i];
            containers[i]->set_sPos(containers[i]->get_sPos() + euclidean_norm(dR));
        }

        Vector_t<double, 3> tmp =
                primary->getRefPartP() * Physics::c / Util::getGamma(primary->getRefPartP());
        double speed = euclidean_norm(tmp);

        if (primary->get_sPos() > stepSizesCopy.getSStop()) {
            ++stepSizesCopy;

            if (stepSizesCopy.reachedEnd()) {
                --stepSizesCopy;
                applyLastStep(stepSizesCopy.getSStop(), speed);
                break;
            }

            selectTimeStep(stepSizesCopy.getdT());
        }

        if (std::abs(primary->get_sPos() - sStart) <= 0.5 * bunch_m.getdT() * speed) {
            applyLastStep(sStart, speed);
            break;
        }
    }

    // An element the reference hits before the start does not end the run.
    hitElement_m = false;
    return true;
}

void BunchReference::tracePrimary(const StepSizeConfig& stepSizes, double minTimeStep) {
    Inform m("BunchReference::tracePrimary");

    // Global spatial bounds: union over all containers
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (bunch_m.getTotalNumAllContainers() > 0) {
        std::vector<const ParticleContainer_t*> containers;
        for (const auto& pc : bunch_m.getParticleContainers()) {
            if (pc) {
                containers.push_back(pc.get());
            }
        }
        // Empty containers contribute nothing; without any particle the box is a point at the
        // origin.
        const auto bounds = opalx::moments::compute(containers, Moment::Bounds);
        rmin              = bounds.minR();
        rmax              = bounds.maxR();
        m << level4 << "Initial bunch bounds:\n"
          << "  rmin = " << rmin << "\n"
          << "  rmax = " << rmax << endl;
    }

    printInitialContainerRefs(m);

    IpplTimings::startTimer(referencePassTimer_m);

    // The primary beam (the first container with a reference) traces its reference particle once
    // through the lattice: autophasing, element-length and reach checks.
    const size_t nContainers = bunch_m.getNumParticleContainers();
    for (size_t ci = 0; ci < nContainers; ++ci) {
        const auto& pc = bunch_m.getParticleContainer(ci);
        if (!pc || !pc->getReference() || bunch_m.isFollower(ci)) {
            continue;  // a decay daughter container follows its root's reference
        }
        ReferenceTracer tracer(
                *pc->getReference(), beamline_m, stepSizes, minTimeStep, autoPhase_m,
                &cavityPhases_m);
        ReferenceTracer::RayState start;
        start.position   = pc->getRefPartR();
        start.momentum   = pc->getRefPartP();
        start.time       = bunch_m.getT();
        start.pathLength = pc->get_sPos();
        tracer.execute(start, -rmin(2));
        break;  // only the primary beam
    }
    m << level4 << "Reference pass done." << endl;
    IpplTimings::stopTimer(referencePassTimer_m);
}

void BunchReference::advance() {
    Inform m("BunchReference::advance");
    if (bunch_m.getT() > 0.0) {
        updateReferenceParticles();
        updateRefToLabCSTrafo();
    }

    const auto& particleContainers = bunch_m.getParticleContainers();
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc || !bunch_m.isPcActive(i)) {
            continue;
        }
        m << level4 << "Current path length (container " << i << ") is " << pc->get_sPos() << "."
          << endl;
    }
}

void BunchReference::copyFollowers() {
    const size_t n = bunch_m.getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!bunch_m.isFollower(i)) {
            continue;
        }
        auto pc     = bunch_m.getParticleContainer(i);
        auto rootPc = bunch_m.getParticleContainer(bunch_m.getPcRoot(i));
        if (pc && rootPc) {
            pc->copyReferenceOf(*rootPc);
        }
    }
}

void BunchReference::updateReferenceParticles() {
    const double dt = std::min(bunch_m.getT(), bunch_m.getdT());

    const size_t n = bunch_m.getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        // A decay daughter container has no reference particle of its own. Its root keeps
        // moving while any of its daughter containers is tracked, even when the root is empty.
        if (bunch_m.isFollower(i) || !bunch_m.referenceMoves(i)) {
            continue;
        }
        auto pcPtr = bunch_m.getParticleContainer(i);
        if (!pcPtr) {
            continue;
        }
        auto& pc                = *pcPtr;
        const PartData& refKick = *pc.getReference();
        const auto advance      = track_reference::advanceInBeamline;
        const auto end =
                advance(beamline_m, refKick, {pc.getRefPartR(), pc.getRefPartP()}, dt,
                        bunch_m.getT(), true);
        pc.getRefPartR() = end.position;
        pc.getRefPartP() = end.momentum;
        if (end.hitMaterial) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            hitElement_m = true;
        }
    }
}

void BunchReference::updateRefToLabCSTrafo() {
    // Map the lab-frame reference position and momentum into the current bunch frame, then
    // move and turn the bunch frame onto them. Momentum is only rotated, not translated.
    const double bunchDT = bunch_m.getdT();
    const size_t n       = bunch_m.getNumParticleContainers();
    std::vector<CoordinateSystemTrafo> updates(n);  // identity unless the container moved
    for (size_t i = 0; i < n; ++i) {
        if (bunch_m.isFollower(i) || !bunch_m.referenceMoves(i)) {
            continue;
        }
        auto pc = bunch_m.getParticleContainer(i);
        if (pc) {
            updates[i] = pc->updateRefToLabCSTrafo(bunchDT);
        }
    }
    // Decay daughter containers take the same frame change as their root, active or not, so
    // daughters created in this step end up in the root's new frame too.
    for (size_t i = 0; i < n; ++i) {
        if (!bunch_m.isFollower(i)) {
            continue;
        }
        const size_t root = bunch_m.getPcRoot(i);
        auto pc           = bunch_m.getParticleContainer(i);
        auto rootPc       = bunch_m.getParticleContainer(root);
        if (pc && rootPc) {
            pc->followFrameOf(*rootPc, updates[root]);
        }
    }
}

void BunchReference::autophaseCavities() {
    const BorisPusher pusher;
    const PartData& ref = *bunch_m.getParticleContainer()->getReference();
    double t            = bunch_m.getT();
    Vector_t<double, 3> nextR =
            bunch_m.getParticleContainer()->getRefPartR() / (Physics::c * bunch_m.getdT());
    pusher.push(nextR, bunch_m.getParticleContainer()->getRefPartP(), bunch_m.getdT());
    nextR *= Physics::c * bunch_m.getdT();

    auto elementSet = beamline_m.getElements(nextR);
    for (auto element : elementSet) {
        // A TravelingWave is an RFCavity.
        if (element->getType() != ElementType::TRAVELINGWAVE
            && element->getType() != ElementType::RFCAVITY) {
            continue;
        }
        const RFCavity* cavity = static_cast<const RFCavity*>(element.get());
        if (!cavity->getAutophaseVeto()) {
            CavityAutophaser ap(ref, element, autoPhase_m, &cavityPhases_m);
            ap.getPhaseAtMaxEnergy(
                    beamline_m.transformToLocalCS(
                            element, bunch_m.getParticleContainer()->getRefPartR()),
                    beamline_m.rotateToLocalCS(
                            element, bunch_m.getParticleContainer()->getRefPartP()),
                    t, bunch_m.getdT());
        }
    }
}

void BunchReference::printInitialContainerRefs(Inform& m) const {
    const auto& particleContainers = bunch_m.getParticleContainers();
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc) {
            m << level3 << "ParallelTrack: container " << i << " is null." << endl;
            continue;
        }
        if (pc->getTotalNum() > 0) {
            // Collective; evaluated on every rank before the rank-filtered log line.
            const double meanPz = opalx::moments::compute(*pc, Moment::MeanMomentum).meanP()(2);
            m << level3 << "ParallelTrack (container " << i << "): momentum z = " << meanPz << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        } else {
            m << level3 << "ParallelTrack: container " << i
              << " empty; total particles = " << bunch_m.getTotalNumAllContainers() << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        }
    }
}
