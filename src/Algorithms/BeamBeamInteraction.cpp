/**
 * @file BeamBeamInteraction.cpp
 * @brief Runtime implementation of the BeamBeam collective element interaction.
 */

#include "Algorithms/BeamBeamInteraction.h"

#include "AbsBeamline/BeamBeam.h"
#include "Algorithms/OrbitThreader.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "Utilities/BeamBeamWindowAnimation.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"
#include "Utility/PAssert.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

extern Inform* gmsg;

BeamBeamInteraction::BeamBeamInteraction(const BeamBeam& element)
    : element_m(element),
      windowTimer_m(IpplTimings::getTimer("BB window total")),
      entryTransitionTimer_m(IpplTimings::getTimer("BB entry trans")),
      meshSetupTimer_m(IpplTimings::getTimer("BB mesh setup")),
      selfFieldTimer_m(IpplTimings::getTimer("BB self fields")),
      transformBackTimer_m(IpplTimings::getTimer("BB transform back")),
      witnessGatherTimer_m(IpplTimings::getTimer("BB witness gather")),
      transitionDumpTimer_m(IpplTimings::getTimer("BB transition dump")),
      state_m(),
      windowAnimation_m(std::make_unique<BeamBeamWindowAnimation>()) {}

BeamBeamInteraction::~BeamBeamInteraction() = default;

ElementInteractionResult BeamBeamInteraction::execute(
        ElementInteractionPhase phase, ElementInteractionContext& context) {
    Inform localMessage("BeamBeamInteraction ");
    Inform& message = context.message != nullptr ? *context.message : localMessage;

    switch (phase) {
        case ElementInteractionPhase::SelfField:
            return ElementInteractionResult{computeSelfFields(context)};
        case ElementInteractionPhase::AfterEmission:
            gatherFieldsToWitnessContainers(context.bunch, message);
            logDiagnostics(context.bunch);
            return {};
        case ElementInteractionPhase::Diagnostics:
            logDiagnostics(context.bunch);
            return {};
    }
    return {};
}

bool BeamBeamInteraction::freezesFieldMesh() const noexcept {
    return state_m.state == BEAMBEAM::WindowState::Active;
}

bool BeamBeamInteraction::suppressesDefaultSelfField() const noexcept {
    return state_m.sourceRetired;
}

bool BeamBeamInteraction::computeSelfFields(ElementInteractionContext& context) {
    if (context.sourceOrbitThreader == nullptr || context.referenceToBeamCSTrafo == nullptr
        || context.beamToReferenceCSTrafo == nullptr) {
        throw OpalException(
                "BeamBeamInteraction::computeSelfFields",
                "The self-field phase requires an orbit threader and both beam-frame "
                "coordinate transforms.");
    }

    Inform localMessage("BeamBeamInteraction::computeSelfFields ");
    Inform& message    = context.message != nullptr ? *context.message : localMessage;
    PartBunch_t& bunch = context.bunch;

    checkInRegion(context);
    if (!state_m.sourceRetired && state_m.geometry.has_value()
        && BEAMBEAM::sourceRetireTimeReached(
                bunch.getT(), state_m.geometry->config.sourceRetireTime)) {
        state_m.sourceRetirementPending = true;
    }
    if (state_m.sourceRetirementPending) {
        referenceToBeamCSTrafo_m.reset();
        if (state_m.state == BEAMBEAM::WindowState::Active) {
            leaveWindow(bunch, message);
        }
        retireSourceContainer(bunch, message);
        return true;
    }

    if (state_m.state != BEAMBEAM::WindowState::Active) {
        referenceToBeamCSTrafo_m.reset();
        return false;
    }

    referenceToBeamCSTrafo_m = *context.referenceToBeamCSTrafo;
    computeWindowSelfFields(*context.beamToReferenceCSTrafo, bunch, message);
    return true;
}

void BeamBeamInteraction::checkInRegion(ElementInteractionContext& context) {
    if (state_m.state == BEAMBEAM::WindowState::Completed) {
        return;
    }

    PartBunch_t& bunch = context.bunch;
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    bunch.calcBeamParameters();
    bunch.get_bounds(rmin, rmax);

    auto source                          = bunch.getParticleContainer();
    const double bunchS                  = source->get_sPos();
    const LongitudinalExtent bunchExtent = computeLongitudinalExtent(bunchS, rmin, rmax);

    std::optional<BEAMBEAM::ActualGeometry> geometry = detectWindow(context, rmin, rmax);
    if (!geometry.has_value() && state_m.state == BEAMBEAM::WindowState::Active
        && state_m.geometry.has_value()) {
        geometry = state_m.geometry;
    }
    if (!geometry.has_value()) {
        diagnostics_m.frameObserved  = false;
        state_m.sourceBunchesOverlap = false;
        return;
    }

    const auto& activeGeometry  = *geometry;
    diagnostics_m.frameObserved = activeGeometry.config.visualize
                                  && bunchExtent.head >= activeGeometry.beginS
                                  && bunchExtent.tail <= activeGeometry.endS;
    const bool copyModelActive =
            BEAMBEAM::copyTimeReached(bunch.getT(), activeGeometry.config.copyTime);
    state_m.sourceBunchesOverlap = copyModelActive
                                   && BEAMBEAM::copiedSourceBunchesOverlap(
                                           bunchExtent.tail, bunchExtent.head, activeGeometry);
    state_m.geometry = activeGeometry;

    const double cellHalfWidth = 0.5 * activeGeometry.length / static_cast<double>(bunch.nr_m[2]);
    const bool leavingWindow   = state_m.state == BEAMBEAM::WindowState::Active
                               && bunchExtent.head > activeGeometry.endS;
    Inform message("BeamBeam ");
    if (leavingWindow) {
        leaveWindow(bunch, message);
    }

    if (diagnostics_m.frameObserved) {
        renderWindowFrame(bunchExtent.tail, bunchExtent.head, activeGeometry);
    }

    if (state_m.state == BEAMBEAM::WindowState::Inactive
        && bunchExtent.tail >= activeGeometry.beginS + cellHalfWidth) {
        enterWindow(activeGeometry, bunch, message);
    }
}

BeamBeamInteraction::LongitudinalExtent BeamBeamInteraction::computeLongitudinalExtent(
        double bunchS, const ippl::Vector<double, 3>& rmin,
        const ippl::Vector<double, 3>& rmax) const {
    return {bunchS + rmin[2], bunchS + rmax[2]};
}

std::optional<BEAMBEAM::ActualGeometry> BeamBeamInteraction::detectWindow(
        ElementInteractionContext& context, const ippl::Vector<double, 3>& rmin,
        const ippl::Vector<double, 3>& rmax) {
    PartBunch_t& bunch            = context.bunch;
    OrbitThreader& sourceThreader = *context.sourceOrbitThreader;
    const double bunchS           = bunch.getParticleContainer()->get_sPos();

    IndexMap::value_t elements;
    try {
        elements = sourceThreader.query(bunchS + 0.5 * (rmax[2] + rmin[2]), rmax[2] - rmin[2]);
    } catch (IndexMap::OutOfBounds&) {
        if (context.endOfLine != nullptr) {
            *context.endOfLine = true;
        }
        return std::nullopt;
    }

    for (const auto& element : elements) {
        if (element.get() != &element_m) {
            continue;
        }

        const double windowLength = element_m.getGeometry().getElementLength();
        if (windowLength <= 0.0) {
            return std::nullopt;
        }

        const IndexMap::key_t range = sourceThreader.getRange(element);
        // The OrbitThreader map is only a discovery aid. Its range can be clipped by the
        // configured tracking horizon (for example, a short MAXSTEPS run that starts inside this
        // element), so using range.end as the physical window end makes the BeamBeam lifecycle
        // depend on the sampled bunch extent. ELEMEDGE and the element length are the authoritative
        // path-length geometry. For 6D-pose placement, where no ELEMEDGE path coordinate exists,
        // retain the threaded entrance as the best available path-length anchor.
        const double windowBeginS =
                element_m.isElementPositionSet() ? element_m.getElementPosition() : range.begin;
        const double windowEndS                  = windowBeginS + windowLength;
        const double interactionPointS =
                BEAMBEAM::interactionPointAtElementMidpoint(windowBeginS, windowEndS);
        std::optional<double> xAperture;
        std::optional<double> yAperture;
        const auto aperture = element_m.getAperture();
        if (element_m.getAttribute("APERTURE_SET") != 0.0 && aperture.second.size() >= 2) {
            xAperture = std::abs(aperture.second[0]);
            yAperture = std::abs(aperture.second[1]);
        }

        std::optional<double> sourceRetireTime;
        const double retireTime = element_m.getAttribute("RETIRE_TIME");
        if (retireTime > 0.0) {
            sourceRetireTime = retireTime;
        }

        std::optional<double> copyTime;
        const double copyTimeValue = element_m.getAttribute("COPY_TIME");
        if (copyTimeValue > 0.0) {
            copyTime = copyTimeValue;
        }

        return BEAMBEAM::ActualGeometry{
                interactionPointS, windowBeginS, windowEndS, windowLength,
                BEAMBEAM::Config{
                        element_m.getAttribute("VISUALIZE") != 0.0,
                        element_m.getAttribute("BBRIGID") != 0.0, copyTime, sourceRetireTime,
                        xAperture, yAperture,
                        BEAMBEAM::decodeWitnessContainerMask(
                                element_m.getAttribute("WITNESS_CONTAINERS_MASK"))}};
    }

    return std::nullopt;
}

void BeamBeamInteraction::enterWindow(
        const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch, Inform& message) {
    state_m.state                        = BEAMBEAM::WindowState::Active;
    state_m.geometry                     = geometry;
    state_m.savedFieldDomain             = bunch.saveFieldDomainState();
    diagnostics_m.entryRhoSnapshotDumped = false;
    bunch.clearBeamBeamWindowVisualizationTail();
    applyWindowConfig(geometry, bunch);

    if (gmsg != nullptr) {
        std::ostringstream diagnostics;
        diagnostics << std::fixed << std::setprecision(3)
                    << "Entering BeamBeam window: interaction_point_s="
                    << geometry.interactionPointS << " m, s_range=(" << geometry.beginS << ", "
                    << geometry.endS << ") m, length=" << geometry.length << " m";
        if (geometry.config.xAperture.has_value() && geometry.config.yAperture.has_value()) {
            diagnostics << ", aperture_half_width=(" << *geometry.config.xAperture << ", "
                        << *geometry.config.yAperture << ") m";
        } else {
            diagnostics << ", aperture_half_width=(not set; preserving transverse field bounds)";
        }
        diagnostics << ", witness_containers=";
        if (geometry.config.witnessContainers.empty()) {
            diagnostics << "NONE";
        } else {
            diagnostics << "(";
            for (std::size_t i = 0; i < geometry.config.witnessContainers.size(); ++i) {
                if (i != 0) {
                    diagnostics << ",";
                }
                diagnostics << geometry.config.witnessContainers[i];
            }
            diagnostics << ")";
        }
        diagnostics << ", retire_time=";
        if (geometry.config.sourceRetireTime.has_value()) {
            diagnostics << *geometry.config.sourceRetireTime << " s";
        } else {
            diagnostics << "NONE";
        }
        diagnostics << ", copy_time=";
        if (geometry.config.copyTime.has_value()) {
            diagnostics << *geometry.config.copyTime << " s";
        } else {
            diagnostics << "NONE";
        }
        diagnostics << ", rigid_source="
                    << (geometry.config.rigidSource ? "TRUE" : "FALSE");
        *gmsg << level2 << diagnostics.str() << endl;
    }
    logDiagnostics(bunch, true);
    message << level5 << "start BeamBeam-window mode" << endl;
}

void BeamBeamInteraction::applyWindowConfig(
        const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch) const {
    const bool copyModel = BEAMBEAM::copyTimeReached(bunch.getT(), geometry.config.copyTime);
    bunch.setBeamBeamWindowConfig(
            geometry.length, geometry.interactionPointS, geometry.beginS, geometry.endS, copyModel,
            geometry.config.xAperture, geometry.config.yAperture);
}

std::optional<double> BeamBeamInteraction::performWindowEntryTransition(
        const BEAMBEAM::ActualGeometry& geometry, const ippl::Vector<double, 3>& physicalRMin,
        const ippl::Vector<double, 3>& physicalRMax, PartBunch_t& bunch) {
    if (diagnostics_m.entryRhoSnapshotDumped) {
        return std::nullopt;
    }
    if (!BEAMBEAM::copyTimeReached(bunch.getT(), geometry.config.copyTime)) {
        diagnostics_m.entryRhoSnapshotDumped = true;
        return std::nullopt;
    }

    IpplTimings::startTimer(entryTransitionTimer_m);
    if (bunch.hasBinning()) {
        // Keep the interaction marker active for the binned diagnostic solve, but force the
        // physical-primary-only model on the pre-enlarged mesh. The binned solver uses this
        // marker to measure deposited charge without imposing that reduction on ordinary
        // self-field steps.
        bunch.setBeamBeamWindowConfig(
                geometry.length, geometry.interactionPointS, geometry.beginS, geometry.endS,
                /*copyModel=*/false, geometry.config.xAperture, geometry.config.yAperture);
    } else {
        // Preserve the legacy pre-enlarge solve exactly; that path already records rho.sum().
        bunch.clearBeamBeamWindowConfig();
    }
    bunch.computeSelfFields();
    if (!bunch.hasLastDepositedChargeBeforeBackground()) {
        IpplTimings::stopTimer(entryTransitionTimer_m);
        throw OpalException(
                "BeamBeamInteraction::performWindowEntryTransition",
                "Missing deposited-charge diagnostics for the pre-enlarge BeamBeam solve.");
    }

    const double referenceCharge = std::abs(bunch.getLastDepositedChargeBeforeBackground());
    dumpTransitionSnapshot("before_interaction_window_mesh_enlarge", bunch);
    applyWindowConfig(geometry, bunch);
    bunch.setPhysicalBounds(physicalRMin, physicalRMax);
    diagnostics_m.entryRhoSnapshotDumped = true;
    IpplTimings::stopTimer(entryTransitionTimer_m);
    return referenceCharge;
}

void BeamBeamInteraction::validateCopiedCharge(double referenceCharge, PartBunch_t& bunch) const {
    if (!bunch.hasLastDepositedChargeBeforeBackground()) {
        throw OpalException(
                "BeamBeamInteraction::validateCopiedCharge",
                "Missing deposited-charge diagnostics for the enlarged BeamBeam solve.");
    }

    const double enlargedCharge = std::abs(bunch.getLastDepositedChargeBeforeBackground());
    const double expectedCharge = 2.0 * referenceCharge;
    const double tolerance      = std::max(1.0e-18, 1.0e-2 * expectedCharge);
    if (std::abs(enlargedCharge - expectedCharge) > tolerance) {
        std::ostringstream message;
        message << "BeamBeam enlarged-domain charge mismatch: expected " << expectedCharge
                << " C from copied bunch, got " << enlargedCharge << " C.";
        throw OpalException("BeamBeamInteraction::validateCopiedCharge", message.str());
    }
}

void BeamBeamInteraction::dumpTransitionSnapshot(
        const std::string& snapshotKind, PartBunch_t& bunch) {
    IpplTimings::startTimer(transitionDumpTimer_m);
    auto headers = bunch.buildScalarDumpHeaders(snapshotKind);
    bunch.getFieldSolver()->dumpScalField("RHO", "collwin_vis", headers);
    IpplTimings::stopTimer(transitionDumpTimer_m);
}

void BeamBeamInteraction::leaveWindow(PartBunch_t& bunch, Inform& message) {
    state_m.state                        = BEAMBEAM::WindowState::Completed;
    diagnostics_m.entryRhoSnapshotDumped = false;

    if (state_m.geometry.has_value()) {
        const auto& geometry = *state_m.geometry;
        bunch.setBeamBeamWindowVisualizationTail(
                geometry.interactionPointS, geometry.beginS, geometry.endS,
                postWindowVisualizationSteps_m);
    }

    if (state_m.savedFieldDomain.has_value()) {
        bunch.restoreFieldDomainState(*state_m.savedFieldDomain);
        bunch.calcBeamParameters();
        state_m.savedFieldDomain.reset();
    }

    bunch.clearBeamBeamWindowConfig();
    logDiagnostics(bunch, true);
    message << level5 << "finished BeamBeam-window mode" << endl;
}

void BeamBeamInteraction::retireSourceContainer(PartBunch_t& bunch, Inform& message) {
    if (!state_m.sourceRetirementPending) {
        return;
    }

    auto source = bunch.getParticleContainer(0);
    if (!source) {
        state_m.sourceRetirementPending = false;
        state_m.sourceRetired           = true;
        return;
    }

    const size_t marked  = source->markAllParticlesInvalid();
    const size_t deleted = source->deleteInvalidParticles();
    bunch.setPcInactive(0);
    state_m.sourceRetirementPending = false;
    state_m.sourceRetired           = true;

    message << level2 << "Retired BeamBeam source container[0] at RETIRE_TIME: marked " << marked
            << ", deleted " << deleted << ", remaining " << source->getTotalNum()
            << ". Witness containers remain active." << endl;
    logDiagnostics(bunch, true);
}

void BeamBeamInteraction::computeWindowSelfFields(
        const CoordinateSystemTrafo& beamToReferenceCSTrafo, PartBunch_t& bunch, Inform& message) {
    IpplTimings::startTimer(windowTimer_m);
    PAssert(state_m.geometry.has_value());
    const auto& geometry = *state_m.geometry;

    const double bunchS                = bunch.getParticleContainer()->get_sPos();
    const double interactionPointBeamZ = geometry.interactionPointS - bunchS;

    IpplTimings::startTimer(meshSetupTimer_m);
    Vector_t<double, 3> physicalRMin(0.0), physicalRMax(0.0);
    bunch.calcBeamParameters();
    bunch.get_bounds(physicalRMin, physicalRMax);
    if (geometry.config.xAperture.has_value()
        && (physicalRMin[0] < -*geometry.config.xAperture
            || physicalRMax[0] > *geometry.config.xAperture)) {
        IpplTimings::stopTimer(meshSetupTimer_m);
        IpplTimings::stopTimer(windowTimer_m);
        std::ostringstream error;
        error << "BeamBeam APERTURE x half-width " << *geometry.config.xAperture
              << " m does not contain the bunch x extent [" << physicalRMin[0] << ", "
              << physicalRMax[0]
              << "] m in the BeamBeam frame. Increase the BeamBeam APERTURE or reduce the "
                 "transverse bunch extent before using the fixed BeamBeam mesh.";
        throw OpalException("BeamBeamInteraction::computeWindowSelfFields", error.str());
    }
    if (geometry.config.yAperture.has_value()
        && (physicalRMin[1] < -*geometry.config.yAperture
            || physicalRMax[1] > *geometry.config.yAperture)) {
        IpplTimings::stopTimer(meshSetupTimer_m);
        IpplTimings::stopTimer(windowTimer_m);
        std::ostringstream error;
        error << "BeamBeam APERTURE y half-width " << *geometry.config.yAperture
              << " m does not contain the bunch y extent [" << physicalRMin[1] << ", "
              << physicalRMax[1]
              << "] m in the BeamBeam frame. Increase the BeamBeam APERTURE or reduce the "
                 "transverse bunch extent before using the fixed BeamBeam mesh.";
        throw OpalException("BeamBeamInteraction::computeWindowSelfFields", error.str());
    }
    IpplTimings::stopTimer(meshSetupTimer_m);

    const std::optional<double> preEnlargePrimaryCharge =
            performWindowEntryTransition(geometry, physicalRMin, physicalRMax, bunch);
    applyWindowConfig(geometry, bunch);

    IpplTimings::startTimer(meshSetupTimer_m);
    bunch.enableBeamBeamWindowMesh(
            interactionPointBeamZ, geometry.length, geometry.config.xAperture,
            geometry.config.yAperture);
    IpplTimings::stopTimer(meshSetupTimer_m);

    IpplTimings::startTimer(selfFieldTimer_m);
    bunch.computeSelfFields();
    IpplTimings::stopTimer(selfFieldTimer_m);
    if (preEnlargePrimaryCharge.has_value()
        && BEAMBEAM::copyTimeReached(bunch.getT(), geometry.config.copyTime)) {
        validateCopiedCharge(*preEnlargePrimaryCharge, bunch);
    }
    if (preEnlargePrimaryCharge.has_value()) {
        dumpTransitionSnapshot("after_interaction_window_mesh_enlarge", bunch);
    }
    bunch.setPhysicalBounds(physicalRMin, physicalRMax);

    IpplTimings::startTimer(transformBackTimer_m);
    transformFieldsToReferenceFrame(beamToReferenceCSTrafo, bunch, message);
    IpplTimings::stopTimer(transformBackTimer_m);
    if (!BEAMBEAM::sourceCollectiveKickEnabled(geometry.config)) {
        auto source = bunch.getParticleContainer(0);
        if (source != nullptr) {
            source->E = 0.0;
            source->B = 0.0;
            Kokkos::fence();
        }
        message << level4
                << "BBRIGID: suppressed the BeamBeam collective kick on source container[0]; "
                   "the solved mesh field remains available to witness containers."
                << endl;
    }
    message << level5 << "Compute self fields on BeamBeam-window mesh done." << endl;
    bunch.calcBeamParameters();
    IpplTimings::stopTimer(windowTimer_m);
}

void BeamBeamInteraction::transformFieldsToReferenceFrame(
        const CoordinateSystemTrafo& beamToReferenceCSTrafo, PartBunch_t& bunch,
        Inform& message) const {
    const size_t nLocal = bunch.getParticleContainer()->getLocalNum();
    beamToReferenceCSTrafo.transformBunchTo(bunch.getParticleContainer()->R.getView(), nLocal);
    message << level5 << "Transform particle positions back to reference coordinate system done."
            << endl;
    beamToReferenceCSTrafo.rotateBunchTo(bunch.getParticleContainer()->E.getView(), nLocal);
    message << level5 << "Rotate E fields back to reference coordinate system done." << endl;
    beamToReferenceCSTrafo.rotateBunchTo(bunch.getParticleContainer()->B.getView(), nLocal);
    message << level5
            << "Rotate B fields back to reference coordinate system done. ComputeSelfFields done."
            << endl;
}

void BeamBeamInteraction::gatherFieldsToWitnessContainers(PartBunch_t& bunch, Inform& message) {
    IpplTimings::startTimer(witnessGatherTimer_m);
    if (state_m.state != BEAMBEAM::WindowState::Active || !state_m.geometry.has_value()) {
        IpplTimings::stopTimer(witnessGatherTimer_m);
        return;
    }

    const auto& witnessContainers = state_m.geometry->config.witnessContainers;
    if (witnessContainers.empty()) {
        IpplTimings::stopTimer(witnessGatherTimer_m);
        return;
    }
    if (!referenceToBeamCSTrafo_m.has_value()) {
        IpplTimings::stopTimer(witnessGatherTimer_m);
        throw OpalException(
                "BeamBeamInteraction::gatherFieldsToWitnessContainers",
                "BeamBeam witness containers are configured, but source-frame transforms are "
                "not available for the current step.");
    }

    auto* solver = bunch.getFieldSolver();
    if (solver == nullptr) {
        IpplTimings::stopTimer(witnessGatherTimer_m);
        throw OpalException(
                "BeamBeamInteraction::gatherFieldsToWitnessContainers",
                "BeamBeam witness containers require an active field solver.");
    }

    const size_t nContainers = bunch.getNumParticleContainers();
    const auto source        = bunch.getParticleContainer();
    const double sourceS     = source->get_sPos();
    for (const size_t containerIndex : witnessContainers) {
        if (containerIndex == 0) {
            IpplTimings::stopTimer(witnessGatherTimer_m);
            throw OpalException(
                    "BeamBeamInteraction::gatherFieldsToWitnessContainers",
                    "container[0] is the BeamBeam source and cannot be a witness container.");
        }
        if (containerIndex >= nContainers) {
            IpplTimings::stopTimer(witnessGatherTimer_m);
            throw OpalException(
                    "BeamBeamInteraction::gatherFieldsToWitnessContainers",
                    "Configured BeamBeam witness container[" + std::to_string(containerIndex)
                            + "] is out of range for the current TRACK BEAMS list.");
        }
        if (!bunch.isPcActive(containerIndex)) {
            continue;
        }
        auto container = bunch.getParticleContainer(containerIndex);
        // IPPL gather updates distributed field halos, so every MPI rank must
        // participate when the witness container is globally nonempty. Some
        // ranks may legitimately own no local witnesses after redistribution.
        if (!container || container->getTotalNum() == 0) {
            continue;
        }

        // R is stored relative to each container's independent reference particle.
        // Translate the witness-local coordinates into the source-local frame before
        // applying the source beam rotation. The longitudinal component continues to
        // use path length, which is the authoritative curvilinear coordinate; x and y
        // come from the lab-frame reference-particle displacement. Omitting these
        // transverse components gathers every offset witness container as if its
        // reference particle were on the source axis (for example, track12 has
        // RefPartR.x = sigma_x and particle-local x approximately zero).
        Vector_t<double, 3> offsetToSourceFrame =
                container->getRefPartR() - source->getRefPartR();
        offsetToSourceFrame[2] =
                BEAMBEAM::longitudinalOffsetToSourceFrame(sourceS, container->get_sPos());
        CoordinateSystemTrafo witnessToBeamCSTrafo(
                -1.0 * offsetToSourceFrame,
                referenceToBeamCSTrafo_m->getRotation());

        const size_t nLocalBeforeRedistribution = container->getLocalNum();
        witnessToBeamCSTrafo.transformBunchTo(
                container->R.getView(), nLocalBeforeRedistribution);
        Kokkos::fence();

        // Timed witnesses can be emitted after the fixed BeamBeam field layout has been
        // initialized. Emission distributes file records by rank capacity, which does not
        // establish the spatial ownership required by IPPL gather. Redistribute only after R
        // has been transformed into the source-field frame; every rank participates because
        // the globally-empty case was rejected above. ParticleContainer::update migrates all
        // registered witness attributes together and leaves the solved source mesh unchanged.
        container->update();
        container->markMomentsDirty();

        solver->gatherCurrentFieldsToContainer(bunch, *container);
        Kokkos::fence();
        const size_t nLocalAfterRedistribution = container->getLocalNum();
        witnessToBeamCSTrafo.transformBunchFrom(
                container->R.getView(), nLocalAfterRedistribution);
        witnessToBeamCSTrafo.rotateBunchFrom(
                container->E.getView(), nLocalAfterRedistribution);
        witnessToBeamCSTrafo.rotateBunchFrom(
                container->B.getView(), nLocalAfterRedistribution);
        Kokkos::fence();

        message << level4 << "Gathered BeamBeam source fields to witness container["
                << containerIndex << "] (" << container->getTotalNum() << " particles)." << endl;
    }
    IpplTimings::stopTimer(witnessGatherTimer_m);
}

void BeamBeamInteraction::logDiagnostics(PartBunch_t& bunch, bool force) {
    if (ippl::Comm->rank() != 0) {
        return;
    }
    if (!force && !state_m.geometry.has_value() && !state_m.sourceRetired) {
        return;
    }

    const size_t nContainers = bunch.getNumParticleContainers();
    size_t activeContainers  = 0;
    for (size_t containerIndex = 0; containerIndex < nContainers; ++containerIndex) {
        if (bunch.isPcActive(containerIndex)) {
            ++activeContainers;
        }
    }

    auto source                  = bunch.getParticleContainer(0);
    const bool sourceActive      = source && bunch.isPcActive(0);
    const bool interactionActive = state_m.state == BEAMBEAM::WindowState::Active;
    const bool copyActive =
            state_m.geometry.has_value()
            && BEAMBEAM::copyTimeReached(bunch.getT(), state_m.geometry->config.copyTime);
    const bool rigidSource =
            state_m.geometry.has_value() && state_m.geometry->config.rigidSource;
    const char* stateName = "Inactive";
    if (state_m.state == BEAMBEAM::WindowState::Active) {
        stateName = "Active";
    } else if (state_m.state == BEAMBEAM::WindowState::Completed) {
        stateName = "Completed";
    }

    std::ostringstream witnessStates;
    bool hasWitnessState = false;
    if (state_m.geometry.has_value() && !state_m.geometry->config.witnessContainers.empty()) {
        for (const size_t containerIndex : state_m.geometry->config.witnessContainers) {
            if (hasWitnessState) {
                witnessStates << ",";
            }
            hasWitnessState = true;

            if (containerIndex >= nContainers) {
                witnessStates << "c" << containerIndex << ":missing";
                continue;
            }

            auto container     = bunch.getParticleContainer(containerIndex);
            const size_t total = container ? container->getTotalNum() : 0;
            const bool active  = bunch.isPcActive(containerIndex);
            witnessStates << "c" << containerIndex << ":" << (active ? "active" : "inactive")
                          << ":n=" << total;
        }
    }

    std::ostringstream signature;
    signature << stateName << "|" << activeContainers << "|" << (state_m.sourceRetired ? 1 : 0)
              << "|" << (hasWitnessState ? witnessStates.str() : "NONE") << "|" << interactionActive
              << "|" << sourceActive << "|" << copyActive << "|" << state_m.sourceRetirementPending
              << "|" << state_m.sourceBunchesOverlap << "|" << rigidSource;
    if (!force && lastDiagnosticSignature_m.has_value()
        && *lastDiagnosticSignature_m == signature.str()) {
        return;
    }
    lastDiagnosticSignature_m = signature.str();

    std::ostringstream line;
    line << std::fixed << std::setprecision(3) << "BB-DIAG BB-state=" << stateName
         << " active_bunches=" << activeContainers
         << " retired_bunches=" << (state_m.sourceRetired ? 1 : 0)
         << " witness_states=" << (hasWitnessState ? witnessStates.str() : "NONE")
         << " rigid_source=" << (rigidSource ? "TRUE" : "FALSE");
    const auto appendBoolIfChanged = [&line](const char* key, bool value,
                                             std::optional<bool>& previous) {
        if (!previous.has_value() || *previous != value) {
            line << " " << key << "=" << (value ? "TRUE" : "FALSE");
            previous = value;
        }
    };
    const auto appendBoolIfChangedAfterInitialFalse = [&line](const char* key, bool value,
                                                              std::optional<bool>& previous) {
        const bool shouldPrint = value || (previous.has_value() && *previous != value);
        if (shouldPrint) {
            line << " " << key << "=" << (value ? "TRUE" : "FALSE");
        }
        previous = value;
    };
    appendBoolIfChanged("BB-active", interactionActive, lastDiagnosticActive_m);
    appendBoolIfChanged("source_active", sourceActive, lastDiagnosticSourceActive_m);
    appendBoolIfChangedAfterInitialFalse("copy_active", copyActive, lastDiagnosticCopyActive_m);
    appendBoolIfChanged(
            "source_retirement_pending", state_m.sourceRetirementPending,
            lastDiagnosticSourceRetirementPending_m);
    appendBoolIfChangedAfterInitialFalse(
            "source_bunches_overlap", state_m.sourceBunchesOverlap, lastDiagnosticSourceOverlap_m);
    std::cout << line.str() << std::endl;
}

void BeamBeamInteraction::renderWindowFrame(
        double bunchTailS, double bunchHeadS, const BEAMBEAM::ActualGeometry& geometry) {
    if (ippl::Comm->rank() != 0) {
        return;
    }

    const bool useFrozenWindowMesh = freezesFieldMesh();
    const double bunchCenterS      = 0.5 * (bunchTailS + bunchHeadS);
    const double meshBeginS        = useFrozenWindowMesh ? geometry.beginS : bunchTailS;
    const double meshEndS          = useFrozenWindowMesh ? geometry.endS : bunchHeadS;

    windowAnimation_m->render(
            bunchCenterS, meshBeginS, meshEndS, geometry.beginS, geometry.endS,
            geometry.interactionPointS, state_m.state == BEAMBEAM::WindowState::Active,
            useFrozenWindowMesh);
}
