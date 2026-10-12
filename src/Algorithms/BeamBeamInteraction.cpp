/**
 * @file BeamBeamInteraction.cpp
 * @brief Runtime implementation of the BeamBeam collective element interaction.
 */

#include "Algorithms/BeamBeamInteraction.h"

#include "AbsBeamline/BeamBeam.h"
#include "Algorithms/OrbitThreader.h"
#include "SpaceCharge/BeamBeamFieldServices.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Utilities/BeamBeamWindowAnimation.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
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
            gatherFieldsToWitnessContainers(context, message);
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
    return state_m.state == BEAMBEAM::WindowState::Completed;
}

bool BeamBeamInteraction::computeSelfFields(ElementInteractionContext& context) {
    if (context.sourceOrbitThreader == nullptr || context.referenceToBeamCSTrafo == nullptr
        || context.beamToReferenceCSTrafo == nullptr || context.spaceChargeSolver == nullptr
        || context.spaceChargeContext == nullptr) {
        throw OpalException(
                "BeamBeamInteraction::computeSelfFields",
                "The self-field phase requires an orbit threader, transforms and solver.");
    }
    auto& bunch       = context.bunch;
    auto& fields      = context.spaceChargeSolver->beamBeamFields();
    const auto source = bunch.getParticleContainer();
    // Bounds are needed in solve axes, but solver calls own R/P/E/B frame conversion.
    context.referenceToBeamCSTrafo->transformBunchTo(source->R.getView(), source->getLocalNum());
    source->markMomentsDirty();
    try {
        checkInRegion(context, fields);
    } catch (...) {
        context.beamToReferenceCSTrafo->transformBunchTo(
                source->R.getView(), source->getLocalNum());
        source->markMomentsDirty();
        throw;
    }
    context.beamToReferenceCSTrafo->transformBunchTo(source->R.getView(), source->getLocalNum());
    source->markMomentsDirty();
    if (state_m.state != BEAMBEAM::WindowState::Active) {
        referenceToBeamCSTrafo_m.reset();
        // Once discovered, the placed element owns the complete self-field policy.
        return state_m.geometry.has_value();
    }
    referenceToBeamCSTrafo_m = *context.referenceToBeamCSTrafo;
    computeWindowSelfFields(context, fields);
    return true;
}

void BeamBeamInteraction::checkInRegion(ElementInteractionContext& context, FieldServices& fields) {
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

    const double fieldBeginS = BEAMBEAM::fieldWindowBegin(activeGeometry.interactionPointS);
    const double cellHalfWidth =
            0.5 * BEAMBEAM::fieldWindowLength / static_cast<double>(bunch.nr_m[2]);
    const bool leavingWindow = state_m.state == BEAMBEAM::WindowState::Active
                               && BEAMBEAM::sourceFullyExitedFieldWindow(
                                       bunchExtent.tail, activeGeometry.interactionPointS);
    Inform message("BeamBeam ");
    if (state_m.state == BEAMBEAM::WindowState::Inactive
        && BEAMBEAM::sourceFullyExitedFieldWindow(
                bunchExtent.tail, activeGeometry.interactionPointS)) {
        state_m.state = BEAMBEAM::WindowState::Completed;
        logDiagnostics(bunch, true);
        return;
    }
    if (leavingWindow) {
        leaveWindow(bunch, message, fields);
    }

    if (diagnostics_m.frameObserved) {
        renderWindowFrame(bunchExtent.tail, bunchExtent.head, activeGeometry);
    }

    if (state_m.state == BEAMBEAM::WindowState::Inactive
        && bunchExtent.tail >= fieldBeginS + cellHalfWidth) {
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

    const auto makeGeometry = [this](double windowBeginS) {
        const double windowLength = element_m.getGeometry().getElementLength();
        const double windowEndS   = windowBeginS + windowLength;
        const double interactionPointS =
                BEAMBEAM::interactionPointAtElementMidpoint(windowBeginS, windowEndS);

        std::optional<double> copyTime;
        const double copyTimeValue = element_m.getAttribute("COPY_TIME");
        if (copyTimeValue > 0.0) {
            copyTime = copyTimeValue;
        }

        return BEAMBEAM::ActualGeometry{
                interactionPointS, windowBeginS, windowEndS, windowLength,
                BEAMBEAM::Config{
                        element_m.getAttribute("VISUALIZE") != 0.0,
                        element_m.getAttribute("BBRIGID") != 0.0, copyTime,
                        BEAMBEAM::decodeWitnessContainerMask(
                                element_m.getAttribute("WITNESS_CONTAINERS_MASK"))}};
    };

    // ELEMEDGE placement is authoritative and lets the interaction suppress ordinary space
    // charge before the source reaches the numerical field window. A 6D-pose placement has no
    // path-length entrance, so retain the orbit-threader discovery fallback below.
    if (element_m.isElementPositionSet()) {
        if (element_m.getGeometry().getElementLength() <= 0.0) {
            return std::nullopt;
        }
        return makeGeometry(element_m.getElementPosition());
    }

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

        if (element_m.getGeometry().getElementLength() <= 0.0) {
            return std::nullopt;
        }

        const IndexMap::key_t range = sourceThreader.getRange(element);
        // The OrbitThreader map is only a discovery aid. Its range can be clipped by the
        // configured tracking horizon (for example, a short MAXSTEPS run that starts inside this
        // element), so using range.end as the physical window end makes the BeamBeam lifecycle
        // depend on the sampled bunch extent. ELEMEDGE and the element length are the authoritative
        // path-length geometry. For 6D-pose placement, where no ELEMEDGE path coordinate exists,
        // retain the threaded entrance as the best available path-length anchor.
        return makeGeometry(range.begin);
    }

    return std::nullopt;
}

void BeamBeamInteraction::enterWindow(
        const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch, Inform& message) {
    state_m.state            = BEAMBEAM::WindowState::Active;
    state_m.geometry         = geometry;
    const auto& domain       = bunch.cartesianDomain();
    state_m.savedFieldDomain = SavedFieldDomainState{
            domain.origin(), domain.lower(), domain.upper(), domain.spacing()};
    transverseMeshLower_m.reset();
    transverseMeshUpper_m.reset();
    diagnostics_m.entryRhoSnapshotDumped = false;

    if (gmsg != nullptr) {
        std::ostringstream diagnostics;
        diagnostics << std::fixed << std::setprecision(3)
                    << "Entering BeamBeam window: interaction_point_s="
                    << geometry.interactionPointS << " m, field_s_range=("
                    << BEAMBEAM::fieldWindowBegin(geometry.interactionPointS) << ", "
                    << BEAMBEAM::fieldWindowEnd(geometry.interactionPointS)
                    << ") m, field_length=" << BEAMBEAM::fieldWindowLength << " m";
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
        diagnostics << ", copy_time=";
        if (geometry.config.copyTime.has_value()) {
            diagnostics << *geometry.config.copyTime << " s";
        } else {
            diagnostics << "NONE";
        }
        diagnostics << ", transverse_domain=(combined non-shrinking envelope with rest-frame "
                       "cell-aspect floor)";
        diagnostics << ", rigid_source=" << (geometry.config.rigidSource ? "TRUE" : "FALSE");
        *gmsg << level2 << diagnostics.str() << endl;
    }
    logDiagnostics(bunch, true);
    message << level5 << "start BeamBeam-window mode" << endl;
}

void BeamBeamInteraction::applyWindowConfig(
        const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch, FieldServices& fields,
        bool captureChargeDensity) const {
    fields.configure(
            opalx::spacecharge::BeamBeamSolvePolicy{
                    true, BEAMBEAM::copyTimeReached(bunch.getT(), geometry.config.copyTime),
                    captureChargeDensity});
}

std::optional<double> BeamBeamInteraction::performWindowEntryTransition(
        ElementInteractionContext& context, FieldServices& fields) {
    if (diagnostics_m.entryRhoSnapshotDumped) {
        return std::nullopt;
    }
    const auto& geometry = *state_m.geometry;
    if (!BEAMBEAM::copyTimeReached(context.bunch.getT(), geometry.config.copyTime)) {
        diagnostics_m.entryRhoSnapshotDumped = true;
        return std::nullopt;
    }
    IpplTimings::startTimer(entryTransitionTimer_m);
    auto& bunch = context.bunch;
    auto source = bunch.getParticleContainer();
    // Establish a primary-only diagnostic domain in solve axes. Retaining a fixed domain
    // avoids migrating reference-frame witnesses through the source-frame mesh.
    context.referenceToBeamCSTrafo->transformBunchTo(source->R.getView(), source->getLocalNum());
    source->computeMinMaxR();
    auto lower = source->getMinR();
    auto upper = source->getMaxR();
    context.beamToReferenceCSTrafo->transformBunchTo(source->R.getView(), source->getLocalNum());
    source->markMomentsDirty();
    const double padding = fields.configuration().grid.boundingBoxIncreasePercent / 100.0;
    std::array<double, 3> low{}, high{};
    for (unsigned d = 0; d < 3; ++d) {
        const double span = std::max(upper[d] - lower[d], 1.0e-6);
        low[d]            = lower[d] - span * padding;
        high[d]           = upper[d] + span * padding;
        if (low[d] >= high[d]) {
            low[d] -= 0.5e-6;
            high[d] += 0.5e-6;
        }
    }
    bunch.getBunchStateHandler()->clearFixedCartesianDomain();
    bunch.getBunchStateHandler()->setFixedCartesianDomain(low, high);
    fields.configure(opalx::spacecharge::BeamBeamSolvePolicy{true, false, true});
    context.spaceChargeSolver->solve(*context.spaceChargeContext);
    const auto charge = fields.depositedCharge();
    if (!charge.has_value()) {
        throw OpalException(
                "BeamBeamInteraction::performWindowEntryTransition",
                "Missing measured deposited charge for the pre-enlarge solve.");
    }
    dumpTransitionSnapshot("before_interaction_window_mesh_enlarge", bunch, fields);
    diagnostics_m.entryRhoSnapshotDumped = true;
    IpplTimings::stopTimer(entryTransitionTimer_m);
    return std::abs(*charge);
}

void BeamBeamInteraction::validateCopiedCharge(
        double referenceCharge, FieldServices& fields) const {
    const auto charge = fields.depositedCharge();
    if (!charge.has_value()) {
        throw OpalException(
                "BeamBeamInteraction::validateCopiedCharge",
                "Missing measured deposited charge for the enlarged solve.");
    }
    const double enlargedCharge = std::abs(*charge);
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
        const std::string& snapshotKind, PartBunch_t& bunch, FieldServices& fields) {
    IpplTimings::startTimer(transitionDumpTimer_m);
    std::vector<std::string> headers;
    const auto add = [&headers](const std::string& key, const auto& value) {
        std::ostringstream line;
        line << std::setprecision(12) << key << "=" << value;
        headers.push_back(line.str());
    };
    const auto& domain = fields.domain();
    const auto source  = bunch.getParticleContainer();
    add("coordinate_frame", "beam_local");
    add("global_step", bunch.getGlobalTrackStep());
    add("time", bunch.getT());
    add("path_length_s", source->get_sPos());
    add("snapshot_kind", snapshotKind);
    add("interaction_window_active", 1);
    add("mesh_origin", domain.origin());
    add("mesh_spacing", domain.spacing());
    add("field_domain_rmin", domain.lower());
    add("field_domain_rmax", domain.upper());
    add("particle_total_num", bunch.getTotalNumAllContainers());
    add("particle_total_charge", source->getTotalCharge());
    if (state_m.geometry.has_value()) {
        add("interaction_point_s", state_m.geometry->interactionPointS);
        add("interaction_point_local_z", state_m.geometry->interactionPointS - source->get_sPos());
    }
    fields.dumpChargeDensity("collwin_vis", headers);
    IpplTimings::stopTimer(transitionDumpTimer_m);
}

void BeamBeamInteraction::leaveWindow(PartBunch_t& bunch, Inform& message, FieldServices& fields) {
    state_m.state                        = BEAMBEAM::WindowState::Completed;
    diagnostics_m.entryRhoSnapshotDumped = false;
    bunch.getBunchStateHandler()->clearFixedCartesianDomain();
    fields.configure(std::nullopt);
    if (state_m.savedFieldDomain.has_value()) {
        const auto& saved = *state_m.savedFieldDomain;
        bunch.cartesianDomain().setGeometry(saved.lower, saved.upper, saved.spacing, saved.origin);
        state_m.savedFieldDomain.reset();
    }
    transverseMeshLower_m.reset();
    transverseMeshUpper_m.reset();
    logDiagnostics(bunch, true);
    message << level5 << "finished BeamBeam-window mode" << endl;
}

void BeamBeamInteraction::transformWitnessPositionsToSourceFrame(
        const CoordinateSystemTrafo& referenceToBeamCSTrafo, PartBunch_t& bunch,
        bool toSourceFrame) const {
    PAssert(state_m.geometry.has_value());
    const auto source = bunch.getParticleContainer(0);
    PAssert(source != nullptr);

    const double sourceS = source->get_sPos();
    for (const size_t containerIndex : state_m.geometry->config.witnessContainers) {
        if (containerIndex == 0 || containerIndex >= bunch.getNumParticleContainers()) {
            continue;
        }
        auto container = bunch.getParticleContainer(containerIndex);
        if (!container || container->getTotalNum() == 0) {
            continue;
        }

        Vector_t<double, 3> offsetToSourceFrame = container->getRefPartR() - source->getRefPartR();
        offsetToSourceFrame[2] =
                BEAMBEAM::longitudinalOffsetToSourceFrame(sourceS, container->get_sPos());
        CoordinateSystemTrafo witnessToBeamCSTrafo(
                -1.0 * offsetToSourceFrame, referenceToBeamCSTrafo.getRotation());
        if (toSourceFrame) {
            witnessToBeamCSTrafo.transformBunchTo(container->R.getView(), container->getLocalNum());
        } else {
            witnessToBeamCSTrafo.transformBunchFrom(
                    container->R.getView(), container->getLocalNum());
        }
        container->markMomentsDirty();
    }
    Kokkos::fence();
}

void BeamBeamInteraction::updateWindowMesh(
        const CoordinateSystemTrafo& referenceToBeamCSTrafo, PartBunch_t& bunch,
        FieldServices& fields) {
    PAssert(state_m.geometry.has_value());
    const auto& geometry               = *state_m.geometry;
    const double sourceS               = bunch.getParticleContainer(0)->get_sPos();
    const double interactionPointBeamZ = geometry.interactionPointS - sourceS;

    // Particle coordinates are container-local. Express passive witnesses in the source frame
    // for distributed min/max operations and BOXINCR padding. Only configured witnesses
    // affect this envelope; they remain excluded from charge deposition.
    const auto source = bunch.getParticleContainer();
    referenceToBeamCSTrafo.transformBunchTo(source->R.getView(), source->getLocalNum());
    source->markMomentsDirty();
    transformWitnessPositionsToSourceFrame(referenceToBeamCSTrafo, bunch, true);
    try {
        Vector_t<double, 3> lower(0.0), upper(0.0);
        bool first         = true;
        const auto include = [&](const auto& container) {
            if (!container || container->getTotalNum() == 0) return;
            container->computeMinMaxR();
            const auto minR = container->getMinR();
            const auto maxR = container->getMaxR();
            for (unsigned d = 0; d < 3; ++d) {
                lower[d] = first ? minR[d] : std::min(lower[d], minR[d]);
                upper[d] = first ? maxR[d] : std::max(upper[d], maxR[d]);
            }
            first = false;
        };
        include(source);
        for (const auto index : geometry.config.witnessContainers) {
            if (index > 0 && index < bunch.getNumParticleContainers()) {
                include(bunch.getParticleContainer(index));
            }
        }
        const double padding = fields.configuration().grid.boundingBoxIncreasePercent / 100.0;
        for (unsigned d = 0; d < 3; ++d) {
            const double span = std::max(upper[d] - lower[d], 1.0e-6);
            lower[d] -= span * padding;
            upper[d] += span * padding;
        }
        const Vector_t<double, 3> aggregateLower = lower;
        const Vector_t<double, 3> aggregateUpper = upper;

        std::ostringstream containerBounds;
        const auto& containers = bunch.getParticleContainers();
        for (std::size_t index = 0; index < containers.size(); ++index) {
            const auto& container = containers[index];
            if (!container || container->getTotalNum() == 0) {
                continue;
            }
            containerBounds << " c" << index << "=(" << container->getMinR() << ", "
                            << container->getMaxR() << ")";
        }

        const auto source = bunch.getParticleContainer(0);
        PAssert(source != nullptr);
        const double gamma = Util::getGamma(source->getRefPartP());
        const double restFrameDz =
                gamma * BEAMBEAM::fieldWindowLength / static_cast<double>(bunch.nr_m[2] - 1);

        for (unsigned dimension = 0; dimension < 2; ++dimension) {
            const double minimumSpan = BEAMBEAM::minimumTransverseSpan(
                    BEAMBEAM::fieldWindowLength, bunch.nr_m[2], bunch.nr_m[dimension], gamma);
            const double span = upper[dimension] - lower[dimension];
            if (span < minimumSpan) {
                const double midpoint = 0.5 * (lower[dimension] + upper[dimension]);
                lower[dimension]      = midpoint - 0.5 * minimumSpan;
                upper[dimension]      = midpoint + 0.5 * minimumSpan;
            }
        }

        if (transverseMeshLower_m.has_value() && transverseMeshUpper_m.has_value()) {
            for (unsigned dimension = 0; dimension < 2; ++dimension) {
                lower[dimension] = std::min(lower[dimension], (*transverseMeshLower_m)[dimension]);
                upper[dimension] = std::max(upper[dimension], (*transverseMeshUpper_m)[dimension]);
            }
        }
        transverseMeshLower_m = lower;
        transverseMeshUpper_m = upper;

        Inform diagnostics("BeamBeamInteraction::updateWindowMesh");
        diagnostics << level4 << std::scientific << std::setprecision(6)
                    << "BeamBeam transverse mesh: aggregate_bounds=(" << aggregateLower << ", "
                    << aggregateUpper << "), retained_bounds=(" << lower << ", " << upper
                    << "), gamma=" << gamma << ", rest_frame_dz=" << restFrameDz << ", aspect_xy=("
                    << restFrameDz
                               / ((upper[0] - lower[0]) / static_cast<double>(bunch.nr_m[0] - 1))
                    << ", "
                    << restFrameDz
                               / ((upper[1] - lower[1]) / static_cast<double>(bunch.nr_m[1] - 1))
                    << "), container_bounds:" << containerBounds.str() << endl;

        lower[2] = interactionPointBeamZ - 0.5 * BEAMBEAM::fieldWindowLength;
        upper[2] = interactionPointBeamZ + 0.5 * BEAMBEAM::fieldWindowLength;
        bunch.getBunchStateHandler()->clearFixedCartesianDomain();
        bunch.getBunchStateHandler()->setFixedCartesianDomain(
                {lower[0], lower[1], lower[2]}, {upper[0], upper[1], upper[2]});
    } catch (...) {
        transformWitnessPositionsToSourceFrame(referenceToBeamCSTrafo, bunch, false);
        referenceToBeamCSTrafo.transformBunchFrom(source->R.getView(), source->getLocalNum());
        source->markMomentsDirty();
        throw;
    }
    transformWitnessPositionsToSourceFrame(referenceToBeamCSTrafo, bunch, false);
    referenceToBeamCSTrafo.transformBunchFrom(source->R.getView(), source->getLocalNum());
    source->markMomentsDirty();
}

void BeamBeamInteraction::computeWindowSelfFields(
        ElementInteractionContext& context, FieldServices& fields) {
    IpplTimings::startTimer(windowTimer_m);
    auto& bunch          = context.bunch;
    const auto& geometry = *state_m.geometry;
    Inform localMessage("BeamBeamInteraction::computeWindowSelfFields ");
    Inform& message            = context.message ? *context.message : localMessage;
    const auto referenceCharge = performWindowEntryTransition(context, fields);
    applyWindowConfig(geometry, bunch, fields, referenceCharge.has_value());
    IpplTimings::startTimer(meshSetupTimer_m);
    updateWindowMesh(*context.referenceToBeamCSTrafo, bunch, fields);
    IpplTimings::stopTimer(meshSetupTimer_m);
    IpplTimings::startTimer(selfFieldTimer_m);
    context.spaceChargeSolver->solve(*context.spaceChargeContext);
    IpplTimings::stopTimer(selfFieldTimer_m);
    if (referenceCharge.has_value()) {
        validateCopiedCharge(*referenceCharge, fields);
        dumpTransitionSnapshot("after_interaction_window_mesh_enlarge", bunch, fields);
    }
    if (!BEAMBEAM::sourceCollectiveKickEnabled(geometry.config)) {
        auto source = bunch.getParticleContainer();
        source->E   = 0.0;
        source->B   = 0.0;
        Kokkos::fence();
        message << level4
                << "BBRIGID: suppressed source collective kick; mesh retained for witnesses."
                << endl;
    }
    bunch.calcBeamParameters();
    IpplTimings::stopTimer(windowTimer_m);
}

void BeamBeamInteraction::gatherFieldsToWitnessContainers(
        ElementInteractionContext& context, Inform& message) {
    auto& bunch = context.bunch;
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

    if (context.spaceChargeSolver == nullptr) {
        IpplTimings::stopTimer(witnessGatherTimer_m);
        throw OpalException(
                "BeamBeamInteraction::gatherFieldsToWitnessContainers",
                "BeamBeam witness containers require an active field solver.");
    }

    auto& fields             = context.spaceChargeSolver->beamBeamFields();
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
        Vector_t<double, 3> offsetToSourceFrame = container->getRefPartR() - source->getRefPartR();
        offsetToSourceFrame[2] =
                BEAMBEAM::longitudinalOffsetToSourceFrame(sourceS, container->get_sPos());
        CoordinateSystemTrafo witnessToBeamCSTrafo(
                -1.0 * offsetToSourceFrame, referenceToBeamCSTrafo_m->getRotation());

        const size_t nLocalBeforeRedistribution = container->getLocalNum();
        witnessToBeamCSTrafo.transformBunchTo(container->R.getView(), nLocalBeforeRedistribution);
        Kokkos::fence();

        // Timed witnesses can be emitted after the fixed BeamBeam field layout has been
        // initialized. Emission distributes file records by rank capacity, which does not
        // establish the spatial ownership required by IPPL gather. Redistribute only after R
        // has been transformed into the source-field frame; every rank participates because
        // the globally-empty case was rejected above. ParticleContainer::update migrates all
        // registered witness attributes together and leaves the solved source mesh unchanged.
        container->updateLayout(fields.domain().layout(), fields.domain().mesh());
        container->update();
        container->markMomentsDirty();

        fields.gatherFields(*container);
        Kokkos::fence();
        const size_t nLocalAfterRedistribution = container->getLocalNum();
        witnessToBeamCSTrafo.transformBunchFrom(container->R.getView(), nLocalAfterRedistribution);
        witnessToBeamCSTrafo.rotateBunchFrom(container->E.getView(), nLocalAfterRedistribution);
        witnessToBeamCSTrafo.rotateBunchFrom(container->B.getView(), nLocalAfterRedistribution);
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
    if (!force && !state_m.geometry.has_value()) {
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
    const bool rigidSource = state_m.geometry.has_value() && state_m.geometry->config.rigidSource;
    const char* stateName  = "Inactive";
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
    signature << stateName << "|" << activeContainers << "|"
              << (hasWitnessState ? witnessStates.str() : "NONE") << "|" << interactionActive << "|"
              << sourceActive << "|" << copyActive << "|" << state_m.sourceBunchesOverlap << "|"
              << rigidSource;
    if (!force && lastDiagnosticSignature_m.has_value()
        && *lastDiagnosticSignature_m == signature.str()) {
        return;
    }
    lastDiagnosticSignature_m = signature.str();

    std::ostringstream line;
    line << std::fixed << std::setprecision(3) << "BB-DIAG BB-state=" << stateName
         << " active_bunches=" << activeContainers
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
    const double meshBeginS        = useFrozenWindowMesh
                                             ? BEAMBEAM::fieldWindowBegin(geometry.interactionPointS)
                                             : bunchTailS;
    const double meshEndS =
            useFrozenWindowMesh ? BEAMBEAM::fieldWindowEnd(geometry.interactionPointS) : bunchHeadS;

    windowAnimation_m->render(
            bunchCenterS, meshBeginS, meshEndS, geometry.beginS, geometry.endS,
            geometry.interactionPointS, state_m.state == BEAMBEAM::WindowState::Active,
            useFrozenWindowMesh);
}
