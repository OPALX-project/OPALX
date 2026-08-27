/**
 * @file BeamBeamInteraction.h
 * @brief Per-run collective behavior of a placed BeamBeam element.
 */

#ifndef OPALX_BeamBeamInteraction_HH
#define OPALX_BeamBeamInteraction_HH

#include "AbsBeamline/BeamBeamDefinitions.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/ElementInteraction.h"
#include "PartBunch/PartBunch.h"
#include "Utility/IpplTimings.h"

#include <memory>
#include <optional>
#include <string>

class BeamBeam;
class BeamBeamWindowAnimation;

/**
 * @brief Stateful, element-scoped execution of the BeamBeam collective interaction.
 *
 * A BeamBeam element owns immutable input and placement data. This object owns
 * the mutable window, mesh, source/witness and diagnostic state for one tracking
 * run. It is created through ElementBase::createInteraction() and is never known
 * by the tracker as a concrete type.
 */
class BeamBeamInteraction final : public ElementInteraction {
public:
    explicit BeamBeamInteraction(const BeamBeam& element);
    ~BeamBeamInteraction() override;

    ElementInteractionResult execute(
            ElementInteractionPhase phase, ElementInteractionContext& context) override;

    bool freezesFieldMesh() const noexcept override;
    bool suppressesDefaultSelfField() const noexcept override;

private:
    struct LongitudinalExtent {
        double tail = 0.0;
        double head = 0.0;
    };

    bool computeSelfFields(ElementInteractionContext& context);
    void checkInRegion(ElementInteractionContext& context);
    std::optional<BEAMBEAM::ActualGeometry> detectWindow(
            ElementInteractionContext& context, const ippl::Vector<double, 3>& rmin,
            const ippl::Vector<double, 3>& rmax);
    LongitudinalExtent computeLongitudinalExtent(
            double bunchS, const ippl::Vector<double, 3>& rmin,
            const ippl::Vector<double, 3>& rmax) const;

    void enterWindow(const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch, Inform& message);
    void leaveWindow(PartBunch_t& bunch, Inform& message);
    void applyWindowConfig(const BEAMBEAM::ActualGeometry& geometry, PartBunch_t& bunch) const;
    std::optional<double> performWindowEntryTransition(
            const BEAMBEAM::ActualGeometry& geometry, const ippl::Vector<double, 3>& physicalRMin,
            const ippl::Vector<double, 3>& physicalRMax, PartBunch_t& bunch);
    void validateCopiedCharge(double referenceCharge, PartBunch_t& bunch) const;
    void dumpTransitionSnapshot(const std::string& snapshotKind, PartBunch_t& bunch);

    void computeWindowSelfFields(
            const CoordinateSystemTrafo& beamToReferenceCSTrafo, PartBunch_t& bunch,
            Inform& message);
    void updateWindowMesh(const CoordinateSystemTrafo& referenceToBeamCSTrafo, PartBunch_t& bunch);
    void transformWitnessPositionsToSourceFrame(
            const CoordinateSystemTrafo& referenceToBeamCSTrafo, PartBunch_t& bunch,
            bool toSourceFrame) const;
    void gatherFieldsToWitnessContainers(PartBunch_t& bunch, Inform& message);
    void transformFieldsToReferenceFrame(
            const CoordinateSystemTrafo& beamToReferenceCSTrafo, PartBunch_t& bunch,
            Inform& message) const;

    void logDiagnostics(PartBunch_t& bunch, bool force = false);
    void renderWindowFrame(
            double bunchTailS, double bunchHeadS, const BEAMBEAM::ActualGeometry& geometry);

    const BeamBeam& element_m;

    IpplTimings::TimerRef windowTimer_m;
    IpplTimings::TimerRef entryTransitionTimer_m;
    IpplTimings::TimerRef meshSetupTimer_m;
    IpplTimings::TimerRef selfFieldTimer_m;
    IpplTimings::TimerRef transformBackTimer_m;
    IpplTimings::TimerRef witnessGatherTimer_m;
    IpplTimings::TimerRef transitionDumpTimer_m;

    BEAMBEAM::Runtime<PartBunch_t::SavedFieldDomainState> state_m;
    BEAMBEAM::Diagnostics diagnostics_m;
    std::optional<CoordinateSystemTrafo> referenceToBeamCSTrafo_m;
    std::optional<bool> lastDiagnosticActive_m;
    std::optional<bool> lastDiagnosticSourceActive_m;
    std::optional<bool> lastDiagnosticCopyActive_m;
    std::optional<bool> lastDiagnosticSourceOverlap_m;
    std::optional<std::string> lastDiagnosticSignature_m;
    static constexpr int postWindowVisualizationSteps_m = 4;
    std::unique_ptr<BeamBeamWindowAnimation> windowAnimation_m;
};

#endif  // OPALX_BeamBeamInteraction_HH
