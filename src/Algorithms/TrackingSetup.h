//
// TrackingSetup: what TrackRun hands to ParallelTracker.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_TRACKING_SETUP_H
#define OPALX_TRACKING_SETUP_H

#include <memory>
#include <optional>
#include <vector>

#include "Algorithms/BeamStop.h"
#include "Algorithms/ClosedOrbitInitialState.h"
#include "Algorithms/Quaternion.hpp"
#include "Algorithms/SpectralTunes.h"
#include "Algorithms/StepSizeConfig.h"
#include "Distribution/SamplingBase.hpp"
#include "OPALTypes.h"
#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "Track/RunDescription.h"

class DataSink;
class ElementBase;

namespace opalx::spacecharge {
    class SpaceChargeSolver;
}  // namespace opalx::spacecharge

/// Position within drift-kick-drift at which the particle self-field is evaluated.
enum class SpaceChargeFieldUpdate {
    MIDPOINT,  ///< Solve after the first half drift (default OPALX ordering).
    PRESTEP    ///< Solve before the first half drift (historical OPAL ordering).
};

/**
 * @brief Everything TrackRun builds for one RUN and ParallelTracker reads.
 *
 * Replaces the 16 constructor arguments and the setters ParallelTracker had. The bunch, the
 * space-charge solver and the data sink stay owned by TrackRun and are borrowed here.
 */
struct TrackingSetup {
    PartBunch_t* bunch = nullptr;  ///< Borrowed particle bunch (one container per beam).
    opalx::spacecharge::SpaceChargeSolver* spaceChargeSolver = nullptr;  ///< Borrowed solver.
    opalx::spacecharge::DirichletPlaneConfig dirichletPlane;  ///< Plane used for particle loss.
    DataSink* dataSink = nullptr;  ///< Borrowed statistics and phase-space output.

    /// The elements of the tracked line, in line order, built from the run description
    /// (opalx::run::buildElements()); ParallelTracker copies them into its OpalBeamline.
    std::vector<std::shared_ptr<ElementBase>> elements;
    /// Start position and direction of the tracked line in the global frame.
    Vector_t<double, 3> lineOrigin{0.0, 0.0, 0.0};
    Quaternion lineDirection{1.0, 0.0, 0.0, 0.0};

    /// Step-size schedule, one entry per DT segment: time step (s), stop path length (m) and
    /// number of steps.
    std::vector<double> dt;
    std::vector<double> sStop;
    std::vector<unsigned long long> maxSteps;
    double sStart = 0;  ///< Path length where tracking starts (m).

    /// Per container, the samplers that emit particles during tracking.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplers;
    /// Per container, the reference momentum [beta*gamma] set on BEAM for a file
    /// distribution; an empty entry uses the bunch mean.
    std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta;
    /// Bunch placed in the solved COF orbit frame.
    std::optional<ClosedOrbitInitialState> initialOrbit;

    /// Restart from a checkpoint: completed global steps, the stored time step and the saved
    /// position in the step-size schedule.
    bool restarting                      = false;
    unsigned long long restartGlobalStep = 0;
    double restartDt                     = 0;
    StepSizeConfig::ResumePosition restartPosition{0, 0};

    double ringPeriod = 0;  ///< One-turn path length of a RING (m), zero for a LINE.
    unsigned long long requestedTurns = 0;  ///< RUN TURNS on a RING; zero if not given.
    double kineticEnergyStop          = 0;  ///< EKINSTOP (eV); zero if not given.
    std::vector<BeamStop> beamStops;        ///< Per container: ZSTOP plane, UNTILMAXSTEPS.

    bool bareTracking = false;  ///< FIELDSOLVER TYPE = NONE.
    /// Ring boundary control may repeat a step: only for solvers that allow that.
    bool allowBoundaryControl                     = true;
    SpaceChargeFieldUpdate spaceChargeFieldUpdate = SpaceChargeFieldUpdate::MIDPOINT;

    /// OPTION values of this RUN (output frequencies, checkpoints, particle loss, ...).
    opalx::run::OptionValues options;

    /// SPECTRALTUNES: start states and settings; empty tuneInitial means ordinary tracking.
    std::vector<double> tuneInitial;
    SpectralTunes::Settings tuneSettings;
};

#endif  // OPALX_TRACKING_SETUP_H
