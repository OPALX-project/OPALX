//
// TrackingSetup: what buildTracking() hands to ParallelTracker.
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
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/Quaternion.hpp"
#include "Algorithms/StepSizeConfig.h"
#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"
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
 * @brief Everything buildTracking() makes for one RUN and ParallelTracker reads.
 *
 * Made by opalx::run::buildTracking(). The bunch, the space-charge solver and the data sink are
 * owned elsewhere (opalx::run::RunObjects, OpalData) and borrowed here.
 */
struct TrackingSetup {
    PartBunch_t* bunch = nullptr;  ///< Borrowed particle bunch (one container per beam).
    opalx::spacecharge::SpaceChargeSolver* spaceChargeSolver = nullptr;  ///< Borrowed solver.
    opalx::spacecharge::DirichletPlaneConfig dirichletPlane;  ///< Plane used for particle loss.
    DataSink* dataSink = nullptr;  ///< Borrowed statistics and phase-space output.

    /// The tracked line: the elements built from the run description, placed, with the 3D
    /// lattice files written and the element mapping built.
    OpalBeamline beamline;

    /// Step-size schedule, one entry per DT segment: time step (s), stop path length (m) and
    /// number of steps.
    std::vector<double> dt;
    std::vector<double> sStop;
    std::vector<unsigned long long> maxSteps;
    double sStart = 0;  ///< Path length where tracking starts (m).
    /// A fresh run starts this long before t = 0 (s): the largest time shift of the samplers,
    /// for flat-top pulses centred on the RF reference time. Zero on a restart.
    double globalPhaseShift = 0;

    /// Per container, the samplers that emit particles during tracking.
    std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplers;
    /// Per container, the reference momentum [beta*gamma] set on BEAM for a file
    /// distribution; an empty entry uses the bunch mean.
    std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta;

    /// Restart from a checkpoint: completed global steps, the stored time step and the saved
    /// position in the step-size schedule.
    bool restarting                      = false;
    unsigned long long restartGlobalStep = 0;
    double restartDt                     = 0;
    StepSizeConfig::ResumePosition restartPosition{0, 0};
    /// Cavity phases the checkpoint stored; autophase adds to them and checkpoints save them.
    CavityPhases cavityPhases;

    std::vector<BeamStop> beamStops;  ///< Per container: ZSTOP plane, UNTILMAXSTEPS.

    SpaceChargeFieldUpdate spaceChargeFieldUpdate = SpaceChargeFieldUpdate::MIDPOINT;

    /// OPTION values of this RUN (output frequencies, checkpoints, particle loss, ...).
    opalx::run::OptionValues options;
};

#endif  // OPALX_TRACKING_SETUP_H
