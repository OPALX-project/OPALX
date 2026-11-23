/**
 * @file SpaceChargeSolver.cpp
 * @brief The space-charge solve of one tracker step.
 */

#include "SpaceCharge/SpaceChargeSolver.h"

#include "Moments/BeamMoments.h"
#include "PartBunch/PartBunch.h"

#include <utility>

namespace opalx::spacecharge {

    using opalx::moments::Moment;

    void SpaceChargeSolver::solveStep(
            const PartBunch_t& bunch, int minBinEmitted,
            const std::function<EmissionProgress()>& emissionProgress) {
        Inform m("SpaceChargeSolver::solveStep");
        const size_t totalParticles = bunch.getTotalNumAllContainers();
        if (totalParticles <= static_cast<size_t>(minBinEmitted)) {
            m << level4 << "Skipping space charge until more than MINBINEMITTED=" << minBinEmitted
              << " particles are present (total=" << totalParticles << ")." << endl;
            return;
        }

        const auto frame    = stepSolveFrame(bunch, solveFrame());
        const auto emission = emissionProgress();

        for (std::size_t index = 0; index < containerActivity_m.size(); ++index) {
            containerActivity_m[index] = bunch.isPcActive(index) ? 1 : 0;
        }

        const std::size_t step = static_cast<std::size_t>(bunch.getGlobalTrackStep());

        // Algorithms borrow the stable containers at construction; this context supplies step
        // state.
        SpaceChargeStepState stepState{step,
                                       bunch.getT(),
                                       bunch.getdT(),
                                       emission.active,
                                       emission.fraction,
                                       ippl::Comm->size(),
                                       frame};
        SpaceChargeSolveContext context(containerActivity_m, std::move(stepState));
        solve(context);
        bunch.requireReferenceFrame("SpaceChargeSolver::solveStep");
        m << level3 << "Compute space-charge fields done." << endl;
    }

    ParticleFrame SpaceChargeSolver::stepSolveFrame(
            const PartBunch_t& bunch, SolveFrame request) const {
        if (request == SolveFrame::Reference) return ParticleFrame();
        const auto primary = bunch.getParticleContainer();
        primary->requireReferenceFrame("SpaceChargeSolver::stepSolveFrame");
        const bool centroid = request == SolveFrame::Centroid;
        const auto moments  = opalx::moments::compute(
                *primary, centroid ? Moment::MeanPosition | Moment::MeanMomentum
                                    : opalx::moments::Request(Moment::MeanMomentum));
        return makeSolveFrame(
                request, centroid ? moments.meanR() : Vector_t<double, 3>(0.0), moments.meanP(),
                primary->getRefPartP(), primary->get_sPos());
    }

}  // namespace opalx::spacecharge
