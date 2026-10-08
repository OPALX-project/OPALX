/**
 * @file FFT2D5Algorithm.cpp
 * @brief Implements the independent reference-path FFT2D5 space-charge algorithm.
 */

#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"

#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include <Kokkos_NumericTraits.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

namespace opalx::spacecharge {

    FFT2D5Algorithm::FFT2D5Algorithm(
            FFT2D5Config config, std::span<ParticleContainer* const> particles,
            std::shared_ptr<const BunchStateHandler> bunchState)
        : config_m(std::move(config)), bunchState_m(std::move(bunchState)) {
        validateSpaceChargeConfig(SpaceChargeConfig(config_m));
        particles_m.reserve(particles.size());
        for (ParticleContainer* particleContainer : particles) {
            if (particleContainer == nullptr) {
                throw OpalException(
                        "FFT2D5Algorithm::FFT2D5Algorithm",
                        "FFT2D5 cannot borrow a null particle container.");
            }
            particles_m.push_back(particleContainer);
        }
        if (particles_m.empty()) {
            throw OpalException(
                    "FFT2D5Algorithm::FFT2D5Algorithm",
                    "FFT2D5 requires at least one particle container.");
        }
        if (bunchState_m == nullptr) {
            throw OpalException(
                    "FFT2D5Algorithm::FFT2D5Algorithm", "The bunch state handler is null.");
        }
    }

    SpaceChargeSolveResult FFT2D5Algorithm::solve(const SpaceChargeSolveContext& context) {
        if (context.trackingActive().size() != particles_m.size()) {
            throw OpalException(
                    "FFT2D5Algorithm::solve", "The container activity set has the wrong size.");
        }
        if (context.stepState().mpiSize != 1 || ippl::Comm->size() != 1) {
            throw OpalException(
                    "FFT2D5Algorithm::solve",
                    "FFT2D5 currently supports only one MPI rank. Distributed fields and ORB "
                    "load balancing are not implemented for this solver.");
        }
        if (bunchState_m->fixedCartesianDomain().has_value()) {
            throw OpalException(
                    "FFT2D5Algorithm::solve", "FFT2D5 does not support a fixed Cartesian domain.");
        }
        ensureInitialized();

        for (std::size_t index = 0; index < particles_m.size(); ++index) {
            if (context.trackingActive()[index] != 0) {
                clearSelfFields(*particles_m[index]);
                enterSolveFrame(context.stepState().frames, *particles_m[index]);
            }
        }
        SpaceChargeSolveResult result;
        doRunSolver(context);
        result.backendSolves = sliceCount();
        for (std::size_t index = 0; index < particles_m.size(); ++index) {
            if (context.trackingActive()[index] == 0) {
                continue;
            }
            ParticleContainer& particles = *particles_m[index];
            leaveSolveFrame(context.stepState().frames, particles);
            particles.markMomentsDirty();
            particles.updateMoments();
        }
        return result;
    }

    std::size_t FFT2D5Algorithm::sliceCount() const {
        return fieldStorage_m == nullptr ? 0 : fieldStorage_m->slices().size();
    }

    const ReferencePath::View& FFT2D5Algorithm::referencePathView() const {
        if (referencePath_m == nullptr) {
            throw OpalException(
                    "FFT2D5Algorithm::referencePathView", "FFT2D5 has not been initialized.");
        }
        return referencePath_m->deviceView();
    }

    void FFT2D5Algorithm::ensureInitialized() {
        if (fieldStorage_m != nullptr) {
            return;
        }

        // Publish dependent storage together only after path loading and allocation succeed.
        auto referencePath =
                std::make_unique<ReferencePath>(ReferencePath::load(config_m.referencePathFile));
        auto fieldStorage = std::make_unique<FFT2D5FieldStorage>(config_m, referencePath->length());
        LineDensityView_t lineDensity(
                "FFT2D5LineDensity", fieldStorage->slices().size() + LineDensityGhostCells);
        LineDensityView_t lineDensityGradient(
                "FFT2D5LineDensityGradient", fieldStorage->slices().size());

        referencePath_m       = std::move(referencePath);
        fieldStorage_m        = std::move(fieldStorage);
        lineDensity_m         = std::move(lineDensity);
        lineDensityGradient_m = std::move(lineDensityGradient);
    }

}  // namespace opalx::spacecharge
