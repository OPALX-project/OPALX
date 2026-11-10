//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
/**
 * @file FFT2D5FieldStorage.cpp
 * @brief Implements persistent FFT2D5 field and slice-solver construction.
 */

#include "SpaceCharge/FFT2D5/FFT2D5FieldStorage.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <array>
#include <numeric>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        ippl::NDIndex<3> makeDomain3(const std::array<std::size_t, 3>& meshSize) {
            ippl::NDIndex<3> domain;
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                domain[dimension] = ippl::Index(static_cast<int>(meshSize[dimension]));
            }
            return domain;
        }

        ippl::NDIndex<2> makeDomain2(const std::array<std::size_t, 3>& meshSize) {
            ippl::NDIndex<2> domain;
            domain[0] = ippl::Index(static_cast<int>(meshSize[0]));
            domain[1] = ippl::Index(static_cast<int>(meshSize[1]));
            return domain;
        }

        FFT2D5FieldStorage::Vector3 makeSize(const FFT2D5Config& config, double pathLength) {
            FFT2D5FieldStorage::Vector3 size;
            size[0] = config.pipeSizeX;
            size[1] = config.pipeSizeY;
            size[2] = pathLength;
            return size;
        }

        FFT2D5FieldStorage::Vector3 makeOrigin(const FFT2D5Config& config) {
            FFT2D5FieldStorage::Vector3 origin;
            origin[0] = -0.5 * config.pipeSizeX;
            origin[1] = -0.5 * config.pipeSizeY;
            origin[2] = 0.0;
            return origin;
        }

        FFT2D5FieldStorage::Vector3 makeSpacing(const FFT2D5Config& config, double pathLength) {
            FFT2D5FieldStorage::Vector3 spacing;
            const auto& meshSize = config.grid.meshSize;
            spacing[0]           = config.pipeSizeX / static_cast<double>(meshSize[0]);
            spacing[1]           = config.pipeSizeY / static_cast<double>(meshSize[1]);
            spacing[2]           = pathLength / static_cast<double>(meshSize[2]);
            return spacing;
        }

    }  // namespace

    FFT2D5FieldStorage::FFT2D5FieldStorage(const FFT2D5Config& config, double pathLength)
        : meshSize_m(config.grid.meshSize),
          spacing_m(makeSpacing(config, pathLength)),
          origin_m(makeOrigin(config)),
          size_m(makeSize(config, pathLength)),
          domain_m(makeDomain3(meshSize_m)),
          mesh_m(domain_m, spacing_m, origin_m),
          layout_m(MPI_COMM_WORLD, domain_m, config.grid.decomposition, config.closedRing),
          sliceDomain_m(makeDomain2(meshSize_m)),
          sliceMesh_m(sliceDomain_m, Vector2(spacing_m[0], spacing_m[1]), Vector2(0, 0)),
          sliceLayout_m(MPI_COMM_SELF, sliceDomain_m, std::array<bool, 2>{false, false}) {
        if (!(pathLength > 0.0)) {
            throw OpalException(
                    "FFT2D5FieldStorage::FFT2D5FieldStorage",
                    "The FFT2D5 reference-path length must be positive.");
        }

        batchSize_m = Options::fft2d5BatchSize;
        if (ippl::Comm->size() > 1 && batchSize_m == 0) batchSize_m = 1;
        if (batchSize_m > 0) {
            *ippl::Info << "FFT2D5: --2D5BatchSize " << batchSize_m
                        << " (independent transverse slices)" << endl;
        }
        chargeDensity_m.initialize(mesh_m, layout_m);
        electricField_m.initialize(mesh_m, layout_m);

        workspace_m.electricField = std::make_unique<VectorField2>(sliceMesh_m, sliceLayout_m);
        workspace_m.chargeDensity = std::make_unique<ScalarField2>(sliceMesh_m, sliceLayout_m);
        workspace_m.solver        = std::make_unique<FFT2D5Poisson>();
    }

    FFT2D5FieldStorage::~FFT2D5FieldStorage() = default;

    void FFT2D5FieldStorage::warmup() {
        Kokkos::deep_copy(chargeDensity_m.getView(), 0.);
        Kokkos::deep_copy(electricField_m.getView(), Vector3(0.));
        Kokkos::deep_copy(workspace_m.chargeDensity->getView(), 0.);
        Kokkos::deep_copy(workspace_m.electricField->getView(), Vector2(0.));
        workspace_m.solver->setWarmupTiming(true);
        try {
            if (batchSize_m > 0) {
                const auto count = std::min(
                        size_t(batchSize_m), size_t(layout_m.getLocalNDIndex()[2].length()));
                std::vector<size_t> planes(count);
                std::iota(planes.begin(), planes.end(), 0);
                // Force a complete zero-source batch despite its empty occupancy.
                solveBatchedSlices(planes);
            } else {
                solveSlice();
            }
        } catch (...) {
            workspace_m.solver->setWarmupTiming(false);
            throw;
        }
        workspace_m.solver->setWarmupTiming(false);
        Kokkos::deep_copy(chargeDensity_m.getView(), 0.);
        Kokkos::deep_copy(electricField_m.getView(), Vector3(0.));
        Kokkos::fence();
        resetBackendSolves();
    }

    void FFT2D5FieldStorage::solveBatchedSlices(
            std::span<const size_t> occupied, ScalarField3* potential) {
        if (ippl::Comm->size() == 1)
            workspace_m.solver->initialiseViews(*workspace_m.chargeDensity);
        workspace_m.solver->solveBatched(
                chargeDensity_m, electricField_m, occupied, batchSize_m, 1 / Physics::epsilon_0,
                potential);
        backendSolves_m += occupied.size();
    }

    void FFT2D5FieldStorage::solveSlice() {
        workspace_m.solver->solve(*workspace_m.chargeDensity, *workspace_m.electricField);
        ++backendSolves_m;
    }

}  // namespace opalx::spacecharge
