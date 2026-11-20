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
 * @file FFT2D5FieldStorage.h
 * @brief Owns persistent 3D staging fields and one reusable FFT2D5 slice solver.
 */

#ifndef OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H
#define OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H

#include "Manager/datatypes.h"
#include "SpaceCharge/FFT2D5/FFT2D5Poisson.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <cstddef>
#include <memory>
#include <span>

namespace opalx::spacecharge {

    /**
     * @brief Persistent fields and a reusable Poisson workspace for one FFT2D5 run.
     *
     * The 3D fields stage deposition and gathering. The default serial path shares one pair
     * of 2D IPPL fields and one FFT plan. MPI always uses local native batches; on one rank,
     * explicit --2D5BatchSize enables bounded native
     * heFFTe batches owned by the same OPALX solver; no longitudinal transform is performed.
     * The Field API cannot borrow a Kokkos subview, so the serial path copies each slice;
     * batching instead packs physical cells directly into contiguous FFT buffers.
     */
    class FFT2D5FieldStorage final {
    public:
        using Vector2      = Vector_t<double, 2>;
        using Vector3      = Vector_t<double, 3>;
        using Mesh2        = Mesh_t<2>;
        using Mesh3        = Mesh_t<3>;
        using Layout2      = FieldLayout_t<2>;
        using Layout3      = FieldLayout_t<3>;
        using ScalarField2 = Field_t<2>;
        using ScalarField3 = Field_t<3>;
        using VectorField2 = VField_t<double, 2>;
        using VectorField3 = VField_t<double, 3>;

        struct Slice final {
            std::unique_ptr<VectorField2> electricField;
            std::unique_ptr<ScalarField2> chargeDensity;
            std::unique_ptr<FFT2D5Poisson> solver;
        };

        FFT2D5FieldStorage(const FFT2D5Config& config, double pathLength);
        /// Destroy native solver fields/plans in the same library translation unit that constructs
        /// them.
        ~FFT2D5FieldStorage();
        /** Prepare plans, Green spectra and FFT execution on a zero source, without
         * touching particles or counting a physical solve. Runtime FFT-phase samples
         * exclude this work without resetting measurements from other instances. */
        void warmup();

        FFT2D5FieldStorage(const FFT2D5FieldStorage&)            = delete;
        FFT2D5FieldStorage& operator=(const FFT2D5FieldStorage&) = delete;
        FFT2D5FieldStorage(FFT2D5FieldStorage&&)                 = delete;
        FFT2D5FieldStorage& operator=(FFT2D5FieldStorage&&)      = delete;

        [[nodiscard]] const Vector3& spacing() const { return spacing_m; }
        [[nodiscard]] const Vector3& origin() const { return origin_m; }
        [[nodiscard]] const Vector3& size() const { return size_m; }
        [[nodiscard]] const std::array<std::size_t, 3>& meshSize() const { return meshSize_m; }
        /// Return the path-coordinate deposition mesh shared by all local z planes.
        [[nodiscard]] Mesh3& mesh() { return mesh_m; }
        [[nodiscard]] Mesh2& sliceMesh() { return sliceMesh_m; }
        [[nodiscard]] Layout2& sliceLayout() { return sliceLayout_m; }
        [[nodiscard]] Layout3& layout() { return layout_m; }
        [[nodiscard]] ScalarField3& chargeDensity() { return chargeDensity_m; }
        [[nodiscard]] VectorField3& electricField() { return electricField_m; }
        /// Return the single reusable native transverse solver workspace.
        [[nodiscard]] Slice& workspace() { return workspace_m; }
        /// Return mathematical occupied-slice solves since the last reset, not FFT launch count.
        [[nodiscard]] std::size_t backendSolves() const { return backendSolves_m; }
        /// Reset the per-physical-solve occupied-slice counter.
        void resetBackendSolves() { backendSolves_m = 0; }

        /**
         * @brief Solve a physical slice using its already scaled charge density.
         * @note Native solver construction, execution, and destruction are compiled together
         * to keep CUDA kernel stubs and their device registration in the same translation unit.
         */
        void solveSlice();
        /// Effective batch capacity: zero keeps serial slices on one rank; MPI promotes zero to
        /// one.
        [[nodiscard]] int batchSize() const { return batchSize_m; }
        /// Solve the occupied physical slices together; counts mathematical slice solves.
        void solveBatchedSlices(
                std::span<const size_t> occupied, ScalarField3* potential = nullptr);

    private:
        std::array<std::size_t, 3> meshSize_m;
        Vector3 spacing_m;
        Vector3 origin_m;
        Vector3 size_m;
        ippl::NDIndex<3> domain_m;
        Mesh3 mesh_m;
        Layout3 layout_m;
        ScalarField3 chargeDensity_m;
        VectorField3 electricField_m;
        ippl::NDIndex<2> sliceDomain_m;
        Mesh2 sliceMesh_m;
        Layout2 sliceLayout_m;
        Slice workspace_m;
        std::size_t backendSolves_m = 0;
        int batchSize_m             = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H
