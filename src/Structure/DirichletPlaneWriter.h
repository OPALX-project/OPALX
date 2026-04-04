//
// Class DirichletPlaneWriter
//   Writes interpolated potential values on a 2D dirichlet plane to ASCII files.
//
// Copyright (c) 2026, Paul Scherrer Institut
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

#ifndef OPAL_DIRICHLET_PLANE_WRITER_H
#define OPAL_DIRICHLET_PLANE_WRITER_H

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "Utilities/OpalException.h"

/**
 * @brief Writes dirichlet-plane potential snapshots to ASCII files.
 *
 * One file is generated per dump call inside a configured output directory.
 * Each file stores metadata and one row per sampled `(x, y)` point.
 */
class DirichletPlaneWriter {
public:
    struct PlaneDiagnostics {
        double mean            = 0.0;
        double variance        = 0.0;
        std::size_t sampleCount = 0;
    };

    /// @brief Create a writer that emits files into @p outputDirectory.
    explicit DirichletPlaneWriter(const std::string& outputDirectory);

    /// @brief Write one plane snapshot.
    void writePlane(long long step,
                    double time,
                    double zPlane,
                    const std::vector<double>& xCoords,
                    const std::vector<double>& yCoords,
                    const std::vector<double>& phiValues,
                    std::size_t nx,
                    std::size_t ny,
                    const std::string& solveTag);

    /**
     * @brief Interpolate a z-plane from a 3D field on device and write it to disk.
     *
     * The interpolation follows IPPL's cell-centered convention used by scatter/gather,
     * i.e. cell center coordinates are `origin + (i + 0.5) * dx`.
     *
     * @tparam FieldType IPPL field type with `getView()`, `getLayout()`, `getNghost()`,
     *                   and mesh accessors.
     * @param step      Global tracking step.
     * @param time      Absolute simulation time (seconds).
     * @param zPlane    Physical z position to sample (meters).
     * @param field     Input scalar field to interpolate from.
     * @param solveTag  Label for output file naming.
     *
     * @return Mean/variance statistics of the sampled plane.
     */
    template <typename FieldType>
    PlaneDiagnostics dumpInterpolatedPlane(long long step,
                                           double time,
                                           double zPlane,
                                           const FieldType& field,
                                           const std::string& solveTag);

private:
    std::filesystem::path outputDirectory_m;
};

template <typename FieldType>
DirichletPlaneWriter::PlaneDiagnostics DirichletPlaneWriter::dumpInterpolatedPlane(
        long long step,
        double time,
        double zPlane,
        const FieldType& field,
        const std::string& solveTag) {
    const auto lDom    = field.getLayout().getLocalNDIndex();
    const int nghost   = field.getNghost();
    const auto spacing = field.get_mesh().getMeshSpacing();
    const auto origin  = field.get_mesh().getOrigin();

    if (spacing[0] <= 0.0 || spacing[1] <= 0.0 || spacing[2] <= 0.0) {
        throw OpalException(
                "DirichletPlaneWriter::dumpInterpolatedPlane",
                "Mesh spacing must be positive in all dimensions.");
    }

    const int nxOwned = lDom[0].length();
    const int nyOwned = lDom[1].length();
    const int nzOwned = lDom[2].length();
    if (nxOwned <= 0 || nyOwned <= 0 || nzOwned <= 0) {
        return PlaneDiagnostics{};
    }

    const std::size_t nx = static_cast<std::size_t>(nxOwned);
    const std::size_t ny = static_cast<std::size_t>(nyOwned);

    const int kMinGlobal = lDom[2].first();
    const int kMaxGlobal = lDom[2].last();

    const double kFloatRaw = (zPlane - origin[2]) / spacing[2] - 0.5;
    const double kFloat = std::min(
            std::max(kFloatRaw, static_cast<double>(kMinGlobal)),
            static_cast<double>(kMaxGlobal));

    int k0Global = static_cast<int>(std::floor(kFloat));
    int k1Global = k0Global + 1;
    k0Global     = std::max(kMinGlobal, std::min(k0Global, kMaxGlobal));
    k1Global     = std::max(kMinGlobal, std::min(k1Global, kMaxGlobal));

    const double alpha = (k1Global == k0Global) ? 0.0 : (kFloat - static_cast<double>(k0Global));

    const int k0Local = k0Global - lDom[2].first() + nghost;
    const int k1Local = k1Global - lDom[2].first() + nghost;

    using exec_space = typename FieldType::execution_space;
    using mem_space  = typename FieldType::view_type::memory_space;

    const std::size_t nxy = nx * ny;
    Kokkos::View<double**, mem_space> planeDevice("dirichlet_plane", nx, ny);
    auto fieldView = field.getView();

    Kokkos::parallel_for(
            "DirichletPlaneWriter::interpolatePlane",
            Kokkos::RangePolicy<exec_space>(0, nxy),
            KOKKOS_LAMBDA(const std::size_t idx) {
                const std::size_t ii = idx / ny;
                const std::size_t jj = idx - ii * ny;

                const int iLocal = static_cast<int>(ii) + nghost;
                const int jLocal = static_cast<int>(jj) + nghost;

                const double phi0 = fieldView(iLocal, jLocal, k0Local);
                const double phi1 = fieldView(iLocal, jLocal, k1Local);
                planeDevice(ii, jj) = (1.0 - alpha) * phi0 + alpha * phi1;
            });

    double sum = 0.0;
    Kokkos::parallel_reduce(
            "DirichletPlaneWriter::sumPlane",
            Kokkos::RangePolicy<exec_space>(0, nxy),
            KOKKOS_LAMBDA(const std::size_t idx, double& acc) {
                const std::size_t ii = idx / ny;
                const std::size_t jj = idx - ii * ny;
                acc += planeDevice(ii, jj);
            },
            sum);

    double sumSq = 0.0;
    Kokkos::parallel_reduce(
            "DirichletPlaneWriter::sumSqPlane",
            Kokkos::RangePolicy<exec_space>(0, nxy),
            KOKKOS_LAMBDA(const std::size_t idx, double& acc) {
                const std::size_t ii = idx / ny;
                const std::size_t jj = idx - ii * ny;
                const double v       = planeDevice(ii, jj);
                acc += v * v;
            },
            sumSq);

    auto planeHost = Kokkos::create_mirror_view(planeDevice);
    Kokkos::deep_copy(planeHost, planeDevice);

    std::vector<double> xCoords(nx, 0.0);
    std::vector<double> yCoords(ny, 0.0);
    std::vector<double> phiValues(nxy, 0.0);

    for (std::size_t i = 0; i < nx; ++i) {
        const int iGlobal = lDom[0].first() + static_cast<int>(i);
        xCoords[i]        = origin[0] + (static_cast<double>(iGlobal) + 0.5) * spacing[0];
    }
    for (std::size_t j = 0; j < ny; ++j) {
        const int jGlobal = lDom[1].first() + static_cast<int>(j);
        yCoords[j]        = origin[1] + (static_cast<double>(jGlobal) + 0.5) * spacing[1];
    }

    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) {
            phiValues[i * ny + j] = planeHost(i, j);
        }
    }

    writePlane(step, time, zPlane, xCoords, yCoords, phiValues, nx, ny, solveTag);

    PlaneDiagnostics diagnostics;
    diagnostics.sampleCount = nxy;
    diagnostics.mean        = sum / static_cast<double>(nxy);
    diagnostics.variance    = sumSq / static_cast<double>(nxy) - diagnostics.mean * diagnostics.mean;
    if (diagnostics.variance < 0.0) {
        diagnostics.variance = 0.0;
    }
    return diagnostics;
}

#endif  // OPAL_DIRICHLET_PLANE_WRITER_H
