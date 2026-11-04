/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "FFT2D5Poisson.h"

#include "OpalException.h"

FFT2D5Poisson::FFT2D5Poisson() {
    params_m.add("use_heffte_defaults", false);
    params_m.add("use_pencils", true);
    params_m.add("use_gpu_aware", true);
    params_m.add("comm", ippl::a2av);
    params_m.add("r2c_direction", 0);
    params_m.add("use_reorder", true);
}

void FFT2D5Poisson::solve(ScalarField2_t& chargeDensity, VectorField2_t& electricField) {
    initialiseViews(chargeDensity);
    placeChargeDensity(chargeDensity);
    performConvolution();
    extractPotential(chargeDensity);
    determineField(chargeDensity, electricField);
}

void FFT2D5Poisson::initialiseViews(const ScalarField2_t& chargeDensity) {
    // Have the dimensions of the charge density array changed?
    const auto mesh    = chargeDensity.get_mesh();
    const Vector2_t dr = mesh.getMeshSpacing();
    const auto nx      = static_cast<size_t>(mesh.getGridsize(X));
    const auto ny      = static_cast<size_t>(mesh.getGridsize(Y));
    const auto nx2     = nx * 2;
    const auto ny2     = ny * 2;
    if (nx2_ != nx2 || ny2_ != ny2 || dr_[X] != dr[X] || dr_[Y] != dr[Y]) {
        const bool resizeFields = nx2_ != 0 && (nx2_ != nx2 || ny2_ != ny2);
        nx2_                    = nx2;
        ny2_                    = ny2;
        dr_                     = dr;
        // Create the doubled mesh
        doubledDomain_m = {nx2, ny2};
        doubledMesh_m.initialize(doubledDomain_m, {dr[X], dr[Y]}, {0, 0});
        doubledLayout_m.initialize(doubledDomain_m, {false, false});
        // Create the complex mesh
        // The x dimension has only n/2+1 as the original fields are fully real.
        complexDomain_m = {nx2 / 2 + 1, ny2};
        complexMesh_m.initialize(complexDomain_m, {dr[X], dr[Y]}, {0, 0});
        complexLayout_m.initialize(complexDomain_m, {false, false});
        // Size the arrays
        doubledRho_m.initialize(doubledMesh_m, doubledLayout_m);
        doubledRhoTr_m.initialize(complexMesh_m, complexLayout_m);
        greensFn_m.initialize(doubledMesh_m, doubledLayout_m);
        greensFnTr_m.initialize(complexMesh_m, complexLayout_m);
        // Field::initialize is deliberately a no-op after first allocation.
        // Reuse its existing updateLayout API when a later solve changes the grid.
        if (resizeFields) {
            doubledRho_m.updateLayout(doubledLayout_m);
            doubledRhoTr_m.updateLayout(complexLayout_m);
            greensFn_m.updateLayout(doubledLayout_m);
            greensFnTr_m.updateLayout(complexLayout_m);
        }
        // Create the FFT
        fft_m = std::make_unique<FFT_t>(doubledLayout_m, complexLayout_m, this->params_m);
        fft_m->warmup(doubledRho_m, doubledRhoTr_m);  // Really?
        // And the Greens function
        makeGreensFn();
    }
}

void FFT2D5Poisson::makeGreensFn() {
    // Generate the Green's function
    const auto greensFn = greensFn_m.getView();
    const auto nx2      = nx2_;
    const auto ny2      = ny2_;
    const auto dr       = dr_;
    const auto nGhost   = greensFn_m.getNghost();
    Kokkos::parallel_for(
            "GenerateGreensFunction", Kokkos::MDRangePolicy({0, 0}, {nx2, ny2}),
            KOKKOS_LAMBDA(const size_t i, const size_t j) {
                // Wrap indices into physical displacement
                const size_t ix = i <= nx2 / 2 ? i : nx2 - i;
                const size_t iy = j <= ny2 / 2 ? j : ny2 - j;
                const double x  = ix * dr[0];
                const double y  = iy * dr[1];
                // Source density is piecewise constant on each rectangular cell.
                // Integrating the logarithmic singularity also makes the self
                // value scale correctly for narrow beams and anisotropic cells.
                greensFn(i + nGhost, j + nGhost) = cellAveragedGreensFunction(x, y, dr[0], dr[1]);
            });
    Kokkos::fence();
    fft_m->transform(ippl::FORWARD, greensFn_m, greensFnTr_m);
    Kokkos::fence();
}

void FFT2D5Poisson::placeChargeDensity(const ScalarField2_t& chargeDensity) {
    const auto doubledRho       = doubledRho_m.getView();
    const auto nGhostDoubledRho = doubledRho_m.getNghost();
    const auto rho              = chargeDensity.getView();
    const auto nGhostRho        = chargeDensity.getNghost();
    Kokkos::deep_copy(doubledRho, 0.0);
    Kokkos::parallel_for(
            "CopyChargeDensity", Kokkos::MDRangePolicy({0, 0}, {nx2_ / 2, ny2_ / 2}),
            KOKKOS_LAMBDA(const size_t i, const size_t j) {
                doubledRho(i + nGhostDoubledRho, j + nGhostDoubledRho) =
                        rho(i + nGhostRho, j + nGhostRho);
            });
    Kokkos::fence();
}

void FFT2D5Poisson::performConvolution() {
    fft_m->transform(ippl::FORWARD, doubledRho_m, doubledRhoTr_m);
    Kokkos::fence();
    const auto greens       = greensFnTr_m.getView();
    const auto nGhostGreens = greensFnTr_m.getNghost();
    const auto rho          = doubledRhoTr_m.getView();
    const auto nGhostRho    = doubledRhoTr_m.getNghost();
    const auto scaling      = nx2_ * ny2_ * dr_[X] * dr_[Y];
    Kokkos::parallel_for(
            "ApplyGreens",
            Kokkos::MDRangePolicy(
                    {0, 0}, {complexDomain_m[0].length(), complexDomain_m[1].length()}),
            KOKKOS_LAMBDA(const size_t i, const size_t j) {
                rho(i + nGhostRho, j + nGhostRho) *=
                        greens(i + nGhostGreens, j + nGhostGreens) * scaling;
            });
    Kokkos::fence();
    fft_m->transform(ippl::BACKWARD, doubledRho_m, doubledRhoTr_m);
    Kokkos::fence();
}

void FFT2D5Poisson::extractPotential(ScalarField2_t& potential) {
    const auto doubledPhi       = doubledRho_m.getView();
    const auto nGhostDoubledPhi = doubledRho_m.getNghost();
    const auto phi              = potential.getView();
    const auto nGhostPhi        = potential.getNghost();
    const auto nx2              = static_cast<int>(nx2_);
    const auto ny2              = static_cast<int>(ny2_);
    // The doubled convolution also contains the free-space potential just
    // outside the physical domain. Preserve it in the halo so centered
    // differences remain valid at the first and last physical cells.
    Kokkos::parallel_for(
            "CopyPotential", Kokkos::MDRangePolicy({0, 0}, {phi.extent(0), phi.extent(1)}),
            KOKKOS_LAMBDA(const int i, const int j) {
                const int ix = (i - nGhostPhi + nx2) % nx2;
                const int iy = (j - nGhostPhi + ny2) % ny2;
                phi(i, j)    = doubledPhi(ix + nGhostDoubledPhi, iy + nGhostDoubledPhi);
            });
    Kokkos::fence();
}

void FFT2D5Poisson::determineField(const ScalarField2_t& potential, VectorField2_t& electricField) {
    const auto phi       = potential.getView();
    const auto nGhostPhi = potential.getNghost();
    const auto e         = electricField.getView();
    const auto nGhostE   = electricField.getNghost();
    const auto nx        = nx2_ / 2;
    const auto ny        = ny2_ / 2;
    const auto inv2dx    = 1.0 / (2.0 * dr_[X]);
    const auto inv2dy    = 1.0 / (2.0 * dr_[Y]);
    if (nGhostPhi == 0) {
        throw OpalException(
                "FFT2D5Poisson::determineField", "Phi ghost cells must be greater than zero");
    }
    if (nx == 0 || ny == 0) {
        throw OpalException("FFT2D5Poisson::determineField", "Domain must contain physical cells");
    }
    // Now fields by central differences
    Kokkos::parallel_for(
            "CentralDifferences", Kokkos::MDRangePolicy({0, 0}, {nx, ny}),
            KOKKOS_LAMBDA(const size_t i, const size_t j) {
                e(i + nGhostE, j + nGhostE).data_m[X] = -(phi(i + 1 + nGhostPhi, j + nGhostPhi)
                                                          - phi(i - 1 + nGhostPhi, j + nGhostPhi))
                                                        * inv2dx;
                e(i + nGhostE, j + nGhostE).data_m[Y] = -(phi(i + nGhostPhi, j + 1 + nGhostPhi)
                                                          - phi(i + nGhostPhi, j - 1 + nGhostPhi))
                                                        * inv2dy;
            });
    Kokkos::fence();
}

void FFT2D5Poisson::solveBatched(
        const ScalarField3_t& density, VectorField3_t& electric, std::span<const size_t> occupied,
        int batchSize, double coupling, ScalarField3_t* potential) {
    if (batchSize <= 0 || ippl::Comm->size() != 1)
        throw OpalException(
                "FFT2D5Poisson::solveBatched", "Positive batch size and one MPI rank required.");
    const auto mesh = density.get_mesh();
    const size_t nx = mesh.getGridsize(0), ny = mesh.getGridsize(1), nz = mesh.getGridsize(2);
    const auto spacing = mesh.getMeshSpacing();
    if (nx2_ != 2 * nx || ny2_ != 2 * ny || dr_[0] != spacing[0] || dr_[1] != spacing[1])
        throw OpalException(
                "FFT2D5Poisson::solveBatched",
                "Initialize the matching transverse workspace before batching.");
    for (size_t z : occupied)
        if (z >= nz)
            throw OpalException("FFT2D5Poisson::solveBatched", "Invalid physical slice index.");
    const auto field = electric.getView();
    Kokkos::deep_copy(field, Vector3_t(0.));
    if (potential) Kokkos::deep_copy(potential->getView(), 0.);
    if (occupied.empty()) return;
    const int capacity = static_cast<int>(std::min(size_t(batchSize), nz));
    const size_t nr = nx2_ * ny2_, ncx = nx2_ / 2 + 1, nc = ncx * ny2_;
    if (!batchPlan_m || batchNx_m != nx2_ || batchNy_m != ny2_ || batchCapacity_m != capacity) {
        const heffte::box3d<long long> inbox(
                {0, 0, 0}, {static_cast<long long>(nx2_ - 1), static_cast<long long>(ny2_ - 1), 0});
        const heffte::box3d<long long> outbox(
                {0, 0, 0}, {static_cast<long long>(ncx - 1), static_cast<long long>(ny2_ - 1), 0});
        batchPlan_m = std::make_unique<BatchPlan>(
                inbox, outbox, 0, ippl::Comm->getCommunicator(),
                ippl::fft::makeHeffteOptions<BatchBackend::backend_t>(params_m));
        batchWorkspace_m = BatchWorkspace(batchPlan_m->size_workspace() * capacity);
        batchReal_m      = decltype(batchReal_m)("FFT2D5BatchReal", nr * capacity);
        batchSpectrum_m  = decltype(batchSpectrum_m)("FFT2D5BatchSpectrum", nc * capacity);
        batchGreen_m     = decltype(batchGreen_m)("FFT2D5BatchGreen", nc);
        batchNx_m        = nx2_;
        batchNy_m        = ny2_;
        batchCapacity_m  = capacity;
    }
    if (batchSlices_m.extent(0) < occupied.size())
        batchSlices_m = decltype(batchSlices_m)("FFT2D5BatchSlices", occupied.size());
    auto hostSlices = Kokkos::create_mirror_view(batchSlices_m);
    for (size_t i = 0; i < occupied.size(); ++i)
        hostSlices(i) = occupied[i];
    Kokkos::deep_copy(batchSlices_m, hostSlices);
    const auto green         = batchGreen_m;
    const auto originalGreen = greensFnTr_m.getView();
    const auto ng            = greensFnTr_m.getNghost();
    // The native outbox is x-fastest, independent of the Field view's layout.
    Kokkos::parallel_for(
            "FFT2D5Batch::green", nc,
            KOKKOS_LAMBDA(size_t p) { green(p) = originalGreen(p % ncx + ng, p / ncx + ng); });
    const auto source   = density.getView();
    const auto real     = batchReal_m;
    const auto spectrum = batchSpectrum_m;
    const auto slices   = batchSlices_m;
    const auto ghost = density.getNghost(), eg = electric.getNghost();
    const auto nx2 = nx2_, ny2 = ny2_;
    const auto scaling  = nx2_ * ny2_ * dr_[0] * dr_[1];
    const double inv2dx = 1. / (2 * dr_[0]), inv2dy = 1. / (2 * dr_[1]);
    using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
    for (size_t first = 0; first < occupied.size(); first += capacity) {
        const int count = static_cast<int>(std::min(size_t(capacity), occupied.size() - first));
        Kokkos::parallel_for(
                "FFT2D5Batch::pack",
                Policy({0, 0, 0},
                       {count, static_cast<long long>(nx2), static_cast<long long>(ny2)}),
                KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                    real(b * nr + i + nx2 * j) =
                            (i < nx && j < ny)
                                    ? source(i + ghost, j + ghost, slices(first + b) + ghost)
                                              * coupling
                                    : 0.;
                });
        Kokkos::fence();
        batchPlan_m->forward(
                count, real.data(), spectrum.data(), batchWorkspace_m.data(), heffte::scale::full);
        Kokkos::fence();
        Kokkos::parallel_for(
                "FFT2D5Batch::applyGreen", size_t(count) * nc,
                KOKKOS_LAMBDA(size_t p) { spectrum(p) *= green(p % nc) * scaling; });
        Kokkos::fence();
        batchPlan_m->backward(
                count, spectrum.data(), real.data(), batchWorkspace_m.data(), heffte::scale::none);
        Kokkos::fence();
        Kokkos::parallel_for(
                "FFT2D5Batch::field",
                Policy({0, 0, 0}, {count, static_cast<long long>(nx), static_cast<long long>(ny)}),
                KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                    const size_t base                             = b * nr;
                    field(i + eg, j + eg, slices(first + b) + eg) = Vector3_t(
                            -(real(base + (i + 1) % nx2 + nx2 * j)
                              - real(base + (i + nx2 - 1) % nx2 + nx2 * j))
                                    * inv2dx,
                            -(real(base + i + nx2 * ((j + 1) % ny2))
                              - real(base + i + nx2 * ((j + ny2 - 1) % ny2)))
                                    * inv2dy,
                            0.);
                });
        if (potential) {
            const auto phi = potential->getView();
            const auto pg  = potential->getNghost();
            Kokkos::parallel_for(
                    "FFT2D5Batch::potential",
                    Policy({0, 0, 0}, {count, static_cast<long long>(phi.extent(0)),
                                       static_cast<long long>(phi.extent(1))}),
                    KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                        const size_t x = (static_cast<int>(i) - pg + static_cast<int>(nx2))
                                         % static_cast<int>(nx2);
                        const size_t y = (static_cast<int>(j) - pg + static_cast<int>(ny2))
                                         % static_cast<int>(ny2);
                        phi(i, j, slices(first + b) + pg) = real(b * nr + x + nx2 * y);
                    });
        }
        Kokkos::fence();
    }
}
