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

#ifndef OPALX_FFT2D5SOLVER_H
#define OPALX_FFT2D5SOLVER_H

// clang-format off
#include "Ippl.h"
#include "FFT/FFT.h"
#include <span>
#include <map>
#include <complex>
#include "FFT/Backend/Heffte.h"
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
// clang-format on

namespace opalx::spacecharge::fft2d5_detail {
    /** Register the same names/order on every rank, including ranks with no occupied planes.
     * IPPL prints timers with collectives; rank-dependent registration can deadlock printing. */
    inline void registerTimers() {
        for (const char* name :
             {"2D5 Initialize",  "2D5 Warmup",   "2D5 Plan",     "2D5 Green Init",
              "2D5 First Solve", "2D5 Solve",    "2D5 Geometry", "2D5 Migration",
              "2D5 Scatter",     "2D5 Rho Halo", "2D5 Poisson",  "2D5 Line density",
              "2D5 E Halo",      "2D5 Gather",   "2D5 Pack",     "2D5 FFT Rho",
              "2D5 Convolution", "2D5 FFT Back", "2D5 Extract",  "2D5 Field"})
            IpplTimings::getTimer(name);
    }

    /** Time completed device work, not asynchronous launch submission. No MPI barrier is
     * introduced. Timers are inclusive: parent phases contain their named subphases.
     * Stop on exceptions too, so a failed diagnostic call cannot contaminate later samples. */
    template <class Work>
    void timed(const char* name, Work&& work) {
        const auto timer = IpplTimings::getTimer(name);
        Kokkos::fence("FFT2D5 timing begin");
        IpplTimings::startTimer(timer);
        try {
            work();
        } catch (...) {
            Kokkos::fence("FFT2D5 timing exception");
            IpplTimings::stopTimer(timer);
            throw;
        }
        Kokkos::fence("FFT2D5 timing end");
        IpplTimings::stopTimer(timer);
    }
}  // namespace opalx::spacecharge::fft2d5_detail

class FFT2D5Poisson {
public:
    using Vector2_t       = ippl::Vector<double, 2>;
    using Mesh2_t         = ippl::UniformCartesian<double, 2>;
    using ScalarField2_t  = ippl::Field<double, 2, Mesh2_t, Cell>;
    using VectorField2_t  = ippl::Field<Vector2_t, 2, Mesh2_t, Cell>;
    using Layout2_t       = ippl::FieldLayout<2>;
    using Vector3_t       = ippl::Vector<double, 3>;
    using FFT_t           = ippl::FFT<ippl::RCTransform, ScalarField2_t>;
    using NDIndex2_t      = ippl::NDIndex<2>;
    using ComplexField2_t = FFT_t::ComplexField;

    // API
    FFT2D5Poisson();
    void solve(ScalarField2_t& chargeDensity, VectorField2_t& electricField);

    /**
     * @brief Free-space potential of a uniformly charged rectangular source cell.
     *
     * Returns the cell average of -log(r)/(2*pi), including its finite value at
     * the origin. The convolution separately multiplies by the source-cell area.
     * Positive hx and hy are required. The arbitrary reference length in log(r)
     * changes only the potential gauge, never the electric field.
     */
    KOKKOS_INLINE_FUNCTION static double cellAveragedGreensFunction(
            const double x, const double y, const double hx, const double hy) {
        // Scale first so physical units do not enter cancellation of the four
        // rectangle-corner primitives. Restore the logarithmic gauge afterward.
        const double scale = hx > hy ? hx : hy;
        const double u     = x / scale;
        const double v     = y / scale;
        const double a     = 0.5 * hx / scale;
        const double b     = 0.5 * hy / scale;
        double average;
        if (u * u + v * v >= 64.0) {
            // Away from the singularity, 4-point Gauss integration in each
            // direction avoids subtraction of nearly equal corner primitives.
            // At r >= 8*max(hx,hy), its log-kernel truncation error is < 1e-13.
            constexpr double nodes[4] = {
                    -0.8611363115940525752, -0.3399810435848562648, 0.3399810435848562648,
                    0.8611363115940525752};
            constexpr double weights[4] = {
                    0.3478548451374538574, 0.6521451548625461426, 0.6521451548625461426,
                    0.3478548451374538574};
            average = 0.0;
            for (unsigned i = 0; i < 4; ++i) {
                for (unsigned j = 0; j < 4; ++j) {
                    const double dx = u + a * nodes[i];
                    const double dy = v + b * nodes[j];
                    average += weights[i] * weights[j] * Kokkos::log(dx * dx + dy * dy) / 8.0;
                }
            }
        } else {
            average = (logPrimitive(u + a, v + b) - logPrimitive(u - a, v + b)
                       - logPrimitive(u + a, v - b) + logPrimitive(u - a, v - b))
                      / (4.0 * a * b);
        }
        return -(Kokkos::log(scale) + average) / (2.0 * Kokkos::numbers::pi_v<double>);
    }

    using Mesh3_t        = ippl::UniformCartesian<double, 3>;
    using ScalarField3_t = ippl::Field<double, 3, Mesh3_t, Cell>;
    using VectorField3_t = ippl::Field<Vector3_t, 3, Mesh3_t, Cell>;

    /**
     * @brief Solve independent transverse slices in bounded batches using existing heFFTe R2C.
     * @details Only x and y are transformed; slice indices are independent batch slots.
     * Reuses the serial solver's normalized cell-averaged Green spectrum with unscaled
     * source/inverse transforms and source-cell-area scaling. The input density is unchanged.
     * Coupling is applied while packing physical cells into zero-padded contiguous buffers.
     * One batched Green multiplication and centered-gradient kernel replace per-slice
     * staging. Negative transverse neighbours wrap into the doubled potential, preserving
     * the free-space potential halo. No longitudinal FFT, interpolation or physics change.
     * Persistent plans/buffers resize only for grid or batch-capacity changes; partial final
     * batches use their actual count. Empty/non-selected slices and all field halos are
     * cleared; the caller subsequently sets longitudinal field ghosts. Nonzero neutral
     * slices must appear in occupied. Memory is O(batchSize*nx*ny), including workspace.
     * The optional potential output includes its transverse halo for numerical verification.
     * Both host and CUDA fields retain buffers in their native memory space; heFFTe selects
     * the matching FFT backend. Packing, Green multiplication and field extraction operate
     * on all batch slots in device kernels on CUDA. Fences before/after heFFTe calls preserve
     * ordering across Kokkos and the backend stream; no host staging of density or fields is
     * performed. Only the occupied-slice index list is copied from the host once per solve.
     * Existing heFFTe executor factories batch all x rows across the selected slices in
     * one R2C operation. A device transpose makes every y column contiguous for one C2C
     * operation; inverse transforms reverse these steps. This avoids the public plan's
     * per-slice executor loop without changing IPPL or heFFTe. Actual-count plans are cached
     * for the full batch and final remainder; geometry/capacity changes invalidate them.
     * MPI ranks own disjoint z slabs and never transform z. The public MPI_COMM_SELF plan
     * is retained solely to initialize the normalized Green spectrum on each rank.
     * @pre On one rank, initialiseViews() has received a matching transverse workspace.
     * @pre Unique occupied local slice indices; complete x-y planes and matching source/output
     * grids.
     */
    void solveBatched(
            const ScalarField3_t& density, VectorField3_t& electric,
            std::span<const size_t> occupied, int batchSize, double coupling = 1.,
            ScalarField3_t* potential = nullptr);

    /// Exclude zero-source warmup FFTs from physical phase counters; plans/Green remain timed.
    void setWarmupTiming(bool value) { warmupTiming_m = value; }

    // Accessors for unit tests
    const ScalarField2_t& getGreensFn() const { return greensFn_m; }
    const ComplexField2_t& getGreensFnTr() const { return greensFnTr_m; }
    const ScalarField2_t& getDoubledRho() const { return doubledRho_m; }

    // Constants
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t Z = 2;

    // Helpers
    void initialiseViews(const ScalarField2_t& chargeDensity);
    void makeGreensFn();
    void placeChargeDensity(const ScalarField2_t& chargeDensity);
    void performConvolution();
    /// Copy the physical potential and its halo from the doubled free-space convolution.
    void extractPotential(ScalarField2_t& potential);
    /// Centered gradient; potential must include a valid one-cell halo (as extractPotential does).
    void determineField(const ScalarField2_t& potential, VectorField2_t& electricField);

private:
    template <class Work>
    void timedPhase(const char* name, Work&& work) {
        if (warmupTiming_m)
            work();
        else
            opalx::spacecharge::fft2d5_detail::timed(name, work);
    }
    bool warmupTiming_m = false;
    // Mixed primitive of log(sqrt(x*x+y*y)), extended continuously to the axes.
    // Using atan2 on absolute coordinates preserves the odd symmetry in each
    // coordinate without branch jumps across the source-cell singularity.
    KOKKOS_INLINE_FUNCTION static double logPrimitive(const double x, const double y) {
        const double a = Kokkos::abs(x);
        const double b = Kokkos::abs(y);
        if (a == 0.0 || b == 0.0) {
            return 0.0;
        }
        const double value = 0.5 * a * b * Kokkos::log(a * a + b * b) - 1.5 * a * b
                             + 0.5 * a * a * Kokkos::atan2(b, a)
                             + 0.5 * b * b * Kokkos::atan2(a, b);
        return (x < 0.0) != (y < 0.0) ? -value : value;
    }

    using MemSpace       = ScalarField2_t::view_type::memory_space;
    using BatchBackend   = ippl::fft::HeffteR2C<double, 2, MemSpace>;
    using BatchPlan      = BatchBackend::heffte_t;
    using BatchWorkspace = BatchBackend::workspace_t;
    using Complex        = Kokkos::complex<double>;
    std::unique_ptr<BatchPlan> batchPlan_m;
    BatchWorkspace batchWorkspace_m;
    using RowExecutor = heffte::one_dim_backend<BatchBackend::backend_t>::executor;
    using RealRowExecutor = heffte::one_dim_backend<BatchBackend::backend_t>::executor_r2c;
    // A full batch and its final remainder each have a native backend plan.
    std::map<int, std::unique_ptr<RealRowExecutor>> batchXRows_m;
    std::map<int, std::unique_ptr<RowExecutor>> batchYRows_m;
    Kokkos::View<Complex*, MemSpace> batchColumns_m, batchRowScratch_m;
    Kokkos::View<double*, MemSpace> batchReal_m;
    Kokkos::View<Complex*, MemSpace> batchSpectrum_m, batchGreen_m;
    Kokkos::View<size_t*, MemSpace> batchSlices_m;
    size_t batchNx_m = 0, batchNy_m = 0;
    int batchCapacity_m    = 0;
    bool localGreenReady_m = false;

    size_t nx2_{};
    size_t ny2_{};
    Vector2_t dr_{};
    NDIndex2_t doubledDomain_m;
    NDIndex2_t complexDomain_m;
    Mesh2_t doubledMesh_m;
    Mesh2_t complexMesh_m;
    Layout2_t doubledLayout_m;
    Layout2_t complexLayout_m;
    ScalarField2_t doubledRho_m;
    ScalarField2_t greensFn_m;
    ComplexField2_t doubledRhoTr_m;
    ComplexField2_t greensFnTr_m;
    std::unique_ptr<FFT_t> fft_m;
    ippl::ParameterList params_m;
};

#endif  // OPALX_FFT2D5SOLVER_H
