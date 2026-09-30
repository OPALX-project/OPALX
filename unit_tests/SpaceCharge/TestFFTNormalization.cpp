#include <gtest/gtest.h>

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "FFT/Backend/Heffte.h"

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime_api.h>
#endif

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

namespace opalx::spacecharge {
    namespace {

        using RealField = Field_t<3>;
        using FFT       = ippl::FFT<ippl::RCTransform, RealField>;
        using Spectrum  = FFT::ComplexField;
        using Complex   = Kokkos::complex<double>;
        using Mesh      = RealField::Mesh_t;
        using Layout    = ippl::FieldLayout<3>;
        using Vector    = ippl::Vector<double, 3>;

        constexpr int physicalSize       = 16;
        constexpr int paddedSize         = 2 * physicalSize;
        constexpr double paddedCellCount = paddedSize * paddedSize * paddedSize;

        // A conservative allowance for double-precision FFT roundoff on this small grid.
        // Expectations are independent amplitudes, so even compensating forward/inverse
        // scaling errors cannot hide behind a successful round trip.
        constexpr double tolerance = 1.0e-12;

        /** @brief Exact samples of cos(2*pi*8*i/32), without trigonometry or random input. */
        double cosineMode(int i) {
            constexpr std::array values{1.0, 0.0, -1.0, 0.0};
            return values[static_cast<std::size_t>(i % values.size())];
        }

        /** @brief Visit owned cells with both local view indices and global mesh indices. */
        template <typename Field, typename Function>
        void forOwnedCells(Field& field, Function function) {
            const auto owned = field.getLayout().getLocalNDIndex();
            const int ghost  = field.getNghost();
            for (int i = 0; i < owned[0].length(); ++i) {
                for (int j = 0; j < owned[1].length(); ++j) {
                    for (int k = 0; k < owned[2].length(); ++k) {
                        function(
                                i + ghost, j + ghost, k + ghost, i + owned[0].first(),
                                j + owned[1].first(), k + owned[2].first());
                    }
                }
            }
        }

        /** @brief Initialize the active execution-space field from deterministic host data. */
        template <typename Field, typename Function>
        void fillField(Field& field, Function value) {
            auto host = field.getHostMirror();
            Kokkos::deep_copy(host, typename Field::value_type(0.0));
            forOwnedCells(field, [&](int i, int j, int k, int x, int y, int z) {
                host(i, j, k) = value(x, y, z);
            });
            Kokkos::deep_copy(field.getView(), host);
        }

        double magnitude(double value) { return std::abs(value); }
        double magnitude(Complex value) { return Kokkos::abs(value); }

        /** @brief Compare all owned values, reporting the worst error without flooding CI logs. */
        template <typename Field, typename Function>
        void expectField(Field& field, Function expectedValue, double scale = 1.0) {
            auto actual = Kokkos::create_mirror(field.getView());
            Kokkos::deep_copy(actual, field.getView());
            double maximumError = 0.0;
            std::array<int, 3> worstCell{};
            forOwnedCells(field, [&](int i, int j, int k, int x, int y, int z) {
                const auto expected = expectedValue(x, y, z);
                const double error  = magnitude(scale * actual(i, j, k) - expected);
                if (!std::isfinite(error) || error > maximumError) {
                    maximumError =
                            std::isfinite(error) ? error : std::numeric_limits<double>::infinity();
                    worstCell = {x, y, z};
                }
            });
            EXPECT_LE(maximumError, tolerance)
                    << "Maximum absolute error at (" << worstCell[0] << ", " << worstCell[1] << ", "
                    << worstCell[2] << ")";
        }

        /**
         * @brief Isolate external FFT calls from asynchronous copies and report GPU errors.
         *
         * Device-wide synchronization is intentional in these diagnostic tests: it removes
         * stream ordering as a variable without changing the existing IPPL-path tests.
         * Check the launch error before another API call can obscure its origin.
         */
        ::testing::AssertionResult synchronizeDirectProbe(const char* stage) {
#if defined(KOKKOS_ENABLE_CUDA)
            const auto launch     = cudaGetLastError();
            const auto completion = cudaDeviceSynchronize();
            if (launch != cudaSuccess || completion != cudaSuccess) {
                return ::testing::AssertionFailure()
                       << stage << ": CUDA launch=" << cudaGetErrorString(launch)
                       << ", synchronization=" << cudaGetErrorString(completion);
            }
#elif defined(KOKKOS_ENABLE_HIP)
            const auto launch     = hipGetLastError();
            const auto completion = hipDeviceSynchronize();
            if (launch != hipSuccess || completion != hipSuccess) {
                return ::testing::AssertionFailure()
                       << stage << ": HIP launch=" << hipGetErrorString(launch)
                       << ", synchronization=" << hipGetErrorString(completion);
            }
#else
            Kokkos::fence(stage);
#endif
            return ::testing::AssertionSuccess();
        }

        /**
         * @brief Check HeFFTe directly, bypassing IPPL's field copies and transform wrapper.
         *
         * Use the same backend, options, default stream, persistent workspace and x-fast
         * layout as the IPPL R2C path. Only buffer initialization and synchronization differ.
         * Independent amplitudes distinguish a transform failure from a scaling failure.
         */
        void checkDirectHeffteForward(heffte::scale scaling) {
            using Backend     = FFT::heffteBackend;
            using DirectFFT   = heffte::fft3d_r2c<Backend, long long>;
            using MemorySpace = RealField::memory_space;
            const char* mode  = scaling == heffte::scale::none ? "none" : "full";
            SCOPED_TRACE(mode);
            ASSERT_TRUE(synchronizeDirectProbe("before direct HeFFTe plan"));
            const heffte::box3d<long long> inbox(
                    {0, 0, 0}, {paddedSize - 1, paddedSize - 1, paddedSize - 1});
            const heffte::box3d<long long> outbox(
                    {0, 0, 0}, {paddedSize / 2, paddedSize - 1, paddedSize - 1});
            DirectFFT plan(
                    inbox, outbox, 0, MPI_COMM_WORLD,
                    ippl::fft::makeHeffteOptions<Backend>(detail::commonFftParameters()));
            ASSERT_EQ(plan.size_inbox(), static_cast<long long>(paddedCellCount));
            ASSERT_EQ(plan.size_outbox(), (paddedSize / 2 + 1) * paddedSize * paddedSize);
            EXPECT_DOUBLE_EQ(plan.get_scale_factor(heffte::scale::full), 1.0 / paddedCellCount);
            Kokkos::View<double*, MemorySpace> input("direct_fft_input", plan.size_inbox());
            Kokkos::View<Complex*, MemorySpace> output("direct_fft_output", plan.size_outbox());
            DirectFFT::buffer_container<Complex> workspace(plan.size_workspace());
            auto hostInput = Kokkos::create_mirror_view(input);
            for (std::size_t i = 0; i < hostInput.extent(0); ++i) {
                hostInput(i) = 2.0 + cosineMode(static_cast<int>(i % paddedSize));
            }
            Kokkos::deep_copy(input, hostInput);
            ASSERT_TRUE(synchronizeDirectProbe("before direct HeFFTe forward"));
            plan.forward(input.data(), output.data(), workspace.data(), scaling);
            ASSERT_TRUE(synchronizeDirectProbe("after direct HeFFTe forward"));
            const auto actual = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
            std::cout << "Direct HeFFTe scale::" << mode
                      << ": plan full factor=" << plan.get_scale_factor(heffte::scale::full)
                      << ", DC=" << actual(0) << ", mode=" << actual(paddedSize / 4) << std::endl;

            // Compare order-one values for both modes. No relative comparison of the
            // two transforms: a shared error must not make this diagnostic pass.
            const double divisor   = scaling == heffte::scale::none ? paddedCellCount : 1.0;
            double maximumError    = 0.0;
            std::size_t worstIndex = 0;
            for (std::size_t i = 0; i < actual.extent(0); ++i) {
                const Complex expected(i == 0 ? 2.0 : (i == paddedSize / 4 ? 0.5 : 0.0), 0.0);
                const double error = magnitude(actual(i) / divisor - expected);
                if (!std::isfinite(error) || error > maximumError) {
                    maximumError =
                            std::isfinite(error) ? error : std::numeric_limits<double>::infinity();
                    worstIndex = i;
                }
            }
            EXPECT_LE(maximumError, tolerance)
                    << "Direct HeFFTe scale::" << mode << " at flat spectrum index " << worstIndex
                    << ", actual=" << actual(worstIndex) << ", comparison divisor=" << divisor;
        }

        /**
         * @brief Reproduce the FFT setup of the rotated CartesianPIC3D algorithm test.
         *
         * OPEN/HOCKNEY pads the physical 16^3 grid to 32^3, shortening x to 17 for R2C.
         * The fixed four-particle cloud spans (0.004, 0.004, 0.006) m. Domain expansion
         * adds 2% per side, the mesh uses 15 intervals, and p_z/(mc)=1 stretches z by
         * gamma=sqrt(2) in the bin rest frame. The origin has no effect on FFT scaling.
         *
         * These tests isolate FFT amplitudes and the Hockney convolution convention;
         * they do not validate the integrated Green function or complete PIC fields.
         */
        struct HockneyFFT {
            const Vector spacing{
                    0.004 * 1.04 / (physicalSize - 1), 0.004 * 1.04 / (physicalSize - 1),
                    0.006 * 1.04 * std::sqrt(2.0) / (physicalSize - 1)};
            const double cellVolume = spacing[0] * spacing[1] * spacing[2];
            const std::array<bool, 3> decomposition{false, false, false};
            ippl::NDIndex<3> realDomain{
                    ippl::Index(paddedSize), ippl::Index(paddedSize), ippl::Index(paddedSize)};
            ippl::NDIndex<3> complexDomain{
                    ippl::Index(paddedSize / 2 + 1), ippl::Index(paddedSize),
                    ippl::Index(paddedSize)};
            Mesh realMesh{realDomain, spacing, Vector(0.0)};
            Mesh complexMesh{complexDomain, spacing, Vector(0.0)};
            Layout realLayout{MPI_COMM_WORLD, realDomain, decomposition};
            Layout complexLayout{MPI_COMM_WORLD, complexDomain, decomposition};
            RealField real{realMesh, realLayout};
            Spectrum spectrum{complexMesh, complexLayout};
            FFT fft{realLayout, complexLayout, detail::commonFftParameters()};
        };

        class FFTNormalizationTest : public ::testing::Test {
        public:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                ASSERT_EQ(ippl::Comm->size(), 1)
                        << "This regression reproduces the one-rank CartesianPIC3D FFT layout.";
            }
        };

        TEST_F(FFTNormalizationTest, DirectHeffteForwardWithoutNormalization) {
            ASSERT_NO_FATAL_FAILURE(checkDirectHeffteForward(heffte::scale::none));
        }

        TEST_F(FFTNormalizationTest, DirectHeffteForwardWithFullNormalization) {
            ASSERT_NO_FATAL_FAILURE(checkDirectHeffteForward(heffte::scale::full));
        }

#if defined(KOKKOS_ENABLE_CUDA) && defined(Heffte_ENABLE_CUDA)
        TEST_F(FFTNormalizationTest, CudaScalingKernelAppliesNormalization) {
            ASSERT_TRUE(synchronizeDirectProbe("before standalone CUDA scaling"));
            // HeFFTe scales a complex spectrum as twice as many real entries. Exercise
            // that exact kernel size independently of FFT plans and field copies.
            // This probes scale_data itself, even when a build enables MAGMA scaling.
            constexpr long long count = 2LL * (paddedSize / 2 + 1) * paddedSize * paddedSize;
            Kokkos::View<double*, Kokkos::CudaSpace> values("cuda_fft_scale_probe", count);
            auto host = Kokkos::create_mirror_view(values);
            for (long long i = 0; i < count; ++i) {
                host(i) = paddedCellCount * (1.0 + i % 4) * (i % 2 == 0 ? 1.0 : -1.0);
            }
            Kokkos::deep_copy(values, host);
            ASSERT_TRUE(synchronizeDirectProbe("before HeFFTe CUDA scale_data launch"));
            heffte::cuda::scale_data<double, long long>(
                    nullptr, count, values.data(), 1.0 / paddedCellCount);
            ASSERT_TRUE(synchronizeDirectProbe("after HeFFTe CUDA scale_data launch"));
            Kokkos::deep_copy(host, values);
            double maximumError  = 0.0;
            long long worstIndex = 0;
            for (long long i = 0; i < count; ++i) {
                const double expected = (1.0 + i % 4) * (i % 2 == 0 ? 1.0 : -1.0);
                const double error    = std::abs(host(i) - expected);
                if (!std::isfinite(error) || error > maximumError) {
                    maximumError =
                            std::isfinite(error) ? error : std::numeric_limits<double>::infinity();
                    worstIndex = i;
                }
            }
            std::cout << "HeFFTe CUDA scale_data: first=" << host(0) << ", last=" << host(count - 1)
                      << ", max error=" << maximumError << std::endl;
            EXPECT_LE(maximumError, tolerance) << "Standalone CUDA scaling at scalar index "
                                               << worstIndex << ", actual=" << host(worstIndex);
        }
#endif

        TEST_F(FFTNormalizationTest, ForwardPreservesAbsoluteFourierAmplitudes) {
            HockneyFFT fields;
            fillField(fields.real, [](int x, int, int) {
                return 2.0 + cosineMode(x);
            });
            fields.fft.transform(ippl::FORWARD, fields.real, fields.spectrum);
            expectField(fields.spectrum, [](int x, int y, int z) {
                if (y == 0 && z == 0) {
                    if (x == 0) return Complex(2.0, 0.0);
                    if (x == paddedSize / 4) return Complex(0.5, 0.0);
                }
                return Complex(0.0, 0.0);
            });
        }

        TEST_F(FFTNormalizationTest, BackwardPreservesIndependentlySeededAmplitudes) {
            HockneyFFT fields;
            fillField(fields.real, [](int, int, int) {
                return 0.0;
            });
            fillField(fields.spectrum, [](int x, int y, int z) {
                if (y == 0 && z == 0) {
                    if (x == 0) return Complex(2.0, 0.0);
                    if (x == paddedSize / 4) return Complex(0.5, 0.0);
                }
                return Complex(0.0, 0.0);
            });
            fields.fft.transform(ippl::BACKWARD, fields.real, fields.spectrum);
            expectField(fields.real, [](int x, int, int) {
                return 2.0 + cosineMode(x);
            });
        }

        TEST_F(FFTNormalizationTest, HockneyConvolutionPreservesUnitCellCharge) {
            HockneyFFT fields;
            Spectrum kernelSpectrum(fields.complexMesh, fields.complexLayout);

            // A unit cell charge in the physical 16^3 domain is zero-padded to 32^3.
            // Density is 1/dV, matching the OPEN solver's charge-per-volume convention.
            fillField(fields.real, [&](int x, int y, int z) {
                return x == 0 && y == 0 && z == 0 ? 1.0 / fields.cellVolume : 0.0;
            });
            fields.fft.transform(ippl::FORWARD, fields.real, fields.spectrum);
            expectField(
                    fields.spectrum,
                    [](int, int, int) {
                        return Complex(1.0, 0.0);
                    },
                    paddedCellCount * fields.cellVolume);

            // An analytic kernel isolates normalization from Green-function discretization.
            fillField(fields.real, [](int x, int, int) {
                return 3.0 + 2.0 * cosineMode(x);
            });
            fields.fft.transform(ippl::FORWARD, fields.real, kernelSpectrum);
            fields.spectrum = fields.spectrum * kernelSpectrum;
            fields.fft.transform(ippl::BACKWARD, fields.real, fields.spectrum);

            // Both forwards include 1/N, while the inverse is unscaled. Hockney applies
            // N*dV, so convolution with unit cell charge must recover the original kernel.
            // Missing both forward factors would inflate this result by N^2 = 2^30.
            fields.real = fields.real * (paddedCellCount * fields.cellVolume);
            expectField(fields.real, [](int x, int, int) {
                return 3.0 + 2.0 * cosineMode(x);
            });
        }

    }  // namespace
}  // namespace opalx::spacecharge
