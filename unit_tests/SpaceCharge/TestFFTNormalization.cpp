#include <gtest/gtest.h>

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include <array>
#include <cmath>
#include <cstddef>
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
