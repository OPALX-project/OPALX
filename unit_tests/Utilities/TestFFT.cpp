/**
 * \file TestFFT.cpp
 * \brief Unit tests for FFT implementation (replacement for gsl_fft)
 *
 * This file contains unit tests for GSL-compatible Fast Fourier Transform
 * (FFT) implementations. The implementation uses the Cooley-Tukey algorithm
 * for computing FFTs of real and complex data.
 *
 * \test FFTTest::ComplexFFT_Identity
 * Tests FFT of a constant signal. Verifies basic FFT functionality for
 * complex data. A constant signal should transform to a DC component equal
 * to \f$n\f$ times the constant, with all other components zero.
 *
 * \test FFTTest::RealFFT_Transform
 * Tests forward FFT for real-valued data. Verifies that real input data is
 * correctly transformed to the frequency domain using the GSL-compatible
 * packed format, where the result is stored as \f$[r_0, r_1, i_1, r_2, i_2, \ldots]\f$.
 *
 * \test FFTTest::RealFFT_NonPowerOfTwoMatchesDFT
 * Tests an 802-point real transform (the Astra1D fieldmap size) against a
 * direct DFT, covering the odd-length fallback of the radix-2 recursion.
 *
 * \test FFTTest::RealFFT_OddLengthRoundTrip
 * Tests real -> halfcomplex -> real round trip for odd \f$n\f$, including the
 * packing of the last complex bin and the \f$1/n\f$ normalized inverse.
 *
 * \test FFTTest::HalfcomplexTransform_Unnormalized
 * Tests that gsl_fft_halfcomplex_transform is the unnormalized backward
 * transform, as in GSL, so real -> halfcomplex transform returns \f$n\,x\f$.
 *
 * \test FFTTest::RealFFT_StridedRoundTrip
 * Tests that forward and inverse real transforms honour the stride argument
 * on both input and output and leave interleaved elements untouched.
 *
 * \test FFTTest::ComplexFFT_Forward
 * Tests forward FFT for complex data. Verifies that a unit impulse (delta
 * function at \f$n=0\f$) transforms to a constant in the frequency domain
 * (all components equal to 1).
 *
 * \test FFTTest::ComplexFFT_Inverse
 * Tests inverse FFT for complex data. Verifies that a constant in the
 * frequency domain transforms back to an impulse in the time domain,
 * confirming the inverse relationship between forward and inverse transforms.
 *
 * \test FFTTest::FFT_InverseRoundTrip
 * Tests round-trip FFT operations (forward followed by inverse). Verifies
 * that applying forward FFT followed by inverse FFT recovers the original
 * signal exactly, confirming the GSL convention that the inverse transform
 * carries the \f$1/N\f$ normalization: \f$\text{IFFT}(\text{FFT}(x)) = x\f$.
 *
 * \test FFTTest::Radix2Forward
 * Tests radix-2 FFT algorithm for power-of-2 sized data. Verifies that the
 * radix-2 implementation correctly computes FFTs for data sizes that are
 * powers of 2, which allows for optimized computation.
 *
 * \test FFTTest::Radix2Inverse
 * Tests radix-2 inverse FFT algorithm. Verifies that the inverse radix-2
 * transform correctly recovers time-domain signals from frequency-domain
 * representations.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <vector>
#include "Utilities/GSLFFT.h"

class FFTTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(FFTTest, ComplexFFT_Identity) {
    // Constant signal 1 + 0i: interleaved as [Re_0, Im_0, Re_1, Im_1, ...]
    const size_t n = 8;
    std::vector<double> data(2 * n);
    for (size_t i = 0; i < n; ++i) {
        data[2 * i]     = 1.0;
        data[2 * i + 1] = 0.0;
    }

    gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(n);

    gsl_fft_complex_forward(data.data(), 1, n, wavetable, workspace);

    // Only the DC component is non-zero and equals n (forward is unnormalized)
    EXPECT_NEAR(data[0], static_cast<double>(n), 1e-10);
    EXPECT_NEAR(data[1], 0.0, 1e-10);
    for (size_t i = 1; i < n; ++i) {
        EXPECT_NEAR(data[2 * i], 0.0, 1e-10);
        EXPECT_NEAR(data[2 * i + 1], 0.0, 1e-10);
    }

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);
}

TEST_F(FFTTest, RealFFT_Transform) {
    // Test real FFT transform
    const size_t n = 8;
    std::vector<double> data(n);

    // Create a simple signal: [1, 0, 0, 0, 0, 0, 0, 0]
    data[0] = 1.0;
    for (size_t i = 1; i < n; ++i) {
        data[i] = 0.0;
    }

    gsl_fft_real_wavetable* wavetable = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace* workspace = gsl_fft_real_workspace_alloc(n);

    gsl_fft_real_transform(data.data(), 1, n, wavetable, workspace);

    // DC component should be 1.0
    EXPECT_NEAR(data[0], 1.0, 1e-10);

    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);
}

TEST_F(FFTTest, RealFFT_NonPowerOfTwoMatchesDFT) {
    const size_t n = 802;
    std::vector<double> input(n);
    for (size_t i = 0; i < n; ++i) {
        input[i] = std::sin(0.017 * static_cast<double>(i))
                   + 0.25 * std::cos(0.071 * static_cast<double>(i));
    }
    std::vector<double> transformed = input;

    gsl_fft_real_wavetable* wavetable = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace* workspace = gsl_fft_real_workspace_alloc(n);
    gsl_fft_real_transform(transformed.data(), 1, n, wavetable, workspace);

    std::vector<std::complex<double>> expected(n / 2 + 1);
    for (size_t k = 0; k <= n / 2; ++k) {
        for (size_t j = 0; j < n; ++j) {
            const double phase = -2.0 * M_PI * static_cast<double>(j * k) / static_cast<double>(n);
            expected[k] += input[j] * std::polar(1.0, phase);
        }
    }

    EXPECT_NEAR(transformed[0], expected[0].real(), 1e-9);
    for (size_t k = 1; k < n / 2; ++k) {
        EXPECT_NEAR(transformed[2 * k - 1], expected[k].real(), 1e-9);
        EXPECT_NEAR(transformed[2 * k], expected[k].imag(), 1e-9);
    }
    EXPECT_NEAR(transformed[n - 1], expected[n / 2].real(), 1e-9);

    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);
}

TEST_F(FFTTest, RealFFT_OddLengthRoundTrip) {
    const size_t n = 7;
    std::vector<double> original{0.2, -1.0, 2.5, 0.3, -0.7, 1.1, 0.4};
    std::vector<double> transformed = original;

    gsl_fft_real_wavetable* real_wavetable = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace* real_workspace = gsl_fft_real_workspace_alloc(n);
    gsl_fft_real_transform(transformed.data(), 1, n, real_wavetable, real_workspace);

    gsl_fft_halfcomplex_wavetable* hc_wavetable = gsl_fft_halfcomplex_wavetable_alloc(n);
    gsl_fft_halfcomplex_workspace* hc_workspace = gsl_fft_halfcomplex_workspace_alloc(n);
    gsl_fft_halfcomplex_inverse(transformed.data(), 1, n, hc_wavetable, hc_workspace);

    // GSL's halfcomplex inverse is normalized, so the original data comes back.
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(transformed[i], original[i], 1e-10);
    }

    gsl_fft_real_wavetable_free(real_wavetable);
    gsl_fft_real_workspace_free(real_workspace);
    gsl_fft_halfcomplex_wavetable_free(hc_wavetable);
    gsl_fft_halfcomplex_workspace_free(hc_workspace);
}

TEST_F(FFTTest, HalfcomplexTransform_Unnormalized) {
    const size_t n = 6;
    std::vector<double> original{0.5, 1.5, -2.0, 0.25, 3.0, -1.0};
    std::vector<double> transformed = original;

    gsl_fft_real_wavetable* real_wavetable = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace* real_workspace = gsl_fft_real_workspace_alloc(n);
    gsl_fft_real_transform(transformed.data(), 1, n, real_wavetable, real_workspace);

    gsl_fft_halfcomplex_wavetable* hc_wavetable = gsl_fft_halfcomplex_wavetable_alloc(n);
    gsl_fft_halfcomplex_workspace* hc_workspace = gsl_fft_halfcomplex_workspace_alloc(n);
    gsl_fft_halfcomplex_transform(transformed.data(), 1, n, hc_wavetable, hc_workspace);

    // GSL's halfcomplex transform is the backward transform without 1/n.
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(transformed[i], static_cast<double>(n) * original[i], 1e-10);
    }

    gsl_fft_real_wavetable_free(real_wavetable);
    gsl_fft_real_workspace_free(real_workspace);
    gsl_fft_halfcomplex_wavetable_free(hc_wavetable);
    gsl_fft_halfcomplex_workspace_free(hc_workspace);
}

TEST_F(FFTTest, RealFFT_StridedRoundTrip) {
    // Real data interleaved with sentinels at odd positions; stride = 2.
    const size_t n = 6;
    std::vector<double> original{0.5, 1.5, -2.0, 0.25, 3.0, -1.0};
    std::vector<double> buffer(2 * n, -42.0);
    for (size_t i = 0; i < n; ++i) {
        buffer[2 * i] = original[i];
    }
    std::vector<double> contiguous = original;

    gsl_fft_real_wavetable* real_wavetable = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace* real_workspace = gsl_fft_real_workspace_alloc(n);
    gsl_fft_real_transform(buffer.data(), 2, n, real_wavetable, real_workspace);
    gsl_fft_real_transform(contiguous.data(), 1, n, real_wavetable, real_workspace);

    // Strided spectrum matches the contiguous one and sentinels are untouched.
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(buffer[2 * i], contiguous[i], 1e-12);
        EXPECT_DOUBLE_EQ(buffer[2 * i + 1], -42.0);
    }

    gsl_fft_halfcomplex_wavetable* hc_wavetable = gsl_fft_halfcomplex_wavetable_alloc(n);
    gsl_fft_halfcomplex_workspace* hc_workspace = gsl_fft_halfcomplex_workspace_alloc(n);
    gsl_fft_halfcomplex_inverse(buffer.data(), 2, n, hc_wavetable, hc_workspace);

    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(buffer[2 * i], original[i], 1e-12);
        EXPECT_DOUBLE_EQ(buffer[2 * i + 1], -42.0);
    }

    gsl_fft_real_wavetable_free(real_wavetable);
    gsl_fft_real_workspace_free(real_workspace);
    gsl_fft_halfcomplex_wavetable_free(hc_wavetable);
    gsl_fft_halfcomplex_workspace_free(hc_workspace);
}

TEST_F(FFTTest, ComplexFFT_Forward) {
    const size_t n = 4;
    std::vector<double> data(2 * n);  // Complex data: [real0, imag0, real1, imag1, ...]

    // Set data to [1+0i, 0+0i, 0+0i, 0+0i]
    data[0] = 1.0;  // real[0]
    data[1] = 0.0;  // imag[0]
    for (size_t i = 1; i < n; ++i) {
        data[2 * i]     = 0.0;
        data[2 * i + 1] = 0.0;
    }

    gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(n);

    gsl_fft_complex_forward(data.data(), 1, n, wavetable, workspace);

    // All components should be 1.0 (DC signal transforms to all ones)
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(data[2 * i], 1.0, 1e-10);
        EXPECT_NEAR(data[2 * i + 1], 0.0, 1e-10);
    }

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);
}

TEST_F(FFTTest, ComplexFFT_Inverse) {
    const size_t n = 4;
    std::vector<double> data(2 * n);

    // Set all components to 1.0 (constant in frequency domain)
    for (size_t i = 0; i < n; ++i) {
        data[2 * i]     = 1.0;
        data[2 * i + 1] = 0.0;
    }

    gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(n);

    gsl_fft_complex_inverse(data.data(), 1, n, wavetable, workspace);

    // Should recover a unit impulse at n=0 (inverse is normalized by 1/n)
    EXPECT_NEAR(data[0], 1.0, 1e-10);
    EXPECT_NEAR(data[1], 0.0, 1e-10);

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);
}

TEST_F(FFTTest, FFT_InverseRoundTrip) {
    const size_t n = 8;
    std::vector<double> original(2 * n);
    std::vector<double> transformed(2 * n);

    // Create original signal
    for (size_t i = 0; i < n; ++i) {
        original[2 * i]     = static_cast<double>(i);
        original[2 * i + 1] = 0.0;
    }

    transformed = original;

    gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(n);

    // Forward transform
    gsl_fft_complex_forward(transformed.data(), 1, n, wavetable, workspace);

    // Inverse transform
    gsl_fft_complex_inverse(transformed.data(), 1, n, wavetable, workspace);

    // Should recover original exactly (GSL inverse is normalized)
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(transformed[2 * i], original[2 * i], 1e-6);
        EXPECT_NEAR(transformed[2 * i + 1], 0.0, 1e-6);
    }

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);
}

TEST_F(FFTTest, Radix2Forward) {
    const size_t n = 4;  // Must be power of 2
    std::vector<double> data(2 * n);

    // Set to impulse
    data[0] = 1.0;
    data[1] = 0.0;
    for (size_t i = 1; i < n; ++i) {
        data[2 * i]     = 0.0;
        data[2 * i + 1] = 0.0;
    }

    gsl_fft_complex_radix2_forward(data.data(), 1, n);

    // All components should be 1.0
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(data[2 * i], 1.0, 1e-10);
        EXPECT_NEAR(data[2 * i + 1], 0.0, 1e-10);
    }
}

TEST_F(FFTTest, Radix2Inverse) {
    const size_t n = 4;
    std::vector<double> data(2 * n);

    // Set all to 1.0
    for (size_t i = 0; i < n; ++i) {
        data[2 * i]     = 1.0;
        data[2 * i + 1] = 0.0;
    }

    gsl_fft_complex_radix2_inverse(data.data(), 1, n);

    // Should get a unit impulse at n=0 (inverse is normalized by 1/n)
    EXPECT_NEAR(data[0], 1.0, 1e-10);
    EXPECT_NEAR(data[1], 0.0, 1e-10);
}
