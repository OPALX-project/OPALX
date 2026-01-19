//
// FFT implementation to replace GSL FFT
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_GSL_FFT_HH
#define OPAL_GSL_FFT_HH

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/fft.html
/// \see Implementation on https://www.gnu.org/software/gsl/
// Helper functions for FFT
namespace {
    // Cooley-Tukey FFT
    void fft(std::vector<std::complex<double>>& x) {
        size_t n = x.size();
        if (n <= 1) return;
        
        // Divide
        std::vector<std::complex<double>> even(n / 2);
        std::vector<std::complex<double>> odd(n / 2);
        for (size_t i = 0; i < n / 2; ++i) {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }
        
        // Conquer
        fft(even);
        fft(odd);
        
        // Combine
        for (size_t i = 0; i < n / 2; ++i) {
            std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / n) * odd[i];
            x[i] = even[i] + t;
            x[i + n / 2] = even[i] - t;
        }
    }
    
    // Inverse FFT
    void ifft(std::vector<std::complex<double>>& x) {
        size_t n = x.size();
        if (n <= 1) return;
        
        // Conjugate input
        for (auto& val : x) {
            val = std::conj(val);
        }
        
        // Forward FFT
        fft(x);
        
        // Conjugate and scale
        for (auto& val : x) {
            val = std::conj(val) / static_cast<double>(n);
        }
    }
}

/// \brief Simple FFT implementation using Cooley-Tukey algorithm.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
namespace FFT {
    
    /// \brief Forward FFT for real data (packed output).
    /// \details Computes the discrete Fourier transform of a real sequence and
    /// stores it in GSL halfcomplex packed format.
    /// \param data Input/Output: real input sequence overwritten with packed spectrum.
    /// \param stride Input: element stride in \p data.
    /// \param n Input: number of samples.
    inline void fft_real_transform(double* data, size_t stride, size_t n) {
        if (n == 0) return;
        
        // Convert real data to complex
        std::vector<std::complex<double>> complex_data(n);
        for (size_t i = 0; i < n; ++i) {
            complex_data[i] = std::complex<double>(data[i * stride], 0.0);
        }
        
        // Perform FFT
        fft(complex_data);
        
        // Convert back to real (packed format like GSL)
        // GSL packs the result: [r0, r1, i1, r2, i2, ..., r(n/2), i(n/2)]
        data[0] = complex_data[0].real();
        if (n > 1) {
            size_t half = n / 2;
            for (size_t i = 1; i < half; ++i) {
                data[2 * i - 1] = complex_data[i].real();
                data[2 * i] = complex_data[i].imag();
            }
            if (n % 2 == 0) {
                data[n - 1] = complex_data[half].real();
            }
        }
    }
}

/// \brief GSL-compatible interface for FFT routines.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_real_wavetable {
    size_t n;
};

/// \brief Workspace for real FFT routines.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_real_workspace {
    size_t n;
    std::vector<double> work;
};

/// \brief Allocate a real FFT wavetable of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: wavetable pointer.
inline gsl_fft_real_wavetable* gsl_fft_real_wavetable_alloc(size_t n) {
    gsl_fft_real_wavetable* w = new gsl_fft_real_wavetable;
    w->n = n;
    return w;
}

/// \brief Allocate a real FFT workspace of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: workspace pointer.
inline gsl_fft_real_workspace* gsl_fft_real_workspace_alloc(size_t n) {
    gsl_fft_real_workspace* w = new gsl_fft_real_workspace;
    w->n = n;
    w->work.resize(n);
    return w;
}

/// \brief Forward real FFT with GSL-compatible packed output.
/// \param data Input/Output: real input sequence overwritten with packed spectrum.
/// \param stride Input: element stride in \p data.
/// \param n Input: number of samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_real_transform(double* data, size_t stride, size_t n,
                                    gsl_fft_real_wavetable* /* wavetable */,
                                    gsl_fft_real_workspace* /* workspace */) {
    FFT::fft_real_transform(data, stride, n);
}

/// \brief Free a real FFT wavetable.
/// \param w Input: wavetable to release (can be null).
inline void gsl_fft_real_wavetable_free(gsl_fft_real_wavetable* w) {
    delete w;
}

/// \brief Free a real FFT workspace.
/// \param w Input: workspace to release (can be null).
inline void gsl_fft_real_workspace_free(gsl_fft_real_workspace* w) {
    delete w;
}

/// \brief Halfcomplex FFT types for inverse transforms.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_halfcomplex_wavetable {
    size_t n;
};

/// \brief Workspace for halfcomplex FFT routines.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_halfcomplex_workspace {
    size_t n;
    std::vector<double> work;
};

/// \brief Allocate a halfcomplex FFT wavetable of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: wavetable pointer.
inline gsl_fft_halfcomplex_wavetable* gsl_fft_halfcomplex_wavetable_alloc(size_t n) {
    gsl_fft_halfcomplex_wavetable* w = new gsl_fft_halfcomplex_wavetable;
    w->n = n;
    return w;
}

/// \brief Allocate a halfcomplex FFT workspace of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: workspace pointer.
inline gsl_fft_halfcomplex_workspace* gsl_fft_halfcomplex_workspace_alloc(size_t n) {
    gsl_fft_halfcomplex_workspace* w = new gsl_fft_halfcomplex_workspace;
    w->n = n;
    w->work.resize(n);
    return w;
}

/// \brief Inverse transform from halfcomplex packed data.
/// \param data Input/Output: packed spectrum overwritten with real sequence.
/// \param stride Input: element stride in \p data.
/// \param n Input: number of samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_halfcomplex_transform(double* data, size_t stride, size_t n,
                                          gsl_fft_halfcomplex_wavetable* /* wavetable */,
                                          gsl_fft_halfcomplex_workspace* /* workspace */) {
    // Inverse FFT for halfcomplex data (packed format)
    // Convert halfcomplex to complex, perform inverse FFT, convert back
    std::vector<std::complex<double>> complex_data(n);
    complex_data[0] = std::complex<double>(data[0], 0.0);
    
    size_t half = n / 2;
    for (size_t i = 1; i < half; ++i) {
        complex_data[i] = std::complex<double>(data[2 * i - 1], data[2 * i]);
        complex_data[n - i] = std::conj(complex_data[i]);
    }
    if (n % 2 == 0) {
        complex_data[half] = std::complex<double>(data[n - 1], 0.0);
    }
    
    // Inverse FFT (GSL convention: scale by n)
    // Our ifft scales by 1/n, so we multiply by n to get n scaling
    ifft(complex_data);
    double scale = static_cast<double>(n);
    
    // Convert back to real
    for (size_t i = 0; i < n; ++i) {
        data[i * stride] = complex_data[i].real() * scale;
    }
}

/// \brief Free a halfcomplex FFT wavetable.
/// \param w Input: wavetable to release (can be null).
inline void gsl_fft_halfcomplex_wavetable_free(gsl_fft_halfcomplex_wavetable* w) {
    delete w;
}

/// \brief Free a halfcomplex FFT workspace.
/// \param w Input: workspace to release (can be null).
inline void gsl_fft_halfcomplex_workspace_free(gsl_fft_halfcomplex_workspace* w) {
    delete w;
}

/// \brief Alias for halfcomplex inverse transform.
/// \param data Input/Output: packed spectrum overwritten with real sequence.
/// \param stride Input: element stride in \p data.
/// \param n Input: number of samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_halfcomplex_inverse(double* data, size_t stride, size_t n,
                                         gsl_fft_halfcomplex_wavetable* wavetable,
                                         gsl_fft_halfcomplex_workspace* workspace) {
    gsl_fft_halfcomplex_transform(data, stride, n, wavetable, workspace);
}

/// \brief Allow real workspace to be used with halfcomplex (compat helper).
/// \param data Input/Output: packed spectrum overwritten with real sequence.
/// \param stride Input: element stride in \p data.
/// \param n Input: number of samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_halfcomplex_inverse(double* data, size_t stride, size_t n,
                                         gsl_fft_halfcomplex_wavetable* /* wavetable */,
                                         gsl_fft_real_workspace* /* workspace */) {
    gsl_fft_halfcomplex_wavetable* hc_w = gsl_fft_halfcomplex_wavetable_alloc(n);
    gsl_fft_halfcomplex_workspace* hc_ws = gsl_fft_halfcomplex_workspace_alloc(n);
    gsl_fft_halfcomplex_transform(data, stride, n, hc_w, hc_ws);
    gsl_fft_halfcomplex_wavetable_free(hc_w);
    gsl_fft_halfcomplex_workspace_free(hc_ws);
}

/// \brief Complex FFT types.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_complex_wavetable {
    size_t n;
};

/// \brief Workspace for complex FFT routines.
/// \see https://www.gnu.org/software/gsl/doc/html/fft.html
struct gsl_fft_complex_workspace {
    size_t n;
    std::vector<std::complex<double>> work;
};

/// \brief Allocate a complex FFT wavetable of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: wavetable pointer.
inline gsl_fft_complex_wavetable* gsl_fft_complex_wavetable_alloc(size_t n) {
    gsl_fft_complex_wavetable* w = new gsl_fft_complex_wavetable();
    w->n = n;
    return w;
}

/// \brief Allocate a complex FFT workspace of size \f$n\f$.
/// \param n Input: transform size.
/// \return Output: workspace pointer.
inline gsl_fft_complex_workspace* gsl_fft_complex_workspace_alloc(size_t n) {
    gsl_fft_complex_workspace* w = new gsl_fft_complex_workspace();
    w->n = n;
    w->work.resize(n);
    return w;
}

/// \brief Forward complex FFT.
/// \details Input/output layout is interleaved: \f$[Re_0, Im_0, Re_1, Im_1, ...]\f$.
/// \param data Input/Output: interleaved complex data overwritten with spectrum.
/// \param stride Input: complex element stride.
/// \param n Input: number of complex samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_complex_forward(double* data, size_t stride, size_t n,
                                   gsl_fft_complex_wavetable* /* wavetable */,
                                   gsl_fft_complex_workspace* /* workspace */) {
    std::vector<std::complex<double>> complex_data(n);
    for (size_t i = 0; i < n; ++i) {
        complex_data[i] = std::complex<double>(data[2 * i * stride], data[(2 * i + 1) * stride]);
    }
    fft(complex_data);
    for (size_t i = 0; i < n; ++i) {
        data[2 * i * stride] = complex_data[i].real();
        data[(2 * i + 1) * stride] = complex_data[i].imag();
    }
}

/// \brief Inverse complex FFT (GSL scaling: multiplies by \f$n\f$).
/// \details Input/output layout is interleaved: \f$[Re_0, Im_0, Re_1, Im_1, ...]\f$.
/// \param data Input/Output: interleaved complex spectrum overwritten with signal.
/// \param stride Input: complex element stride.
/// \param n Input: number of complex samples.
/// \param wavetable Input: wavetable (unused).
/// \param workspace Input: workspace (unused).
inline void gsl_fft_complex_inverse(double* data, size_t stride, size_t n,
                                    gsl_fft_complex_wavetable* /* wavetable */,
                                    gsl_fft_complex_workspace* /* workspace */) {
    std::vector<std::complex<double>> complex_data(n);
    for (size_t i = 0; i < n; ++i) {
        complex_data[i] = std::complex<double>(data[2 * i * stride], data[(2 * i + 1) * stride]);
    }
    // GSL convention: inverse FFT scales by n (not 1/n like standard FFT)
    // Our ifft scales by 1/n, so we multiply by n to get n scaling
    ifft(complex_data);
    double scale = static_cast<double>(n);
    for (size_t i = 0; i < n; ++i) {
        data[2 * i * stride] = complex_data[i].real() * scale;
        data[(2 * i + 1) * stride] = complex_data[i].imag() * scale;
    }
}

/// \brief Free a complex FFT wavetable.
/// \param w Input: wavetable to release (can be null).
inline void gsl_fft_complex_wavetable_free(gsl_fft_complex_wavetable* w) {
    delete w;
}

/// \brief Free a complex FFT workspace.
/// \param w Input: workspace to release (can be null).
inline void gsl_fft_complex_workspace_free(gsl_fft_complex_workspace* w) {
    delete w;
}

/// \brief Radix-2 forward complex FFT (allocates temporary wavetable/workspace).
/// \param data Input/Output: interleaved complex data overwritten with spectrum.
/// \param stride Input: complex element stride.
/// \param n Input: number of complex samples.
inline void gsl_fft_complex_radix2_forward(double* data, size_t stride, size_t n) {
    gsl_fft_complex_wavetable* w = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* ws = gsl_fft_complex_workspace_alloc(n);
    gsl_fft_complex_forward(data, stride, n, w, ws);
    gsl_fft_complex_wavetable_free(w);
    gsl_fft_complex_workspace_free(ws);
}

/// \brief Radix-2 inverse complex FFT (GSL scaling, allocates temporaries).
/// \param data Input/Output: interleaved complex spectrum overwritten with signal.
/// \param stride Input: complex element stride.
/// \param n Input: number of complex samples.
inline void gsl_fft_complex_radix2_inverse(double* data, size_t stride, size_t n) {
    gsl_fft_complex_wavetable* w = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* ws = gsl_fft_complex_workspace_alloc(n);
    gsl_fft_complex_inverse(data, stride, n, w, ws);
    gsl_fft_complex_wavetable_free(w);
    gsl_fft_complex_workspace_free(ws);
}

#endif // OPAL_GSL_FFT_HH

