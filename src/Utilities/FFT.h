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

#ifndef OPAL_FFT_HH
#define OPAL_FFT_HH

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

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

// Simple FFT implementation using Cooley-Tukey algorithm
namespace FFT {
    
    // Forward FFT for real data
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

// GSL-compatible interface
struct gsl_fft_real_wavetable {
    size_t n;
};

struct gsl_fft_real_workspace {
    size_t n;
    std::vector<double> work;
};

inline gsl_fft_real_wavetable* gsl_fft_real_wavetable_alloc(size_t n) {
    gsl_fft_real_wavetable* w = new gsl_fft_real_wavetable;
    w->n = n;
    return w;
}

inline gsl_fft_real_workspace* gsl_fft_real_workspace_alloc(size_t n) {
    gsl_fft_real_workspace* w = new gsl_fft_real_workspace;
    w->n = n;
    w->work.resize(n);
    return w;
}

inline void gsl_fft_real_transform(double* data, size_t stride, size_t n,
                                    gsl_fft_real_wavetable* /* wavetable */,
                                    gsl_fft_real_workspace* /* workspace */) {
    FFT::fft_real_transform(data, stride, n);
}

inline void gsl_fft_real_wavetable_free(gsl_fft_real_wavetable* w) {
    delete w;
}

inline void gsl_fft_real_workspace_free(gsl_fft_real_workspace* w) {
    delete w;
}

// Halfcomplex FFT (for inverse transforms)
struct gsl_fft_halfcomplex_wavetable {
    size_t n;
};

struct gsl_fft_halfcomplex_workspace {
    size_t n;
    std::vector<double> work;
};

inline gsl_fft_halfcomplex_wavetable* gsl_fft_halfcomplex_wavetable_alloc(size_t n) {
    gsl_fft_halfcomplex_wavetable* w = new gsl_fft_halfcomplex_wavetable;
    w->n = n;
    return w;
}

inline gsl_fft_halfcomplex_workspace* gsl_fft_halfcomplex_workspace_alloc(size_t n) {
    gsl_fft_halfcomplex_workspace* w = new gsl_fft_halfcomplex_workspace;
    w->n = n;
    w->work.resize(n);
    return w;
}

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
    
    // Inverse FFT
    ifft(complex_data);
    
    // Convert back to real
    for (size_t i = 0; i < n; ++i) {
        data[i * stride] = complex_data[i].real();
    }
}

inline void gsl_fft_halfcomplex_wavetable_free(gsl_fft_halfcomplex_wavetable* w) {
    delete w;
}

inline void gsl_fft_halfcomplex_workspace_free(gsl_fft_halfcomplex_workspace* w) {
    delete w;
}

// Alias for inverse transform
inline void gsl_fft_halfcomplex_inverse(double* data, size_t stride, size_t n,
                                         gsl_fft_halfcomplex_wavetable* wavetable,
                                         gsl_fft_halfcomplex_workspace* workspace) {
    gsl_fft_halfcomplex_transform(data, stride, n, wavetable, workspace);
}

// Allow real workspace to be used with halfcomplex (they're compatible)
inline void gsl_fft_halfcomplex_inverse(double* data, size_t stride, size_t n,
                                         gsl_fft_halfcomplex_wavetable* /* wavetable */,
                                         gsl_fft_real_workspace* /* workspace */) {
    gsl_fft_halfcomplex_wavetable* hc_w = gsl_fft_halfcomplex_wavetable_alloc(n);
    gsl_fft_halfcomplex_workspace* hc_ws = gsl_fft_halfcomplex_workspace_alloc(n);
    gsl_fft_halfcomplex_transform(data, stride, n, hc_w, hc_ws);
    gsl_fft_halfcomplex_wavetable_free(hc_w);
    gsl_fft_halfcomplex_workspace_free(hc_ws);
}

// Complex FFT
struct gsl_fft_complex_wavetable {
    size_t n;
};

struct gsl_fft_complex_workspace {
    size_t n;
    std::vector<std::complex<double>> work;
};

inline gsl_fft_complex_wavetable* gsl_fft_complex_wavetable_alloc(size_t n) {
    gsl_fft_complex_wavetable* w = new gsl_fft_complex_wavetable();
    w->n = n;
    return w;
}

inline gsl_fft_complex_workspace* gsl_fft_complex_workspace_alloc(size_t n) {
    gsl_fft_complex_workspace* w = new gsl_fft_complex_workspace();
    w->n = n;
    w->work.resize(n);
    return w;
}

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

inline void gsl_fft_complex_inverse(double* data, size_t stride, size_t n,
                                    gsl_fft_complex_wavetable* /* wavetable */,
                                    gsl_fft_complex_workspace* /* workspace */) {
    std::vector<std::complex<double>> complex_data(n);
    for (size_t i = 0; i < n; ++i) {
        complex_data[i] = std::complex<double>(data[2 * i * stride], data[(2 * i + 1) * stride]);
    }
    ifft(complex_data);
    for (size_t i = 0; i < n; ++i) {
        data[2 * i * stride] = complex_data[i].real();
        data[(2 * i + 1) * stride] = complex_data[i].imag();
    }
}

inline void gsl_fft_complex_wavetable_free(gsl_fft_complex_wavetable* w) {
    delete w;
}

inline void gsl_fft_complex_workspace_free(gsl_fft_complex_workspace* w) {
    delete w;
}

// Radix-2 FFT (simplified version without wavetable/workspace)
inline void gsl_fft_complex_radix2_forward(double* data, size_t stride, size_t n) {
    gsl_fft_complex_wavetable* w = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* ws = gsl_fft_complex_workspace_alloc(n);
    gsl_fft_complex_forward(data, stride, n, w, ws);
    gsl_fft_complex_wavetable_free(w);
    gsl_fft_complex_workspace_free(ws);
}

inline void gsl_fft_complex_radix2_inverse(double* data, size_t stride, size_t n) {
    gsl_fft_complex_wavetable* w = gsl_fft_complex_wavetable_alloc(n);
    gsl_fft_complex_workspace* ws = gsl_fft_complex_workspace_alloc(n);
    gsl_fft_complex_inverse(data, stride, n, w, ws);
    gsl_fft_complex_wavetable_free(w);
    gsl_fft_complex_workspace_free(ws);
}

#endif // OPAL_FFT_HH

