//
// Histogram implementation to replace GSL histogram
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

#ifndef OPAL_GSL_HISTOGRAM_HH
#define OPAL_GSL_HISTOGRAM_HH

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/histogram.html
/// \see Implementation on https://www.gnu.org/software/gsl/
// 1D Histogram
/// \brief 1D histogram with explicit bin edges.
class gsl_histogram {
public:
    gsl_histogram() : n_(0), range_(), bin_(0) {}

    /// \brief Set uniform bin edges over \f$[x_{\min}, x_{\max}]\f$.
    /// \param xmin Input: lower bound (inclusive).
    /// \param xmax Input: upper bound (exclusive).
    void set_ranges_uniform(double xmin, double xmax) {
        if (xmin >= xmax) {
            throw std::invalid_argument("gsl_histogram: xmin must be < xmax");
        }
        range_.resize(n_ + 1);
        double bin_width = (xmax - xmin) / n_;
        for (size_t i = 0; i <= n_; ++i) {
            range_[i] = xmin + i * bin_width;
        }
    }

    /// \brief Increment the bin containing \p x by 1.
    /// \param x Input: sample value.
    void increment(double x) {
        size_t bin = find_bin(x);
        if (bin < bin_.size()) {
            bin_[bin]++;
        }
    }

    /// \brief Get the bin count for index \p i.
    /// \param i Input: bin index.
    /// \return Output: bin count.
    double get(size_t i) const {
        if (i >= bin_.size()) {
            throw std::out_of_range("gsl_histogram: bin index out of range");
        }
        return bin_[i];
    }

    // Make members accessible for GSL compatibility functions
    size_t n_;
    std::vector<double> range_;
    std::vector<double> bin_;

private:
    size_t find_bin(double x) const {
        if (range_.empty()) return bin_.size();

        // Handle out-of-range values
        if (x < range_[0]) return bin_.size();
        if (x >= range_.back()) return bin_.size();

        // Binary search for the bin
        auto it = std::upper_bound(range_.begin(), range_.end(), x);
        return std::distance(range_.begin(), it) - 1;
    }
};

// GSL-compatible functions for 1D histogram
/// \brief Allocate a 1D histogram with \p n bins.
/// \param n Input: number of bins.
/// \return Output: histogram pointer.
inline gsl_histogram* gsl_histogram_alloc(size_t n) {
    gsl_histogram* h = new gsl_histogram();
    h->n_            = n;
    h->bin_.resize(n, 0.0);
    return h;
}

/// \brief Set uniform bin edges over \f$[x_{\min}, x_{\max}]\f$.
/// \param h Input/Output: histogram to configure.
/// \param xmin Input: lower bound (inclusive).
/// \param xmax Input: upper bound (exclusive).
inline void gsl_histogram_set_ranges_uniform(gsl_histogram* h, double xmin, double xmax) {
    h->set_ranges_uniform(xmin, xmax);
}

/// \brief Increment the bin containing \p x by 1.
/// \param h Input/Output: histogram to update.
/// \param x Input: sample value.
inline void gsl_histogram_increment(gsl_histogram* h, double x) { h->increment(x); }

/// \brief Get the bin count for index \p i.
/// \param h Input: histogram.
/// \param i Input: bin index.
/// \return Output: bin count.
inline double gsl_histogram_get(const gsl_histogram* h, size_t i) { return h->get(i); }

/// \brief Free a histogram allocated by \c gsl_histogram_alloc.
/// \param h Input: histogram to release (can be null).
inline void gsl_histogram_free(gsl_histogram* h) { delete h; }

#endif  // OPAL_GSL_HISTOGRAM_HH
