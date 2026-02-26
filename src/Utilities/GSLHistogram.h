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
#include <cstdio>
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

    /// \brief Set bin edges from an explicit array.
    /// \param range Input: array of \f$n+1\f$ bin edges.
    /// \param n Input: number of edges (must be \f$n_+1\f$).
    void set_ranges(const double* range, size_t n) {
        if (n != n_ + 1) {
            throw std::invalid_argument("gsl_histogram: range size mismatch");
        }
        range_.assign(range, range + n);
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

    /// \brief Number of bins.
    /// \return Output: bin count.
    size_t n() const { return n_; }
    /// \brief Bin edge array.
    /// \return Output: reference to edges of size \f$n+1\f$.
    const std::vector<double>& range() const { return range_; }
    /// \brief Bin count array.
    /// \return Output: reference to counts of size \f$n\f$.
    const std::vector<double>& bin() const { return bin_; }

    // Make members accessible for GSL compatibility functions
    size_t n_;
    std::vector<double> range_;
    std::vector<double> bin_;

private:
    size_t find_bin(double x) const {
        if (range_.empty())
            return bin_.size();

        // Handle out-of-range values
        if (x < range_[0])
            return bin_.size();
        if (x >= range_.back())
            return bin_.size();

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

/// \brief Number of bins in the histogram.
/// \param h Input: histogram.
/// \return Output: bin count.
inline size_t gsl_histogram_bins(const gsl_histogram* h) { return h->n(); }

/// \brief Get the \f$[x_i, x_{i+1})\f$ range for a bin.
/// \param h Input: histogram.
/// \param i Input: bin index.
/// \param lower Output: lower edge.
/// \param upper Output: upper edge.
/// \return Output: 0 on success, -1 on invalid index.
inline int gsl_histogram_get_range(const gsl_histogram* h, size_t i, double* lower, double* upper) {
    if (i >= h->n()) {
        return -1;  // Error: bin index out of range
    }
    if (h->range().empty() || i >= h->range().size() - 1) {
        return -1;
    }
    *lower = h->range()[i];
    *upper = h->range()[i + 1];
    return 0;
}

/// \brief Set explicit bin edges for a histogram.
/// \param h Input/Output: histogram to configure.
/// \param range Input: array of \f$n\f$ edges (must be \f$n_+1\f$).
/// \param n Input: number of edges.
/// \return Output: 0 on success, -1 on error.
inline int gsl_histogram_set_ranges(gsl_histogram* h, const double* range, size_t n) {
    h->set_ranges(range, n);
    return 0;
}

// 2D Histogram
/// \brief 2D histogram with explicit x/y bin edges.
class gsl_histogram2d {
public:
    gsl_histogram2d() : nx_(0), ny_(0), xrange_(), yrange_(), bin_() {}

    /// \brief Set uniform x/y bin edges over \f$[x_{\min},x_{\max}]\f$ and
    /// \f$[y_{\min},y_{\max}]\f$. \param xmin Input: x lower bound (inclusive). \param xmax Input:
    /// x upper bound (exclusive). \param ymin Input: y lower bound (inclusive). \param ymax Input:
    /// y upper bound (exclusive).
    void set_ranges_uniform(double xmin, double xmax, double ymin, double ymax) {
        if (xmin >= xmax || ymin >= ymax) {
            throw std::invalid_argument("gsl_histogram2d: invalid range");
        }
        xrange_.resize(nx_ + 1);
        yrange_.resize(ny_ + 1);
        double xbin_width = (xmax - xmin) / nx_;
        double ybin_width = (ymax - ymin) / ny_;
        for (size_t i = 0; i <= nx_; ++i) {
            xrange_[i] = xmin + i * xbin_width;
        }
        for (size_t i = 0; i <= ny_; ++i) {
            yrange_[i] = ymin + i * ybin_width;
        }
    }

    /// \brief Accumulate a weighted sample into the 2D histogram.
    /// \param x Input: x sample value.
    /// \param y Input: y sample value.
    /// \param weight Input: weight to add to the bin.
    void accumulate(double x, double y, double weight) {
        size_t i = find_x_bin(x);
        size_t j = find_y_bin(y);
        if (i < nx_ && j < ny_) {
            bin_[i * ny_ + j] += weight;
        }
    }

    /// \brief Increment the bin containing \p x,\p y by 1.
    /// \param x Input: x sample value.
    /// \param y Input: y sample value.
    void increment(double x, double y) { accumulate(x, y, 1.0); }

    /// \brief Set explicit x/y bin edges.
    /// \param xrange Input: x edges array (\f$n_x+1\f$).
    /// \param nx Input: number of x edges.
    /// \param yrange Input: y edges array (\f$n_y+1\f$).
    /// \param ny Input: number of y edges.
    void set_ranges(const double* xrange, size_t nx, const double* yrange, size_t ny) {
        if (nx != nx_ + 1 || ny != ny_ + 1) {
            throw std::invalid_argument("gsl_histogram2d: range size mismatch");
        }
        xrange_.assign(xrange, xrange + nx);
        yrange_.assign(yrange, yrange + ny);
    }

    /// \brief Get the bin count for index \f$(i,j)\f$.
    /// \param i Input: x bin index.
    /// \param j Input: y bin index.
    /// \return Output: bin count.
    double get(size_t i, size_t j) const {
        if (i >= nx_ || j >= ny_) {
            throw std::out_of_range("gsl_histogram2d: bin index out of range");
        }
        return bin_[i * ny_ + j];
    }

    /// \brief Number of x bins.
    /// \return Output: x bin count.
    size_t nx() const { return nx_; }
    /// \brief Number of y bins.
    /// \return Output: y bin count.
    size_t ny() const { return ny_; }
    /// \brief X-axis bin edges.
    /// \return Output: reference to x edges.
    const std::vector<double>& xrange() const { return xrange_; }
    /// \brief Y-axis bin edges.
    /// \return Output: reference to y edges.
    const std::vector<double>& yrange() const { return yrange_; }
    /// \brief Bin count array.
    /// \return Output: reference to counts of size \f$n_x n_y\f$.
    const std::vector<double>& bin() const { return bin_; }

    // Make members accessible for GSL compatibility functions
    size_t nx_, ny_;
    std::vector<double> xrange_;
    std::vector<double> yrange_;
    std::vector<double> bin_;

private:
    size_t find_x_bin(double x) const {
        if (xrange_.empty() || x < xrange_[0] || x >= xrange_.back()) {
            return nx_;
        }
        auto it = std::upper_bound(xrange_.begin(), xrange_.end(), x);
        return std::distance(xrange_.begin(), it) - 1;
    }

    size_t find_y_bin(double y) const {
        if (yrange_.empty() || y < yrange_[0] || y >= yrange_.back()) {
            return ny_;
        }
        auto it = std::upper_bound(yrange_.begin(), yrange_.end(), y);
        return std::distance(yrange_.begin(), it) - 1;
    }
};

// 2D Histogram PDF for sampling
/// \brief 2D PDF derived from histogram for sampling.
class gsl_histogram2d_pdf {
public:
    gsl_histogram2d_pdf() : nx_(0), ny_(0), xrange_(), yrange_(), pdf_(), cumsum_() {}

    /// \brief Initialize a PDF from a 2D histogram.
    /// \param h Input: histogram providing bin weights.
    void init(const gsl_histogram2d* h) {
        nx_     = h->nx();
        ny_     = h->ny();
        xrange_ = h->xrange();
        yrange_ = h->yrange();

        const auto& bin = h->bin();
        double total    = 0.0;
        for (double val : bin) {
            total += val;
        }

        pdf_.resize(bin.size());
        if (total > 0) {
            for (size_t i = 0; i < bin.size(); ++i) {
                pdf_[i] = bin[i] / total;
            }
        }

        // Compute cumulative distribution for sampling
        cumsum_.resize(pdf_.size());
        double cum = 0.0;
        for (size_t i = 0; i < pdf_.size(); ++i) {
            cum += pdf_[i];
            cumsum_[i] = cum;
        }
    }

    /// \brief Sample \f$(x,y)\f$ given two uniform variates.
    /// \param u Input: uniform variate in \f$[0,1]\f$ for x-bin selection.
    /// \param v Input: uniform variate in \f$[0,1]\f$ for y position within the x-bin.
    /// \param x Output: sampled x value.
    /// \param y Output: sampled y value.
    void sample(double u, double v, double* x, double* y) {
        // Find bin using u (for x) and v (for y within that x slice)
        size_t bin_idx = 0;
        double target  = u * cumsum_.back();
        for (size_t i = 0; i < cumsum_.size(); ++i) {
            if (cumsum_[i] >= target) {
                bin_idx = i;
                break;
            }
        }

        size_t i = bin_idx / ny_;
        size_t j = bin_idx % ny_;

        if (i < nx_ && j < ny_) {
            double xmin = xrange_[i];
            double xmax = xrange_[i + 1];
            double ymin = yrange_[j];
            double ymax = yrange_[j + 1];
            *x          = xmin + (xmax - xmin) * v;
            *y          = ymin + (ymax - ymin) * u;
        } else {
            *x = 0.0;
            *y = 0.0;
        }
    }

private:
    size_t nx_, ny_;
    std::vector<double> xrange_;
    std::vector<double> yrange_;
    std::vector<double> pdf_;
    std::vector<double> cumsum_;
};

// GSL-compatible functions for 2D histogram
/// \brief Allocate a 2D histogram with \p nx by \p ny bins.
/// \param nx Input: number of x bins.
/// \param ny Input: number of y bins.
/// \return Output: histogram pointer.
inline gsl_histogram2d* gsl_histogram2d_alloc(size_t nx, size_t ny) {
    gsl_histogram2d* h = new gsl_histogram2d();
    h->nx_             = nx;
    h->ny_             = ny;
    h->bin_.resize(nx * ny, 0.0);
    return h;
}

/// \brief Set uniform x/y bin edges.
/// \param h Input/Output: histogram to configure.
/// \param xmin Input: x lower bound (inclusive).
/// \param xmax Input: x upper bound (exclusive).
/// \param ymin Input: y lower bound (inclusive).
/// \param ymax Input: y upper bound (exclusive).
inline void gsl_histogram2d_set_ranges_uniform(
    gsl_histogram2d* h, double xmin, double xmax, double ymin, double ymax
) {
    h->set_ranges_uniform(xmin, xmax, ymin, ymax);
}

/// \brief Set explicit x/y bin edges.
/// \param h Input/Output: histogram to configure.
/// \param xrange Input: x edges array (\f$n_x+1\f$).
/// \param nx Input: number of x edges.
/// \param yrange Input: y edges array (\f$n_y+1\f$).
/// \param ny Input: number of y edges.
/// \return Output: 0 on success, -1 on error.
inline int gsl_histogram2d_set_ranges(
    gsl_histogram2d* h, const double* xrange, size_t nx, const double* yrange, size_t ny
) {
    h->set_ranges(xrange, nx, yrange, ny);
    return 0;
}

/// \brief Number of x bins.
/// \param h Input: histogram.
/// \return Output: x bin count.
inline size_t gsl_histogram2d_nx(const gsl_histogram2d* h) { return h->nx(); }

/// \brief Number of y bins.
/// \param h Input: histogram.
/// \return Output: y bin count.
inline size_t gsl_histogram2d_ny(const gsl_histogram2d* h) { return h->ny(); }

/// \brief Increment the bin containing \p x,\p y by 1.
/// \param h Input/Output: histogram to update.
/// \param x Input: x sample value.
/// \param y Input: y sample value.
inline void gsl_histogram2d_increment(gsl_histogram2d* h, double x, double y) {
    h->increment(x, y);
}

/// \brief Get the bin count for index \f$(i,j)\f$.
/// \param h Input: histogram.
/// \param i Input: x bin index.
/// \param j Input: y bin index.
/// \return Output: bin count.
inline double gsl_histogram2d_get(const gsl_histogram2d* h, size_t i, size_t j) {
    return h->get(i, j);
}

/// \brief Accumulate a weighted sample into the 2D histogram.
/// \param h Input/Output: histogram to update.
/// \param x Input: x sample value.
/// \param y Input: y sample value.
/// \param weight Input: weight to add.
inline void gsl_histogram2d_accumulate(gsl_histogram2d* h, double x, double y, double weight) {
    h->accumulate(x, y, weight);
}

/// \brief Print the 2D histogram in a simple text table.
/// \param fh Input: output file handle.
/// \param h Input: histogram.
/// \param xformat Input: printf-style format string for x values.
/// \param yformat Input: printf-style format string for y values.
inline void gsl_histogram2d_fprintf(
    FILE* fh, const gsl_histogram2d* h, const char* xformat, const char* yformat
) {
    fprintf(fh, "# xrange: %zu bins from %g to %g\n", h->nx(), h->xrange()[0], h->xrange().back());
    fprintf(fh, "# yrange: %zu bins from %g to %g\n", h->ny(), h->yrange()[0], h->yrange().back());
    const auto& bin = h->bin();
    for (size_t i = 0; i < h->nx(); ++i) {
        for (size_t j = 0; j < h->ny(); ++j) {
            fprintf(fh, xformat, h->xrange()[i]);
            fprintf(fh, " ");
            fprintf(fh, yformat, h->yrange()[j]);
            fprintf(fh, " %g\n", bin[i * h->ny() + j]);
        }
    }
}

/// \brief Free a 2D histogram.
/// \param h Input: histogram to release (can be null).
inline void gsl_histogram2d_free(gsl_histogram2d* h) { delete h; }

/// \brief Allocate a 2D histogram PDF.
/// \param nx Input: number of x bins (unused).
/// \param ny Input: number of y bins (unused).
/// \return Output: PDF pointer.
inline gsl_histogram2d_pdf* gsl_histogram2d_pdf_alloc(size_t /* nx */, size_t /* ny */) {
    return new gsl_histogram2d_pdf();
}

/// \brief Initialize a PDF from a histogram.
/// \param p Input/Output: PDF to initialize.
/// \param h Input: source histogram.
inline void gsl_histogram2d_pdf_init(gsl_histogram2d_pdf* p, const gsl_histogram2d* h) {
    p->init(h);
}

/// \brief Sample from the PDF using two uniform variates.
/// \param p Input: PDF.
/// \param u Input: uniform variate in \f$[0,1]\f$.
/// \param v Input: uniform variate in \f$[0,1]\f$.
/// \param x Output: sampled x value.
/// \param y Output: sampled y value.
inline void gsl_histogram2d_pdf_sample(
    const gsl_histogram2d_pdf* p, double u, double v, double* x, double* y
) {
    const_cast<gsl_histogram2d_pdf*>(p)->sample(u, v, x, y);
}

/// \brief Free a 2D histogram PDF.
/// \param p Input: PDF to release (can be null).
inline void gsl_histogram2d_pdf_free(gsl_histogram2d_pdf* p) { delete p; }

#endif  // OPAL_GSL_HISTOGRAM_HH
