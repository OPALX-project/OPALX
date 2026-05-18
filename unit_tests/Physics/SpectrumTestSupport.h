/**
 * @file SpectrumTestSupport.h
 * @brief Lightweight histogram + CSV utilities shared by decay spectrum tests.
 */

#ifndef OPALX_UNIT_TESTS_PHYSICS_SPECTRUM_TEST_SUPPORT_H
#define OPALX_UNIT_TESTS_PHYSICS_SPECTRUM_TEST_SUPPORT_H

#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace opalx::test {

    struct Histogram1D {
        double low                 = 0.0;
        double high                = 0.0;
        std::size_t nBins          = 0;
        std::vector<double> counts;   ///< raw sample counts per bin (out-of-range excluded)
        std::vector<double> density;  ///< counts / total / binWidth (PDF estimate)

        double binWidth() const { return (high - low) / static_cast<double>(nBins); }

        double center(std::size_t i) const {
            return low + (static_cast<double>(i) + 0.5) * binWidth();
        }

        double edgeLow(std::size_t i) const {
            return low + static_cast<double>(i) * binWidth();
        }

        double edgeHigh(std::size_t i) const {
            return low + static_cast<double>(i + 1) * binWidth();
        }
    };

    inline Histogram1D makeHistogram(
            const std::vector<double>& samples, double low, double high, std::size_t nBins) {
        if (nBins == 0 || !(high > low)) {
            throw std::invalid_argument("makeHistogram: invalid range or nBins");
        }
        Histogram1D h;
        h.low   = low;
        h.high  = high;
        h.nBins = nBins;
        h.counts.assign(nBins, 0.0);
        h.density.assign(nBins, 0.0);

        const double width = h.binWidth();
        double totalIn     = 0.0;
        for (double s : samples) {
            if (s < low || s >= high) {
                continue;
            }
            const std::size_t idx =
                    static_cast<std::size_t>((s - low) / width);
            if (idx < nBins) {
                h.counts[idx] += 1.0;
                totalIn += 1.0;
            }
        }
        if (totalIn > 0.0) {
            for (std::size_t i = 0; i < nBins; ++i) {
                h.density[i] = h.counts[i] / (totalIn * width);
            }
        }
        return h;
    }

    /// L1 distance between sampled density and an analytic PDF on the same binning.
    /// Both inputs are interpreted as densities; the integrand is |a - b| * binWidth.
    inline double l1Distance(const Histogram1D& sampled, const std::vector<double>& analyticPdf) {
        if (sampled.density.size() != analyticPdf.size()) {
            throw std::invalid_argument("l1Distance: size mismatch");
        }
        const double width = sampled.binWidth();
        double sum         = 0.0;
        for (std::size_t i = 0; i < sampled.density.size(); ++i) {
            sum += std::abs(sampled.density[i] - analyticPdf[i]) * width;
        }
        return sum;
    }

    /// Riemann sum of the histogram density (should be ~1 if all samples were in range).
    inline double histogramArea(const Histogram1D& h) {
        const double width = h.binWidth();
        double sum         = 0.0;
        for (double d : h.density) {
            sum += d * width;
        }
        return sum;
    }

    inline double sampleMean(const std::vector<double>& samples) {
        if (samples.empty()) {
            return 0.0;
        }
        double s = 0.0;
        for (double x : samples) {
            s += x;
        }
        return s / static_cast<double>(samples.size());
    }

    /// Write a CSV consumable by tools/spectrum_plots/plot_spectrum.py.
    /// Header: "# x_label: <xLabel>\n# columns: bin_low,bin_high,bin_center,density,count,analytic_pdf\n"
    inline void writeSpectrumCsv(
            const std::string& path, const Histogram1D& h,
            const std::vector<double>& analyticPdf, const std::string& xLabel) {
        if (analyticPdf.size() != h.nBins) {
            throw std::invalid_argument("writeSpectrumCsv: analyticPdf size mismatch");
        }
        std::ofstream out(path);
        if (!out) {
            throw std::runtime_error("writeSpectrumCsv: cannot open " + path);
        }
        out.precision(17);
        out << "# x_label: " << xLabel << "\n";
        out << "# columns: bin_low,bin_high,bin_center,density,count,analytic_pdf\n";
        for (std::size_t i = 0; i < h.nBins; ++i) {
            out << h.edgeLow(i) << ',' << h.edgeHigh(i) << ',' << h.center(i) << ','
                << h.density[i] << ',' << h.counts[i] << ',' << analyticPdf[i] << '\n';
        }
    }

}  // namespace opalx::test

#endif  // OPALX_UNIT_TESTS_PHYSICS_SPECTRUM_TEST_SUPPORT_H
