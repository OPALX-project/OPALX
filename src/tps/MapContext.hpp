#ifndef TPS_MAP_CONTEXT_HPP
#define TPS_MAP_CONTEXT_HPP
#include <Kokkos_Array.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace tps {
    /** @brief Meaning of the coordinates, declared by the adapter (not inferred). */
    enum class CoordinateConvention { unspecified, canonical, mechanical };
    /** @brief One local coordinate; physical value = reference + scale * local value. */
    struct Coordinate {
        std::string name;
        std::string unit;
        double scale                             = 1.0;
        bool operator==(const Coordinate&) const = default;
    };
    /** @brief Right-handed orthonormal frame; axes are row-major, origin uses metres. */
    struct ReferenceFrame {
        std::string name;
        Kokkos::Array<double, 3> origin{};
        Kokkos::Array<double, 9> axes{1, 0, 0, 0, 1, 0, 0, 0, 1};
        bool operator==(const ReferenceFrame&) const = default;
    };
    /** @brief A named physical parameter defining a chart, included in exact boundary matching.
     * For example, normalized magnetic-field components specify a local canonical gauge.
     * Parameters describe coordinate meaning, not the map's algorithm or integration settings.
     */
    struct ChartParameter {
        std::string name;
        std::string unit;
        double value                                 = 0;
        bool operator==(const ChartParameter&) const = default;
    };
    /** @brief A labelled expansion chart, including reference values and normalization.
     * @details reference is in the corresponding Coordinate::unit. referenceTime is
     * in seconds; referenceMomentum is in momentumUnit (no automatic conversion).
     * canonicalPairs lists (q,p) coordinate indices when convention is canonical.
     * Declaring canonical coordinates does not prove the supplied map is symplectic.
     */
    struct PhaseSpaceChart {
        std::string name;
        std::vector<Coordinate> coordinates;
        std::vector<double> reference;
        ReferenceFrame frame;
        std::optional<double> referenceTime;
        std::optional<double> referenceMomentum;
        std::string momentumUnit;
        CoordinateConvention convention = CoordinateConvention::unspecified;
        std::vector<std::pair<std::size_t, std::size_t>> canonicalPairs;
        std::vector<ChartParameter> parameters;  ///< Additional coordinate-defining values.
        bool operator==(const PhaseSpaceChart&) const = default;
    };
    /** @brief Validated, value-owned map boundaries, finite order and segment provenance.
     * @details No unit conversion, reference shift or frame transformation is implicit.
     * Boundary matching is exact metadata equality. Frame validation uses an absolute
     * tolerance of 1e-12 for orthonormality; this is not a tracking error tolerance.
     */
    class MapContext {
    public:
        MapContext(
                PhaseSpaceChart input, PhaseSpaceChart output, int order,
                std::vector<std::string> segments = {})
            : input_(std::move(input)),
              output_(std::move(output)),
              order_(order),
              segments_(std::move(segments)) {
            validate(input_);
            validate(output_);
            if (order < 0 || order > std::numeric_limits<int>::max() - 2)
                throw std::invalid_argument("MapContext: finite nonnegative order required");
            for (std::size_t i = 0; i < segments_.size(); ++i) {
                if (segments_[i].empty())
                    throw std::invalid_argument("MapContext: empty segment ID");
                for (std::size_t j = 0; j < i; ++j)
                    if (segments_[i] == segments_[j])
                        throw std::invalid_argument("MapContext: duplicate segment ID");
            }
        }
        [[nodiscard]] const PhaseSpaceChart& input() const noexcept { return input_; }
        [[nodiscard]] const PhaseSpaceChart& output() const noexcept { return output_; }
        [[nodiscard]] int order() const noexcept { return order_; }
        [[nodiscard]] const std::vector<std::string>& segments() const noexcept {
            return segments_;
        }

    private:
        static void validate(const PhaseSpaceChart& chart) {
            const auto n = chart.coordinates.size();
            if (chart.name.empty() || n == 0 || n > 99 || chart.reference.size() != n
                || chart.frame.name.empty())
                throw std::invalid_argument(
                        "MapContext: invalid chart name, dimension, reference or frame");
            for (std::size_t i = 0; i < n; ++i) {
                const auto& coordinate = chart.coordinates[i];
                if (coordinate.name.empty() || coordinate.unit.empty()
                    || !std::isfinite(coordinate.scale) || coordinate.scale <= 0
                    || !std::isfinite(chart.reference[i]))
                    throw std::invalid_argument("MapContext: invalid coordinate metadata");
                for (std::size_t j = 0; j < i; ++j)
                    if (coordinate.name == chart.coordinates[j].name)
                        throw std::invalid_argument("MapContext: duplicate coordinate name");
            }
            for (double v : chart.frame.origin)
                if (!std::isfinite(v)) throw std::invalid_argument("MapContext: nonfinite frame");
            for (double v : chart.frame.axes)
                if (!std::isfinite(v)) throw std::invalid_argument("MapContext: nonfinite axes");
            const auto& a = chart.frame.axes;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double dot = 0;
                    for (int k = 0; k < 3; ++k)
                        dot += a[3 * i + k] * a[3 * j + k];
                    if (std::abs(dot - (i == j ? 1.0 : 0.0)) > 1e-12)
                        throw std::invalid_argument("MapContext: nonorthonormal axes");
                }
            const double determinant = a[0] * (a[4] * a[8] - a[5] * a[7])
                                       - a[1] * (a[3] * a[8] - a[5] * a[6])
                                       + a[2] * (a[3] * a[7] - a[4] * a[6]);
            if (determinant <= 0) throw std::invalid_argument("MapContext: left-handed frame");
            if (chart.referenceTime && !std::isfinite(*chart.referenceTime))
                throw std::invalid_argument("MapContext: nonfinite time");
            if (chart.referenceMomentum) {
                if (!std::isfinite(*chart.referenceMomentum) || *chart.referenceMomentum <= 0
                    || chart.momentumUnit.empty())
                    throw std::invalid_argument("MapContext: invalid reference momentum");
            } else if (!chart.momentumUnit.empty())
                throw std::invalid_argument("MapContext: momentum unit without reference");
            if (chart.convention == CoordinateConvention::canonical) {
                if (chart.canonicalPairs.size() * 2 != n)
                    throw std::invalid_argument("MapContext: canonical pairs must cover chart");
                std::vector<bool> used(n, false);
                for (auto [q, p] : chart.canonicalPairs) {
                    if (q >= n || p >= n || q == p || used[q] || used[p])
                        throw std::invalid_argument("MapContext: invalid canonical pairing");
                    used[q] = used[p] = true;
                }
            } else if (!chart.canonicalPairs.empty())
                throw std::invalid_argument("MapContext: pairs require canonical convention");
            for (std::size_t i = 0; i < chart.parameters.size(); ++i) {
                const auto& parameter = chart.parameters[i];
                if (parameter.name.empty() || parameter.unit.empty()
                    || !std::isfinite(parameter.value))
                    throw std::invalid_argument("MapContext: invalid chart parameter");
                for (std::size_t j = 0; j < i; ++j)
                    if (parameter.name == chart.parameters[j].name)
                        throw std::invalid_argument("MapContext: duplicate chart parameter");
            }
        }
        PhaseSpaceChart input_, output_;
        int order_;
        std::vector<std::string> segments_;
    };
}  // namespace tps
#endif
