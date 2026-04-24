#ifndef OPALX_ReferencePathSegment_HH
#define OPALX_ReferencePathSegment_HH

#include <optional>

/**
 * @class ReferencePathSegment
 * @brief Reference-coordinate interval with optional legacy ELEMEDGE anchor.
 *
 * This class belongs to the reference and ordering layer, not the geometric
 * placement layer. It can therefore preserve `ELEMEDGE` for backward
 * compatibility without making it the primary placement mechanism.
 */
class ReferencePathSegment {
public:
    ReferencePathSegment(
            double begin = 0.0, double end = 0.0,
            const std::optional<double>& legacyElementEdge = std::nullopt)
        : begin_m(begin), end_m(end), legacyElementEdge_m(legacyElementEdge) {}

    double getBegin() const { return begin_m; }
    double getEnd() const { return end_m; }
    bool hasLegacyElementEdge() const { return legacyElementEdge_m.has_value(); }
    const std::optional<double>& getLegacyElementEdge() const { return legacyElementEdge_m; }

private:
    double begin_m;
    double end_m;
    std::optional<double> legacyElementEdge_m;
};

#endif
