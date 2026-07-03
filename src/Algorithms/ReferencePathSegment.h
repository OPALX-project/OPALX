#ifndef OPALX_ReferencePathSegment_HH
#define OPALX_ReferencePathSegment_HH

#include "AbsBeamline/Component.h"

#include <optional>
#include <set>

/**
 * @class ReferencePathSegment
 * @brief Reference-coordinate interval with optional legacy ELEMEDGE anchor.
 *
 * This class belongs to the reference and ordering layer, not the geometric
 * placement layer. It can therefore preserve `ELEMEDGE` for backward
 * compatibility without making it the primary placement mechanism. In the
 * first redesign stage the segment can also carry the active element set for
 * the interval \f$[s_\mathrm{begin}, s_\mathrm{end}]\f$.
 */
class ReferencePathSegment {
public:
    using element_set_t = std::set<std::shared_ptr<Component>>;

    ReferencePathSegment(
            double begin = 0.0, double end = 0.0, const element_set_t& activeElements = {},
            const std::optional<double>& legacyElementEdge = std::nullopt)
        : begin_m(begin),
          end_m(end),
          activeElements_m(activeElements),
          legacyElementEdge_m(legacyElementEdge) {}

    double getBegin() const { return begin_m; }
    double getEnd() const { return end_m; }
    const element_set_t& getActiveElements() const { return activeElements_m; }
    bool hasLegacyElementEdge() const { return legacyElementEdge_m.has_value(); }
    const std::optional<double>& getLegacyElementEdge() const { return legacyElementEdge_m; }

private:
    double begin_m;
    double end_m;
    element_set_t activeElements_m;
    std::optional<double> legacyElementEdge_m;
};

#endif
