#ifndef OPALX_ReferencePathModel_HH
#define OPALX_ReferencePathModel_HH

#include "Algorithms/ReferencePathSegment.h"

#include <vector>

/**
 * @class ReferencePathModel
 * @brief Ordered collection of reference-path segments.
 *
 * The first redesign stage keeps this as a lightweight container that can later
 * be populated from IndexMap and legacy `ELEMEDGE` information. In the
 * terminology of the placement note, this belongs to the reference base
 * \f$\mathcal{B}\f$ and its reporting coordinate \f$s\f$, not to the
 * laboratory-space geometry layer.
 */
class ReferencePathModel {
public:
    using container_type = std::vector<ReferencePathSegment>;

    void addSegment(const ReferencePathSegment& segment) { segments_m.push_back(segment); }
    void clear() { segments_m.clear(); }
    size_t size() const { return segments_m.size(); }
    bool empty() const { return segments_m.empty(); }
    const container_type& getSegments() const { return segments_m; }
    const ReferencePathSegment& operator[](size_t idx) const { return segments_m[idx]; }

private:
    container_type segments_m;
};

#endif
