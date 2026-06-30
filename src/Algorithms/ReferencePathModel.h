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

    /// @name Short-bunch s ⇄ local-Z seam
    /// The reference path is treated as locally straight, so a particle's longitudinal local-Z
    /// maps linearly to the path length s. **These are exact only for a straight reference path**
    /// and are the single home of that approximation — revisit for curved beamlines (bends).
    ///@{

    /// Path length s of a point at local-Z @p localZ when the reference/centroid is at @p sPos.
    /// Used by the APPLY stage to turn the bunch's local-Z bounds into an s-keyed IndexMap query.
    static double pathLengthFromLocalZ(double sPos, double localZ) { return sPos + localZ; }

    /// Global s of the element edge recovered from a particle's local entry-Z and momentum,
    /// projecting back along the trajectory. Used by the INDEX stage (OrbitThreader) to recover
    /// the legacy ELEMEDGE anchor. @p localPNorm is |P|, @p localPz is the local z-momentum.
    static double elementEdgeFromLocalEntry(
            double start, double localEntryZ, double localPNorm, double localPz) {
        return start - localEntryZ * localPNorm / localPz;
    }
    ///@}

private:
    container_type segments_m;
};

#endif
