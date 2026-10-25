#include "AbsBeamline/CyclotronSector.h"
#include <cmath>
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

void CyclotronSector::configure(
        std::shared_ptr<const CyclotronSectorFieldMap> input, int symmetry, double low, double high,
        double factor, const std::vector<CyclotronTrimCoil>& models) {
    if (!input || symmetry < 2) {
        throw OpalException("CyclotronSector::configure", "A map and SYM >= 2 are required.");
    }
    map   = std::move(input);
    angle = 2 * std::acos(-1.0) / symmetry;
    if (symmetry < 2 || std::abs(map->thetaMin) > 1e-12
        || std::abs(map->nt * map->dtheta - angle) > 1e-8 || !(low < high)
        || !std::isfinite(low + high + factor))
        throw OpalException("CyclotronSector::configure", "Invalid sector bounds or map symmetry.");
    radius = map->rmin + 0.5 * (map->nr - 1) * map->dr;
    vmin   = low;
    vmax   = high;
    scale  = factor;
    // This arc is a placement/geometry representation, not an equilibrium orbit.
    geometry  = Geometry::makeSBend(radius * angle, 1 / radius);
    coils     = Kokkos::View<CyclotronTrimCoil*>("cyclotron_trim_coils", models.size());
    hostCoils = Kokkos::create_mirror_view(coils);
    for (size_t i = 0; i < models.size(); ++i)
        hostCoils(i) = models[i];
    Kokkos::deep_copy(coils, hostCoils);
}
void CyclotronSector::accept(BeamlineVisitor& visitor) const {
    visitor.visitCyclotronSector(*this);
}
void CyclotronSector::initialise(PartBunch_t* bunch) {
    if (!map) throw OpalException("CyclotronSector::initialise", "FMAPFN is required.");
    RefPartBunch_m = bunch;
}
bool CyclotronSector::isInside(const Vector_t<double, 3>& p) const {
    if (!map) return false;
    const double r = std::hypot(p[0] + radius, p[2]), theta = std::atan2(p[2], p[0] + radius);
    return r >= map->rmin && r <= map->rmin + (map->nr - 1) * map->dr && theta >= 0 && theta < angle
           && p[1] >= vmin && p[1] <= vmax;
}
BoundingBox CyclotronSector::getBoundingBoxInLabCoords() const {
    // Conservative box of the entire annulus; valid for arbitrary Mode-A rotations.
    std::vector<Vector_t<double, 3>> corners;
    const double outer = map->rmin + (map->nr - 1) * map->dr;
    for (double x : {-outer, outer})
        for (double z : {-outer, outer})
            for (double y : {vmin, vmax})
                corners.push_back(
                        getCSTrafoGlobal2Local().transformFrom(
                                Vector_t<double, 3>(x - radius, y, z)));
    return BoundingBox::getBoundingBox(corners);
}
void CyclotronSector::apply(
        const Vector_t<double, 3>& r, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>&, Vector_t<double, 3>& b) {
    cyclotron_field::field(
            map->host, hostCoils, map->rmin, map->dr, map->dtheta, map->nr, map->nt, radius, angle,
            vmin, vmax, scale, r, b);
}
void CyclotronSector::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    const auto grid = map->data;
    const auto trim = coils;
    const auto r = pc->R.getView(), b = pc->B.getView();
    const double rmin = map->rmin, dr = map->dr, dt = map->dtheta, rm = radius, a = angle,
                 lo = vmin, hi = vmax, f = scale;
    const int nr = map->nr, nt = map->nt;
    Kokkos::parallel_for(
            "CyclotronSector::apply", pc->getLocalNum(), KOKKOS_LAMBDA(size_t i) {
                cyclotron_field::field(
                        grid, trim, rmin, dr, dt, nr, nt, rm, a, lo, hi, f, r(i), b(i));
            });
}
