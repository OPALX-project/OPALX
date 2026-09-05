#ifndef OPALX_CYCLOTRON_SECTOR_H
#define OPALX_CYCLOTRON_SECTOR_H
#include "AbsBeamline/ElementBase.h"
#include "Fields/CyclotronSectorFieldMap.h"

/** Legacy PSI mirrored trim profile, expressed in metres, tesla and inverse metres.
 * RMIN/RMAX define its shape: the profile has radial tails. The radial correction
 * preserves the legacy derivative convention on both halves of the profile.
 */
struct CyclotronTrimCoil {
    double rmin = 0, rmax = 0, bmax = 0, slope = 0;
    KOKKOS_INLINE_FUNCTION void add(double r, double y, double& br, double& by) const {
        if (bmax == 0) return;
        const double x = (2 * r < rmin + rmax ? r - rmin : rmax - r) * slope;
        const double a = Kokkos::pow(10.0, (x - 4) * 0.03);
        const double b = Kokkos::pow(10.0, -(x - 8) * 0.2);
        const double u = 1 / (1 + a), v = 1 / (1 + b);
        // Use u*(1-u) to avoid overflow in the tails.
        const double deriv =
                bmax / 2.78 * slope * Kokkos::log(10.0) * (-0.09 * u * (1 - u) + v * (1 - v));
        by += bmax / 2.78 * (-3 + 3 * u + 5 * v);
        br += deriv * y;
    }
};

/** Annular magnetic sector in the X-Z plane (Y vertical).
 * The pose is the entrance tangent at the map midpoint radius. Cylindrical
 * interpolation and the first-order off-plane extension use immutable map data.
 */
class CyclotronSector : public ElementBase {
public:
    explicit CyclotronSector(const std::string& name) : ElementBase(name) {}
    CyclotronSector* clone() const override { return new CyclotronSector(*this); }
    void accept(BeamlineVisitor&) const override;
    ElementType getType() const override { return ElementType::CYCLOTRONSECTOR; }
    Geometry& getGeometry() override { return geometry; }
    const Geometry& getGeometry() const override { return geometry; }
    void initialise(PartBunch_t* bunch) override;
    void finalise() override { RefPartBunch_m = nullptr; }
    void getFieldExtent(double& begin, double& end) const override {
        begin = 0;
        end   = geometry.getArcLength();
    }
    bool isInside(const Vector_t<double, 3>& r) const override;
    bool isInsideBody(const Vector_t<double, 3>& r) const override { return isInside(r); }
    BoundingBox getBoundingBoxInLabCoords() const override;
    void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    void apply(
            const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, const double& t,
            Vector_t<double, 3>& e, Vector_t<double, 3>& b) override;
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, const double& t,
            Vector_t<double, 3>& e, Vector_t<double, 3>& b) override {
        apply(r, p, t, e, b);
        return false;
    }
    void configure(
            std::shared_ptr<const CyclotronSectorFieldMap> map, int symmetry, double vmin,
            double vmax, double scale, const std::vector<CyclotronTrimCoil>& coils);

private:
    Geometry geometry;
    std::shared_ptr<const CyclotronSectorFieldMap> map;
    Kokkos::View<CyclotronTrimCoil*> coils;
    decltype(Kokkos::create_mirror_view(
            std::declval<Kokkos::View<CyclotronTrimCoil*>>())) hostCoils;
    double radius = 0, angle = 0, vmin = 0, vmax = 0, scale = 1;
};
#endif
