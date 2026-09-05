#ifndef OPALX_CYCLOTRON_SECTOR_H
#define OPALX_CYCLOTRON_SECTOR_H
#include "AbsBeamline/ElementBase.h"
#include "Fields/CyclotronSectorFieldMap.h"

/**
 * @brief Device-copyable parameters for OPAL 2022.1's PSI-BFIELD-MIRRORED profile.
 *
 * RMIN/RMAX define the profile shape, not hard cutoffs; the profile has radial
 * tails. This first model covers the full azimuth, with no azimuthal gate or
 * base-field threshold. Its strength is independent of the sector's BSCALE.
 *
 * @warning Compatibility deliberately preserves the legacy radial correction on
 * both halves: its sign is not reversed on the outer half when mirroring the
 * profile. It is not the derivative of the mirrored vertical field there.
 * See OPAL 2022.1 src/Classic/TrimCoils/TrimCoilMirrored.cpp for the source model.
 */
struct CyclotronTrimCoil {
    double rmin = 0, rmax = 0; ///< Inner/outer shape radii [m], with 0 < rmin < rmax.
    double bmax = 0; ///< Signed profile strength [T]; zero disables the contribution.
    double slope = 0; ///< Positive radial profile slope [1/m].
    /** @brief Add the legacy coil correction in the OPALX cylindrical frame.
     * @param r Distance from the sector's ring centre [m].
     * @param y Vertical displacement from the median plane [m].
     * @param[in,out] br Radial magnetic field [T], incremented in place.
     * @param[in,out] by Vertical magnetic field [T], incremented in place.
     * No azimuthal field is added. For q=slope*min(r-rmin,rmax-r), define
     * u=1/(1+10^(0.03*(q-4))) and v=1/(1+10^(-0.2*(q-8))). Then
     * @f$\Delta B_y=B_{max}(-3+3u+5v)/2.78@f$ and
     * @f$\Delta B_r=y B_{max} S\ln(10)[-0.09u(1-u)+v(1-v)]/2.78@f$,
     * where S is slope. Parameters must be validated before device execution.
     */
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

/**
 * @brief Static annular magnetic sector in the X-Z plane, with Y vertical.
 *
 * The element pose is the entrance tangent at Rm=(RMIN+RMAX)/2. In its local
 * entrance frame the ring centre is (-Rm,0,0), r=hypot(x+Rm,z), and
 * theta=atan2(z,x+Rm). Support is RMIN<=r<=RMAX, 0<=theta<2*pi/SYM, and
 * VMIN<=y<=VMAX. Half-open angular ownership prevents double counting at seams.
 * Poses are composed with the enclosing beamline frame by PlacementResolver.
 *
 * With stored median-plane map F [T] and scale s, the base field is
 * @f$B_y=sF,\ B_r=sy\partial_rF,\ B_\theta=sy\partial_\theta F/r@f$.
 * The cylindrical components are rotated into the local Cartesian frame and
 * accumulated; OpalBeamline applies the local/lab rotation. The proper rotation
 * from old OPAL coordinates is (x,y,z)_old -> (x,-z,y)_OPALX for both polar and
 * axial vectors. This explains the opposite vertical sign relative to old Bz.
 *
 * @note Only the legacy first-order off-plane extension is implemented. This is
 * not a general 3D field map. The geometry arc length Rm*2*pi/SYM is nominal and
 * must not be used to determine a physical cyclotron turn.
 */
class CyclotronSector : public ElementBase {
public:
    explicit CyclotronSector(const std::string& name) : ElementBase(name) {}
    /// Clone the occurrence, sharing map and read-only coil views; no grid copy.
    CyclotronSector* clone() const override { return new CyclotronSector(*this); }
    void accept(BeamlineVisitor&) const override;
    ElementType getType() const override { return ElementType::CYCLOTRONSECTOR; }
    Geometry& getGeometry() override { return geometry; }
    const Geometry& getGeometry() const override { return geometry; }
    /// Attach the bunch; throws OpalException if configure() has not loaded FMAPFN.
    void initialise(PartBunch_t* bunch) override;
    /// Detach the bunch without releasing the shared field data.
    void finalise() override { RefPartBunch_m = nullptr; }
    /// Nominal geometry interval [0,Rm*2*pi/SYM] in metres, not a Cartesian z bound.
    void getFieldExtent(double& begin, double& end) const override {
        begin = 0;
        end   = geometry.getArcLength();
    }
    /// True if an element-local position [m] lies within the annular field support.
    bool isInside(const Vector_t<double, 3>& r) const override;
    /// Nominal body and field support coincide for this hard-edge sector.
    bool isInsideBody(const Vector_t<double, 3>& r) const override { return isInside(r); }
    /// Conservative lab-frame box of the complete annulus, used for orbit threading.
    BoundingBox getBoundingBoxInLabCoords() const override;
    /** @brief Accumulate B [T] for local particles with R in the element frame [m].
     * Captures device views/parameters by value and launches a Kokkos kernel in
     * the default execution space. It leaves E unchanged, performs no MPI
     * reductions, and does not transfer map data or particle data to the host.
     * The caller owns coordinate transforms and kernel synchronization.
     */
    void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    /** @brief Scalar counterpart of the device kernel; adds to B, leaves E unchanged.
     * @param r Element-local position [m]; outside support contributes zero.
     * @param p Mechanical momentum [beta*gamma], unused for this static field.
     * @param t Time [s], unused for this static field.
     * @param[in,out] e Electric field [V/m], unchanged.
     * @param[in,out] b Element-local magnetic field [T], incremented in place.
     */
    void apply(
            const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, const double& t,
            Vector_t<double, 3>& e, Vector_t<double, 3>& b) override;
    /// Same scalar field; returns false (no element-local loss decision).
    /// Support selection and global tracking bounds are handled by the caller.
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, const double& t,
            Vector_t<double, 3>& e, Vector_t<double, 3>& b) override {
        apply(r, p, t, e, b);
        return false;
    }
    /** @brief Configure geometry and coil views on the host before tracking.
     * @param map Shared SI map; Kokkos must already be initialized.
     * @param symmetry Integer >=2, matching the map's one-sector angular span.
     * @param vmin Lower inclusive vertical limit [m].
     * @param vmax Upper inclusive vertical limit [m], greater than vmin.
     * @param scale Finite base-map multiplier (not applied to trim coils).
     * @param coils Validated full-azimuth coil parameter records, copied to device.
     * @throws OpalException If the map is absent, its angular origin is nonzero,
     * its span does not match symmetry, or vertical/scaling parameters are invalid.
     * No mutation or reconfiguration is allowed while field kernels are in flight.
     */
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
