#ifndef OPALX_SBend_HH
#define OPALX_SBend_HH

#include "AbsBeamline/BendFieldModel.h"
#include "AbsBeamline/ElementBase.h"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <string>
#include <vector>

/**
 * @class SBend
 * @brief Analytic sector bending magnet.
 *
 * The body is a circular arc (curvature and arc length from the Geometry); the
 * field is a vertical dipole with an Enge fringe, evaluated in arc coordinates
 * (radial offset, vertical, arc length) so the field support follows the curved
 * design orbit at any bend angle.
 */
class SBend : public ElementBase {
public:
    SBend();
    explicit SBend(const std::string& name);
    SBend(const SBend&);
    ~SBend() override;

    void accept(BeamlineVisitor& visitor) const override;
    ElementType getType() const override;

    void initialise(PartBunch_t* bunch) override;
    void finalise() override;

    bool apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    bool apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// @brief Calculate the SBend field extent in the element's coordinatesystem
    /// @param zBegin Where the field begins (negative for fringe fields)
    /// @param zEnd Where the field ends (larger than element length for fringe fields)
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// Containment in arc-length coordinates (so the bend stays selected as the
    /// reference orbit curves through it), not the straight-frame z.
    bool isInside(const Vector_t<double, 3>& r) const override;

    /// Full vertical pole gap (GAP), scaling the Enge fringe profile.
    void setFullGap(double gap);

    /// Half gap (HGAP) used by the Enge fringe profile and pole-face focusing.
    void setFringeHalfGap(double halfGap);

    /// Pole-face field integral (FINT) used by the vertical edge focusing.
    void setFringeIntegral(double fringeIntegral);

    void setDesignEnergy(const double& energy, bool changeable = true) override;
    double getDesignEnergy() const override;

    /// Store the normal/skew multipole coefficients into the device views read
    /// by apply()/computeFieldHost(). Values are taken as-is (already scaled by
    /// the caller).
    void setFieldComponents(const std::vector<double>& normal, const std::vector<double>& skew);

    /// Dipole normal component (the "BY" channel attribute), backed by the
    /// coefficient views.
    double getB() const;
    void setB(double B);

private:
    void computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const;

    /// Build the pure-value field inputs (coefficients, gap, pole-face projections,
    /// edge-focusing coefficients) shared by the device and host field paths.
    BendFieldModel::FieldInputs makeFieldInputs() const;

    /// Convert a local point to bend coordinates (radial x, y, arc-length s).
    Vector_t<double, 3> bendCoords(const Vector_t<double, 3>& r) const;

    /// Longitudinal containment: arc-length s within the field extent.
    bool isInsideArc(const Vector_t<double, 3>& arc) const;

    /// Normal/skew multipole coefficients on the device, read by tracking.
    Kokkos::View<double*> normalComponents_m;
    Kokkos::View<double*> skewComponents_m;
    int maxNormal_m = 0;
    int maxSkew_m   = 0;

    double gap_m;
    double fringeHalfGap_m;
    double fringeIntegral_m;
    double designEnergy_m;
    bool designEnergyChangeable_m;
};

inline void SBend::setFullGap(double gap) { gap_m = std::abs(gap); }

inline void SBend::setFringeHalfGap(double halfGap) { fringeHalfGap_m = std::abs(halfGap); }

inline void SBend::setFringeIntegral(double fringeIntegral) {
    fringeIntegral_m = std::abs(fringeIntegral);
}

inline void SBend::setDesignEnergy(const double& energy, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m           = std::abs(energy) * 1e6;
        designEnergyChangeable_m = changeable;
    }
}

inline double SBend::getDesignEnergy() const { return designEnergy_m; }

#endif  // OPALX_SBend_HH
