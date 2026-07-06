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
    /// @brief Default constructor.
    SBend();

    /// @brief Constructor with given name.
    /// @param name The element name.
    explicit SBend(const std::string& name);

    /// @brief Copy constructor.
    SBend(const SBend&);

    /// @brief Destructor.
    ~SBend() override;

    /// @brief Apply a visitor.
    /// @param visitor The visitor to apply.
    void accept(BeamlineVisitor& visitor) const override;

    /// @brief Get the element type.
    /// @return The element type (ElementType::SBEND).
    ElementType getType() const override;

    /// @brief Set up the element for tracking.
    /// @param bunch The reference particle bunch.
    void initialise(PartBunch_t* bunch) override;

    /// @brief Clean up after tracking.
    void finalise() override;

    /// @brief Apply the field to all particles.
    /// @param pc The particle container.
    /// @return True if a particle is out-of-bounds (lost), false otherwise.
    bool apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /// @brief Apply the field to a particle with position R and momentum P.
    /// @param R Position.
    /// @param P Momentum.
    /// @param t Time.
    /// @param E Electric field.
    /// @param B Magnetic field.
    /// @return True if the particle is out-of-bounds (lost), false otherwise.
    bool apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// @brief Apply the field to the reference particle with position R and momentum P.
    /// @param R Position.
    /// @param P Momentum.
    /// @param t Time.
    /// @param E Electric field.
    /// @param B Magnetic field.
    /// @return True if the particle is out-of-bounds (lost), false otherwise.
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// @brief Calculate the SBend field extent in the element's coordinatesystem
    /// @param zBegin Where the field begins (negative for fringe fields)
    /// @param zEnd Where the field ends (larger than element length for fringe fields)
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// @brief Containment in arc-length coordinates.
    /// @note Uses arc-length (so the bend stays selected as the reference orbit curves
    ///       through it), not the straight-frame z.
    /// @param r The point to test.
    /// @return True if r is inside.
    bool isInside(const Vector_t<double, 3>& r) const override;

    /// @brief Set the full vertical pole gap (GAP), scaling the Enge fringe profile.
    /// @param gap The full vertical pole gap.
    void setFullGap(double gap);

    /// @brief Set the pole-face field integral (FINT) used by the vertical edge focusing.
    /// @param fringeIntegral The pole-face field integral.
    void setFringeIntegral(double fringeIntegral);

    /// @brief Set the design energy.
    /// @param energy The design energy.
    /// @param changeable Whether the design energy can be changed later.
    void setDesignEnergy(const double& energy, bool changeable = true) override;

    /// @brief Get the design energy.
    /// @return The design energy.
    double getDesignEnergy() const override;

    /// @brief Store the normal/skew multipole coefficients into the device views
    ///        read by apply()/computeFieldHost().
    /// @note Values are taken as-is (already scaled by the caller).
    /// @param normal The normal multipole coefficients.
    /// @param skew The skew multipole coefficients.
    void setFieldComponents(const std::vector<double>& normal, const std::vector<double>& skew);

    /// @brief Get the dipole normal component (the "BY" channel attribute).
    /// @note Backed by the coefficient views.
    /// @return The dipole normal component.
    double getB() const;

    /// @brief Set the dipole normal component (the "BY" channel attribute).
    /// @param B The dipole normal component.
    void setB(double B);

private:
    /// @brief Compute the field on the host at position R.
    /// @param R Position.
    /// @param B Magnetic field (output).
    void computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const;

    /// @brief Build the pure-value field inputs shared by the device and host field paths.
    /// @note Includes coefficients, gap, pole-face projections and edge-focusing coefficients.
    /// @return The field inputs.
    BendFieldModel::FieldInputs makeFieldInputs() const;

    /// @brief Convert a local point to bend coordinates (radial x, y, arc-length s).
    /// @param r The local point.
    /// @return The point in bend coordinates.
    Vector_t<double, 3> bendCoords(const Vector_t<double, 3>& r) const;

    /// @brief Longitudinal containment: arc-length s within the field extent.
    /// @param arc The point in bend coordinates.
    /// @return True if the arc-length s is within the field extent.
    bool isInsideArc(const Vector_t<double, 3>& arc) const;

    /// Normal/skew multipole coefficients on the device, read by tracking.
    Kokkos::View<double*> normalComponents_m;
    Kokkos::View<double*> skewComponents_m;

    /// Host mirrors of the coefficient views, filled once in setFieldComponents so the
    /// per-apply makeFieldInputs() reads them without a device->host copy.
    Kokkos::View<double*>::host_mirror_type normalComponentsHost_m;
    Kokkos::View<double*>::host_mirror_type skewComponentsHost_m;

    /// Number of stored normal/skew coefficients (sizes of the views above).
    int maxNormal_m = 0;
    int maxSkew_m   = 0;

    /// Full vertical pole gap (GAP).
    double gap_m;

    /// Pole-face field integral (FINT).
    double fringeIntegral_m;

    /// Design energy [eV].
    double designEnergy_m;

    /// Whether the design energy may still be changed.
    bool designEnergyChangeable_m;
};

inline void SBend::setFullGap(double gap) { gap_m = std::abs(gap); }


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
