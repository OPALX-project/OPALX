#ifndef OPALX_RBend_HH
#define OPALX_RBend_HH

#include "AbsBeamline/BendFieldModel.h"
#include "AbsBeamline/ElementBase.h"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <string>
#include <vector>

/**
 * @class RBend
 * @brief Common OPALX interface for analytic horizontal bending magnets.
 *
 * This class implements a rectangular element with a dipole field.
 */
class RBend : public ElementBase {
public:
    /// @brief Default constructor.
    RBend();

    /// @brief Constructor with given name.
    /// @param name The element name.
    explicit RBend(const std::string& name);

    /// @brief Copy constructor.
    RBend(const RBend&);

    /// @brief Destructor.
    ~RBend() override;

    /// @brief Apply a visitor.
    /// @param visitor The visitor to apply.
    void accept(BeamlineVisitor& visitor) const override;

    /// @brief Get the element type.
    /// @return The element type (ElementType::RBEND).
    ElementType getType() const override;

    /// @brief Set up the element for tracking.
    /// @param bunch The reference particle bunch.
    void initialise(PartBunch_t* bunch) override;

    /// @brief Clean up after tracking.
    void finalise() override;

    /// @brief Apply the field to all particles.
    /// @param pc The particle container.
    void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /// @brief Apply the field to a particle with position R and momentum P.
    /// @param R Position.
    /// @param P Momentum.
    /// @param t Time.
    /// @param E Electric field.
    /// @param B Magnetic field.
    void apply(
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

    /// @brief Calculate the RBend field extent in the element's coordinatesystem
    /// @param zBegin Where the field begins (negative for fringe fields)
    /// @param zEnd Where the field ends (larger than element length for fringe fields)
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// @brief Containment in the straight box frame.
    /// @note Box z within the field extent and inside the transverse aperture.
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
    /// @param energy The design energy in eV.
    /// @param changeable Whether the design energy can be changed later.
    void setDesignEnergy(const double& energy, bool changeable = true) override;

    /// @brief Get the design energy.
    /// @return The design energy in eV.
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

    /// @brief Build the pure-value field inputs for the kernel launch.
    /// @note Includes coefficients, gap, pole-face projections and edge-focusing coefficients.
    /// @return The field inputs.
    BendFieldModel::FieldInputs makeFieldInputs() const;

    /// @brief Design-orbit curvature 1/rho = angle / arc length = 2 sin(angle/2) / L.
    /// @note Used only to scale the pole-face edge-focusing kick; the box field itself
    ///       is uniform.
    /// @return The design-orbit curvature.
    double referenceCurvature() const;

    /// @brief Angle between the design orbit and the entrance pole face (angle/2 + E1),
    ///        for the edge focusing.
    /// @note The box faces are perpendicular to the box axis, so the orbit meets them at
    ///       half the bend angle plus any explicit pole-face rotation.
    /// @return The entrance edge angle.
    double edgeAngleEntrance() const;

    /// @brief Angle between the design orbit and the exit pole face (angle/2 + E2).
    /// @return The exit edge angle.
    double edgeAngleExit() const;

    /// Normal/skew multipole coefficients on the device, read by tracking.
    Kokkos::View<double*> normalComponents_m;
    Kokkos::View<double*> skewComponents_m;

    /// Host mirrors of the coefficient views. Required for kernel launch.
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

inline void RBend::setFullGap(double gap) { gap_m = std::abs(gap); }

inline void RBend::setFringeIntegral(double fringeIntegral) {
    fringeIntegral_m = std::abs(fringeIntegral);
}

inline void RBend::setDesignEnergy(const double& energy, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m           = std::abs(energy);
        designEnergyChangeable_m = changeable;
    }
}

inline double RBend::getDesignEnergy() const { return designEnergy_m; }

#endif  // OPALX_RBend_HH
