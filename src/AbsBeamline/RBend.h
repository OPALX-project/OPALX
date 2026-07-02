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
    /// @brief  Constructurs, Destructors, ...
    RBend();
    explicit RBend(const std::string& name);
    RBend(const RBend&);
    ~RBend() override;

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

    /// @brief Calculate the RBend field extent in the element's coordinatesystem
    /// @param zBegin Where the field begins (negative for fringe fields)
    /// @param zEnd Where the field ends (larger than element length for fringe fields)
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// Containment in the straight box frame: box z within the field extent and
    /// inside the transverse aperture.
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

    /// Design-orbit curvature 1/rho = angle / arc length = 2 sin(angle/2) / L. Used only to
    /// scale the pole-face edge-focusing kick; the box field itself is uniform.
    double referenceCurvature() const;
    /// Angle between the design orbit and the entrance pole face (angle/2 + E1), for the
    /// edge focusing. The box faces are perpendicular to the box axis, so the orbit meets
    /// them at half the bend angle plus any explicit pole-face rotation.
    double edgeAngleEntrance() const;
    /// Angle between the design orbit and the exit pole face (angle/2 + E2).
    double edgeAngleExit() const;

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

inline void RBend::setFullGap(double gap) { gap_m = std::abs(gap); }

inline void RBend::setFringeHalfGap(double halfGap) { fringeHalfGap_m = std::abs(halfGap); }

inline void RBend::setFringeIntegral(double fringeIntegral) {
    fringeIntegral_m = std::abs(fringeIntegral);
}

inline void RBend::setDesignEnergy(const double& energy, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m           = std::abs(energy) * 1e6;
        designEnergyChangeable_m = changeable;
    }
}

inline double RBend::getDesignEnergy() const { return designEnergy_m; }

#endif  // OPALX_RBend_HH
