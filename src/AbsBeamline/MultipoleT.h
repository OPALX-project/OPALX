//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef ABSBEAMLINE_MULTIPOLET_H
#define ABSBEAMLINE_MULTIPOLET_H

#include <memory>
#include <string>
#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/MultipoleTFieldModel.h"
#include "Algorithms/AbstractTimeDependence.h"

class AbstractTimeDependence;

/**
 * @class MultipoleT
 * @brief Combined-function multipole with a tanh fringe, straight or curved.
 *
 * The mid-plane field is @f$B_y = T(x)\,S(s)@f$ with the transverse profile
 * @f$T(x) = \sum_k b_k x^k@f$ (the @c TP coefficients: dipole, gradient, and so
 * on, in tesla per metre to the power k) and the tanh fringe @f$S(s)@f$ of the
 * two fringe lengths. The off-plane terms follow from the scalar potential; the
 * math lives in MultipoleTFieldModel.
 *
 * A non-zero bend angle makes the body a planar arc: the geometry is a sector
 * bend, the field is evaluated in arc coordinates, and the design orbit turns
 * with the centre of curvature on the -x side, exactly as for an SBEND. The
 * dipole coefficient @c TP[0] is the physical mid-plane field including its sign,
 * so the deck has to keep it consistent with the bend angle.
 *
 * The geometry, clone() and the channel table live in MultipoleTRep.
 */
class MultipoleT : public ElementBase {
public:
    /// @brief Default constructor.
    MultipoleT();

    /// @brief Constructor with given name.
    /// @param name The element name.
    explicit MultipoleT(const std::string& name);

    /// @brief Copy constructor.
    MultipoleT(const MultipoleT&);

    /// @brief Destructor.
    ~MultipoleT() override;

    /// @brief Apply a visitor.
    /// @note Resolves the scaling model first, so a clone taken by the visitor
    ///       carries it.
    /// @param visitor The visitor to apply.
    void accept(BeamlineVisitor& visitor) const override;

    /// @brief Get the element type.
    /// @return ElementType::MULTIPOLET.
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
    /// @param P Momentum (not used).
    /// @param t Time, for the scaling model.
    /// @param E Electric field (never changed; this element carries no E field).
    /// @param B Magnetic field.
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// @brief Apply the field to the reference particle with position R and momentum P.
    /// @param R Position.
    /// @param P Momentum (not used).
    /// @param t Time, for the scaling model.
    /// @param E Electric field (never changed).
    /// @param B Magnetic field.
    /// @return True if the particle is out of bounds (lost), false otherwise.
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// @brief Mark particles outside the transverse aperture in pc->InvalidMask.
    /// @note Overrides the base version to measure the z-window and the aperture
    ///       in arc-length coordinates, matching isInside().
    size_t markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /// @brief Field extent in the element's coordinate system.
    /// @note The body plus the reach of the tanh fringe on each side. A zero
    ///       fringe length is a hard edge, giving the plain body extent.
    /// @param zBegin Where the field begins (negative with an entrance fringe).
    /// @param zEnd Where the field ends (past the body with an exit fringe).
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// @brief Containment in arc-length coordinates.
    /// @note Uses the arc length, so a curved body stays selected as the
    ///       reference orbit curves through it.
    /// @param r The point to test.
    /// @return True if r is inside.
    bool isInside(const Vector_t<double, 3>& r) const override;

    /// @brief Set the transverse profile coefficients b_k [T m^-k].
    /// @note Zero-padded to the number of poles the field model supports.
    /// @param profile The coefficients, dipole first.
    void setTransverseProfile(const std::vector<double>& profile);

    /// @brief Get the transverse profile coefficients.
    /// @return The coefficients, zero-padded.
    const Kokkos::Array<double, MultipoleTFieldModel::NumPoles>& getTransverseProfile() const;

    /// @brief Set the two tanh fringe lengths [m].
    /// @note Zero is a hard edge on that side.
    /// @param left Entrance-side fringe length.
    /// @param right Exit-side fringe length.
    void setFringeLengths(double left, double right);

    /// @brief Get the entrance-side fringe length [m].
    double getFringeLeft() const;

    /// @brief Get the exit-side fringe length [m].
    double getFringeRight() const;

    /// @brief Set the number of terms in the vertical field expansion.
    /// @note Rebuilds the tanh derivative table when the order changes.
    /// @param order The number of terms, 1 to MultipoleTFieldModel::MaxFOrder.
    void setMaxFOrder(unsigned int order);

    /// @brief Get the number of terms in the vertical field expansion.
    unsigned int getMaxFOrder() const;

    /// @brief Set the name of the time dependence scaling the whole field.
    /// @note Stored upper case; an empty name means no scaling.
    /// @param name The name of the time dependence object.
    void setScalingName(const std::string& name);

    /// @brief Get the name of the time dependence scaling the whole field.
    std::string getScalingName() const;

    /// @brief Get the dipole coefficient b_0 (the "BY" channel attribute).
    double getB() const;

    /// @brief Set the dipole coefficient b_0 (the "BY" channel attribute).
    /// @param B The mid-plane dipole field [T].
    void setB(double B);

private:
    /// @brief Compute the field on the host at position R.
    /// @param R Position in the element frame.
    /// @param t Time, for the scaling model.
    /// @param B Magnetic field (accumulated).
    void computeFieldHost(const Vector_t<double, 3>& R, double t, Vector_t<double, 3>& B) const;

    /// @brief Build the pure-value field inputs shared by the device and host paths.
    /// @param t Time, for the scaling model.
    /// @return The field inputs, profile already scaled.
    MultipoleTFieldModel::FieldInputs makeFieldInputs(double t) const;

    /// @brief Convert a local point to bend coordinates (radial x, y, arc length s).
    /// @param r The local point.
    /// @return The point in bend coordinates.
    Vector_t<double, 3> bendCoords(const Vector_t<double, 3>& r) const;

    /// @brief Longitudinal containment: arc length s within the field extent.
    /// @param arc The point in bend coordinates.
    /// @return True if the arc length is within the field extent.
    bool isInsideArc(const Vector_t<double, 3>& arc) const;

    /// @brief The scaling factor of the time dependence, 1 without one.
    /// @param t Time.
    double getScaling(double t) const;

    /// @brief Look up the time dependence named by setScalingName().
    void initialiseTimeDependencies() const;

    /// @brief Size and fill the tanh derivative table for the current order.
    void rebuildTanhTable();

    /// Transverse profile coefficients b_k [T m^-k], dipole first.
    Kokkos::Array<double, MultipoleTFieldModel::NumPoles> profile_m{};

    /// Tanh derivative coefficients on the device, read by tracking.
    Kokkos::View<double**> tanhTable_m;

    /// Host mirror of the table, so the host field paths need no device copy.
    Kokkos::View<double**>::host_mirror_type tanhTableHost_m;

    /// Fringe lengths [m]; zero is a hard edge on that side.
    double fringeLeft_m;
    double fringeRight_m;

    /// Number of terms in the vertical field expansion (MAXFORDER).
    unsigned int maxFOrder_m;

    /// Name of the time dependence scaling the field; empty means no scaling.
    std::string scalingName_m;

    /// The time dependence itself, resolved from the name in accept().
    mutable std::shared_ptr<AbstractTimeDependence> scalingTD_m;
};

inline const Kokkos::Array<double, MultipoleTFieldModel::NumPoles>&
MultipoleT::getTransverseProfile() const {
    return profile_m;
}

inline double MultipoleT::getFringeLeft() const {
    return fringeLeft_m;
}

inline double MultipoleT::getFringeRight() const {
    return fringeRight_m;
}

inline unsigned int MultipoleT::getMaxFOrder() const {
    return maxFOrder_m;
}

inline std::string MultipoleT::getScalingName() const {
    return scalingName_m;
}

inline double MultipoleT::getB() const {
    return profile_m[0];
}

inline void MultipoleT::setB(double B) {
    profile_m[0] = B;
}

#endif  // ABSBEAMLINE_MULTIPOLET_H
