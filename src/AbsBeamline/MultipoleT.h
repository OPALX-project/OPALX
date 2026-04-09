//
// Cubic Spline Interpolation to replace GSL spline
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
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

/** ---------------------------------------------------------------------
 *
 * MultipoleT defines a straight or curved combined function magnet (up
 * to arbitrary multipole component) with Fringe Fields
 *
 * ---------------------------------------------------------------------
 *
 * Class category: AbsBeamline \n
 * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
 *
 * ---------------------------------------------------------------------
 *
 * The field is obtained from the scalar potential \n
 *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!} +
 *     ...  @f] \n
 *     (x,z,s) -> Frenet-Serret local coordinates along the magnet \n
 *     z -> vertical component \n
 *     assume mid-plane symmetry \n
 *     set field on mid-plane -> @f$ B_z = f_0(x,s) = T(x) \cdot S(s) @f$ \n
 *     T(x) -> transverse profile; this is a polynomial describing
 *             the field expansion on the mid-plane inside the magnet
 *             (not in the fringe field);
 *             1st term is the dipole strength, 2nd term is the
 *             quadrupole gradient * x, etc. \n
 *          -> when setting the magnet, one gives the multipole
 *             coefficients of this polynomial (i.e. dipole strength,
 *             quadrupole gradient, etc.) \n
 * \n
 * ------------- example ----------------------------------------------- \n
 *     Setting a combined function magnet with dipole, quadrupole and
 *     sextupole components: \n
 *     @f$ T(x) = B_0 + B_1 \cdot x + B_2 \cdot x^2 @f$\n
 *     user gives @f$ B_0, B_1, B_2 @f$ \n
 * ------------- example end ------------------------------------------- \n
 * \n
 *     S(s) -> fringe field \n
 *     recursion -> @f$ f_n (x,s) = (-1)^n \cdot \sum_{i=0}^{n} C_n^i
 *     \cdot T^{(2i)} \cdot S^{(2n-2i)} @f$ \n
 *     for curved magnets the above recursion is more complicated \n
 *     @f$ C_n^i @f$ -> binomial coeff;
 *     @f$ T^{(n)} @f$ -> n-th derivative
 *
 * ---------------------------------------------------------------------
 */
#include <vector>
#include "AbsBeamline/Component.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Fields/BMultipoleField.h"
#include "MultipoleTBase.h"
#include "MultipoleTConfig.h"

class MultipoleT : public Component {
public:
    /** Constructor
     *  \param name -> User-defined name
     */
    explicit MultipoleT(const std::string& name);
    /** Copy constructor */
    MultipoleT(const MultipoleT& right);
    /** Destructor */
    ~MultipoleT() override = default;
    /** Inheritable copy constructor */
    ElementBase* clone() const override;
    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;
    /** Return the cell geometry */
    BGeometryBase& getGeometry() override;
    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const override;
    /** Return a dummy field value */
    EMField& getField() override { return dummy; }
    /** Return a dummy field value */
    const EMField& getField() const override { return dummy; }
    /** Calculate the field for all particles */
    bool apply() override;
    /** Calculate the field at some arbitrary position \n
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param R -> Position in the lab coordinate system of the multipole
     *  \param P -> Not used
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    /** Calculate the field at the position of the ith particle
     *  \param i -> Index of the particle event; field is calculated at this
     *  position
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B)
            override;
    /** Initialise the MultipoleT
     *  \param bunch -> Bunch the global bunch object
     *  \param startField -> Not used
     *  \param endField -> Not used
     */
    void initialise(PartBunch_t* bunch, double& startField, double& endField) override;
    /** Finalise the MultipoleT - sets bunch to nullptr */
    void finalise() override;
    /** Return true if dipole component not zero */
    bool bends() const override;
    /** Get the number of terms used in calculation of field components */
    size_t getMaxFOrder() const { return config_m.maxFOrder_m; }
    size_t getMaxXOrder() const { return config_m.maxXOrder_m; }
    /** Set the number of terms used in calculation of field components \n
     *  \param orderZ -> Number of terms in expansion in z
     *  \param orderX -> Number of terms in expansion in x
     */
    void setMaxOrder(size_t orderZ, size_t orderX);
    /** Get the maximum order in the given transverse profile */
    std::size_t getTransMaxOrder() const { return config_m.transverseProfileMaxOrder_m; }
    /** Set the the transverse profile
     *  \param profile -> Multipole field profile
     */
    void setTransProfile(const std::vector<double>& profile);
    /** Get all terms of transverse profile */
    const Kokkos::Array<double, MultipoleTConfig::NumPoles>& getTransProfile() const {
        return config_m.transverseProfile_m;
    }

    /** Set fringe field model \n
     *  Tanh model used here \n
     *  @f[ 1/2 * \left [tanh \left( \frac{s + s_0}{\lambda_{left}} \right)
     *  - tanh \left( \frac{s - s_0}{\lambda_{right}} \right) \right] @f]
     *  \param s0 -> Centre field length and
     *  \param lambda_left -> Left end field length
     *  \param lambda_right -> Right end field length
     */
    void setFringeField(const double& s0, const double& lambda_left, const double& lambda_right);
    /** Get the fringe field model
     * @return {s0, leftFringe, rightFringe}
     */
    std::tuple<double, double, double> getFringeField() const;
    /** Set the entrance angle
     *  \param entranceAngle -> Entrance angle
     */
    void setEntranceAngle(double entranceAngle);
    /** Set the offset of the entry point from the standard position
     * @param offset positive for further away from the center
     */
    void setEntryOffset(double offset);
    /** Get the offset of the entry point from the standard position */
    double getEntryOffset() const { return config_m.entryOffset_m; }
    /** Get the variable radius of the magnet */
    bool getVariableRadius() const { return config_m.variableRadius_m; }
    /** Set the bending angle of the magnet */
    void setBendAngle(double angle, bool variableRadius);
    /** Get the bending angle of the magnet */
    double getBendAngle() const { return config_m.bendAngle_m; }
    /** Get the entrance angle */
    double getEntranceAngle() const { return config_m.entranceAngle_m; }
    /** Set the length of the magnet
     * If straight-> Actual length
     * If curved -> Arc length
     */
    void setElementLength(double length) override;
    /** Get the length of the magnet */
    double getLength() const { return config_m.length_m; }
    /** Set the aperture dimensions \n
     * This element only supports a rectangular aperture
     * \param vertAp -> Vertical aperture length
     * \param horizAp -> Horisontal aperture length
     */
    void setAperture(const double& vertAp, const double& horizAp);
    /** Get the aperture dimensions
     * @return {verticalApert, horizontalApert}
     */
    std::tuple<double, double> getAperture() {
        return {config_m.verticalAperture_m, config_m.horizontalAperture_m};
    }

    /** Set the angle of rotation of the magnet around its axis \n
     *  To make skew components
     *  \param rot -> Angle of rotation
     */
    void setRotation(double rot);
    /** Get the angle of rotation of the magnet around its axis */
    double getRotation() const { return config_m.rotation_m; }
    /** Get the bounding box size */
    double getBoundingBoxLength() const { return config_m.boundingBoxLength_m; }
    /** Set the bounding box size.  This controls the region for
     *  \param boundingBoxLength -> Distance between centre and entrance
     */
    void setBoundingBoxLength(double boundingBoxLength);
    /** Not implemented */
    void getDimensions(double& /*zBegin*/, double& /*zEnd*/) const override {}

    void setScalingName(const std::string& name);
    std::string getScalingName() const { return scalingName_m; }
    void initialiseTimeDependencies() const;

    MultipoleTConfig& getConfig() { return config_m; }

protected:
    /** The magnet configuration */
    MultipoleTConfig config_m;

    void chooseImplementation();
    double getScaling(double t) const;
    void validateConfiguration() const;

    /** Not implemented */
    BMultipoleField dummy;

    // Time dependence
    std::string scalingName_m;
    mutable std::shared_ptr<AbstractTimeDependence> scalingTD_m;

    // The object that does the work
    std::unique_ptr<MultipoleTBase> implementation_{};
};
#endif
