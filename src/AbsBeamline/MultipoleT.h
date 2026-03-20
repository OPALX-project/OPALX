/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
*  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CLASSIC_MULTIPOLET_H
#define CLASSIC_MULTIPOLET_H

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
    bool apply(const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
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
    bool apply(const size_t& i, const double& t, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) override;
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
    void setFringeField(const double& s0, const double& lambda_left,
            const double& lambda_right);
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
    void getDimensions(double&/*zBegin*/, double&/*zEnd*/) const override {}
    /** Returns the value of the fringe field n-th derivative at s
     *  \param n -> nth derivative
     *  \param s -> Coordinate s
     */
    double getFringeDeriv(const std::size_t& n, const double& s);
    /** Calculate partial derivative of fn wrt x using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivX(const std::size_t& n,
            const double& x,
            const double& s);
    /** Calculate partial derivative of fn wrt s using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivS(const std::size_t& n,
            const double& x,
            const double& s);
    /** Returns the value of the transverse field n-th derivative at x \n
     *  Transverse field is a polynomial in x, differentiation follows
     *  usual polynomial rules of differentiation
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     */
    double getTransDeriv(const std::size_t& n, const double& x) const;

    Vector_t<double, 3> toMagnetCoords(const Vector_t<double, 3>& R);
    Vector_t<double, 3> getField(const Vector_t<double, 3>& magnetCoords);
    Vector_t<double, 3> localCartesianToOpalCartesian(const Vector_t<double, 3>& r);
    double localCartesianRotation();

    void setScalingName(const std::string& name);
    void setScalingModel(const std::shared_ptr<AbstractTimeDependence>& td) { scalingTD_m = td; }
    std::string getScalingName() const { return scalingName_m; }
    void initialiseTimeDepencencies() const;

    MultipoleTConfig& getConfig() { return config_m; }

protected:
    /** The magnet configuration */
    MultipoleTConfig config_m;

    /** Rotate the frame to account for the rotation and entry angles.
     * @param R -> coordinate to rotate
     * @return -> rotated coordinate
     */
    Vector_t<double, 3> rotateFrame(const Vector_t<double, 3>& R) const;
    bool insideAperture(const Vector_t<double, 3>& R) const;
    bool insideBoundingBox(const Vector_t<double, 3>& R) const;
    void chooseImplementation();
    double getScaling(double t) const;

    /** Not implemented */
    BMultipoleField dummy;

    // Time dependence
    std::string scalingName_m;
    mutable std::shared_ptr<AbstractTimeDependence> scalingTD_m;

    // The object that does the work
    std::unique_ptr<MultipoleTBase> implementation_{};

    /** This one is here for test purposes.
     * It is required as it is too difficult to mock a PartBunch. */
    void apply(const Kokkos::View<Vector_t<double, 3>*>& R, Kokkos::View<Vector_t<double, 3>*>& E,
            Kokkos::View<Vector_t<double, 3>*>& B, double t) const;
};
#endif
