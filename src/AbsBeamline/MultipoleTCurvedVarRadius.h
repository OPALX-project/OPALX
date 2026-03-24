/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
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

#ifndef CLASSIC_MULTIPOLET_CURVED_VAR_RADIUS_H
#define CLASSIC_MULTIPOLET_CURVED_VAR_RADIUS_H

/** ---------------------------------------------------------------------
 *
 * MultipoleTCurvedVarRadius defines a curved combined function magnet with\n
 * variable radius of curvature (up to arbitrary multipole component) \n
 * with fringe fields
 *
 * ---------------------------------------------------------------------
 *
 * Class category: AbsBeamline \n
 * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
 *
 * ---------------------------------------------------------------------
 *
 * The field is obtained from the scalar potential \n
 *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!}
 *      + ...  @f] \n
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

#include "BeamlineGeometry/VarRadiusGeometry.h"
#include <vector>
#include "AbsBeamline/MultipoleTBase.h"
#include "AbsBeamline/MultipoleTFunctions/RecursionRelationTwo.h"
#include "Tanh.h"

class MultipoleTCurvedVarRadius : public MultipoleTBase {
public:
    /** Constructor */
    explicit MultipoleTCurvedVarRadius(MultipoleT* element);
    /** Initialise the element */
    void initialise() override;
    /** Return the cell geometry */
    BGeometryBase* getGeometry() override { return &varRadiusGeometry_m; }
    /** Return the cell geometry */
    const BGeometryBase* getGeometry() const override { return &varRadiusGeometry_m; }
    /** Return the field for an array of points */
    void getField(
        const Kokkos::View<Vector_t<double, 3>*>& R, Kokkos::View<Vector_t<double, 3>*>& E,
        Kokkos::View<Vector_t<double, 3>*>& B, double scaling) override;
    bool getField(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B,
        double scaling) override;
    /** Transform B-field from Frenet-Serret coordinates to lab coordinates */
    void transformBField(Vector_t<double, 3>& /*B*/, const Vector_t<double, 3>& /*R*/) override;
    /** Returns the scale factor @f$ h_s = 1@f$
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getScaleFactor(double x, double s) override;
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFn(unsigned int n, double x, double s) override;
    /** Set the number of terms used in calculation of field components \n
     *  Maximum power of z in Bz is 2 * maxOrder_m
     *  \param orderZ -> Number of terms in expansion in z
     *  \param orderX -> Number of terms in expansion in x
     */
    void setMaxOrder(size_t orderZ, size_t orderX) override;

    Vector_t<double, 3> localCartesianToCurvilinear(const Vector_t<double, 3>& r);
    Vector_t<double, 3> curvilinearToLocalCartesian(const Vector_t<double, 3>& r);
    Vector_t<double, 3> localCartesianToOpalCartesian(const Vector_t<double, 3>& r) override;
    double reverseTransformResidual(
        const Vector_t<double, 3>& r, const Vector_t<double, 3>& target);
    static constexpr size_t ReverseTransformMaxIterations = 1000;
    static constexpr double ReverseTransformTolerance     = 1e-6;
    static constexpr double TangentStep                   = 1e-3;
    double localCartesianRotation() override { return localCartesianRotation_; }

private:
    /** Objects for storing differential operator acting on Fn */
    std::vector<polynomial::RecursionRelationTwo> recursion_m;
    /** Geometry */
    VarRadiusGeometry varRadiusGeometry_m;
    Vector_t<double, 3> localCartesianEntryPoint_;
    double localCartesianRotation_{};
    // End fields
    endfieldmodel::Tanh fringeField_l{endfieldmodel::Tanh()};  // Left
    endfieldmodel::Tanh fringeField_r{endfieldmodel::Tanh()};  // Right

    // Helper functions
    Vector_t<double, 3> toMagnetCoords(
        const Vector_t<double, 3>& R, const MultipoleTConfig& config);
    Vector_t<double, 3> rotateFrame(
        const Vector_t<double, 3>& R, const MultipoleTConfig& config) const;
    /** Transform to Frenet-Serret coordinates for sector magnets */
    void transformCoords(Vector_t<double, 3>& /*R*/) override;
    bool computeBField(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& B, double scaling,
        const MultipoleTConfig& config);
    /** Returns the radial component of the field \n
     *  Returns zero far outside fringe field
     *  @f$ Bx = sum_n z^(2n+1) / (2n+1)! * \partial_x f_n @f$
     */
    double getBx(const Vector_t<double, 3>& R, const MultipoleTConfig& config);
    /** Returns the vertical field component \n
     *  Returns zero far outside fringe field
     *  @f$ Bz = sum_n  f_n * z^(2n) / (2n)! @f$
     */
    double getBz(const Vector_t<double, 3>& R, const MultipoleTConfig& config);
    /** Returns the component of the field along the central axis \n
     *  Returns zero far outside fringe field
     * @f$ Bs = sum_n z^(2n+1) / (2n+1)! \partial_s f_n / h_s @f$
     */
    double getBs(const Vector_t<double, 3>& R, const MultipoleTConfig& config);
    /** Calculate partial derivative of fn wrt x using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivX(
        const std::size_t& n, const double& x, const double& s, const MultipoleTConfig& config);
    /** Returns the value of the transverse field n-th derivative at x \n
     *  Transverse field is a polynomial in x, differentiation follows
     *  usual polynomial rules of differentiation
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     */
    double getTransDeriv(
        const std::size_t& n, const double& x, const MultipoleTConfig& config) const;
    /** Returns the value of the fringe field n-th derivative at s
     *  \param n -> nth derivative
     *  \param s -> Coordinate s
     */
    double getFringeDeriv(const std::size_t& n, const double& s);
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFn(size_t n, double x, double s);
    /** Calculate partial derivative of fn wrt s using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivS(
        const std::size_t& n, const double& x, const double& s, const MultipoleTConfig& config);
};

#endif
