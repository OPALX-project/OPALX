/*
 *  Copyright (c) 2017, Titus Dascalu
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
  * $Author: Titus Dascalu, Chris Rogers
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

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "gsl/gsl_sf.h"
#include <vector>

class MultipoleT: public Component {
 public: 
    /** Constructor: &name -> user-defined name */
    explicit MultipoleT(const std::string &name);
    
    /** Copy constructor */
    MultipoleT(const MultipoleT &right);

    /** Destructor */ 
    ~MultipoleT();

    /** Inheritable copy constructor */
    ElementBase* clone() const;

    /** Return a dummy (0.) field value */
    EMField &getField();

    /** Return a dummy (0.) field value */
    const EMField &getField() const;

    /** Not implemented */
    void getDimensions(double &zBegin, double &zEnd) const;

    /** Calculate the field at some arbitrary position
     *
     *  \param R -> position in the local coordinate system of the multipole
     *  \param P -> not used
     *  \param t -> time at which the field is to be calculated
     *  \param E -> calculated electric field - always 0 (no E-field)
     *  \param B -> calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B);

    /** Calculate the field at the position of the ith particle
     *
     *  \param i -> index of the particle event; field is calculated at this
     *         position
     *  \param t -> time at which the field is to be calculated
     *  \param E -> calculated electric field - always 0 (no E-field)
     *  \param B -> calculated magnetic field
     *  \returns true if particle is outside the field map
     */
    bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

     /** Initialise the MultipoleT
      *
      *  \param -> bunch the global bunch object
      *  \param -> startField not used
      *  \param -> endField not used
      */
    void initialise(PartBunchBase<double, 3>*, double &startField, double &endField);

     /** Finalise the MultipoleT - sets bunch to NULL */
    void finalise();

     /** Return true if dipole component not zero */
    bool bends() const;

    /** Return the cell geometry */
    PlanarArcGeometry& getGeometry();

    /** Return the cell geometry */
    const PlanarArcGeometry& getGeometry() const;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const;
    
    /** Get the dipole constant B_0 */
    double getDipoleConstant() const;

    /** Set the dipole constant B_0 */
    void setDipoleConstant(double B0);
    
     /** Get the number of terms used in calculation of field components
       * max power of z in Bz is 2 * maxOrder
       */
    unsigned int getMaxOrder() const;

     /** Set the number of terms used in calculation of field components
       * max power of z in Bz is 2 * maxOrder
       */
    void setMaxOrder(unsigned int maxOrder);

    /** Get the maximum order in the given transverse profile
      * -> highest power in given mid-plane field expansion
      */
    unsigned int getTransMaxOrder() const;

    /** Set the maximum order in the given transverse profile
      * -> highest power in given mid-plane field expansion 
      */
    void setTransMaxOrder(unsigned int transMaxOrder);
    
    /** Set transverse profile T(x)
      *
      * T(x) = B_0 + B1 x + B2 x^2 + B3 x^3 + ...
      * n -> the order of the term (d^n/dx^n) to be set
      * dTn -> value of n-th derivative
      */
    void setTransProfile(unsigned int n, double Bn);

    /** Get transverse profile n-th term */
    double getTransProfile(int n) const;

    /** Get all term of transverse profile */
    std::vector<double> getTransProfile() const;

    /** Set fringe field model
      *
      * Tanh model used here
      * @f[ 1/2 * \left [tanh \left( \frac{s + s_0}{\lambda_{left}} \right) 
      * - tanh \left( \frac{s - s_0}{\lambda_{right}} \right) \right] @f] 
      * s0 -> centre field length and
      * @f$ \lambda_{left} @f$ -> left end field length
      * @f$ \lambda_{right} @f$ -> right end field length
      * max_index -> (default 10) used to set up for differentiation - cannot     
      * calculate higher differentials than exist in max_index (this is always
      * set to be 1 + twice the maxOrder_m when setting the fringe field) 
      * return true after fringe field is set
      */
    bool setFringeField(double s0, double lambda_left, double lambda_right);
    
    /** Return vector of 2 doubles
      * [left fringe length, right fringelength]  
      */
    std::vector<double> getFringeLength() const;

    /** Set the bending angle of the magnet */
    void setBendAngle(double angle);

    /** Get the bending angle of the magnet */
    double getBendAngle() const;
    
    /** Set the entrance angle */
    void setEntranceAngle(double entranceAngle);

    /** Get the entrance angle */
    double getEntranceAngle() const;

    /** Get the bending radius 
      * not used
      * when needed radius is found from length_m / angle_m
      */
    double getBendRadius() const;

    /** Set the length of the magnet
      * if straight-> actual length
      * if curved -> arc length
      */
    void setLength(double length);

    /** Get the length of the magnet */
    double getLength() const;

    /** not used */
    double getChordLength() const;

    /** Set the aperture dimensions
      * this element only supports a rectangular aperture
      */
    void setAperture(double vertAp, double horizAp);

    /** Get the aperture dimensions
      * returns a vector of 2 doubles
      */
    std::vector<double> getAperture() const;

    /** Set the angle of rotation of the magnet around its axis
      * -> enables to make skew components
      */
    void setRotation(double rot);

    /** Get the angle of rotation of the magnet around its axis */
    double getRotation() const;
    
    /** Set variable radius flag to true */
    void setVarRadius();

    /** Get the value of variableRadius_m */
    bool getVarRadius() const;

    /** Set the step used in changing coord systems for variable radius 
      * -> units used: mm
      * -> default: 50
      * -> has a considerable effect on tracking time
      */
    void setVarStep(double step);

    /** Get the step used in changing coord system for variable radius */
    double getVarStep() const;

 private:
    MultipoleT operator=(const MultipoleT& rhs);
    
    // End fields
    endfieldmodel::Tanh fringeField_l; // left
    endfieldmodel::Tanh fringeField_r; // right
    
    /** Field expansion parameters */

    // number of terms used in calculating field components
    unsigned int maxOrder_m = 0;

    // highest power in given mid-plane field
    unsigned int transMaxOrder_m = 0; 

    std::vector<double> transProfile_m;
    
    // Geometry
    PlanarArcGeometry planarArcGeometry_m;
    
    /** Rotate frame for skew elements
      * consecutive rotations:
      * 1st -> about central axis
      * 2nd -> azymuthal rotation
      */
    Vector_t rotateFrame(const Vector_t &R);

    /** Inverse of the 1st rotation in rotateFrame() method
      * -> used to rotate B field back to global coord system
      */
    Vector_t rotateFrameInverse(Vector_t &B);

    /** Transform to Frenet-Serret coordinates for sector magnets
      * -> if the radius is variable, the coord system is rotated 
      *    step by step along the reference trajectory
      * -> step size is set by setVarStep() method
      */
    Vector_t transformCoords(const Vector_t &R);

    // Step size for changing coordinates along trajectory 
    // when radius is variable
    double varStep_m;

    // Magnet parameters
    double length_m;
    double angle_m;
    double entranceAngle_m;
    double rotation_m;

    // Variable radius flag;
    bool variableRadius_m;

    // Get field component methods
    double getBx (const Vector_t &R);
    double getBz (const Vector_t &R);
    double getBs (const Vector_t &R);

    // Assume rectangular aperture with dimensions given by
    double verticalApert_m;
    double horizApert_m;

    // Not implemented
    BMultipoleField dummy;
    
    // Returns the value of the fringe field n-th derivative at s
    double getFringeDeriv(int n, double s);

    // Returns the value of the transverse field n-th derivative at x
    double getTransDeriv(unsigned int n, double x);
    
    // Tests if inside the magnet
    bool insideAperture(const Vector_t &R);

    double getRadius(double s);
    /** The radius is calculated to be proportional to the field on the 
      * refrence trajectory (central axis of magnet) 
      * @f$ \rho (s) = \rho(0) * S(0) / S(s) @f$
      * -> returns -1 if radius is infinite
      */

    double getRadiusFirstDeriv(double s);
    /** @return the derivative of the function which describes the radius
      * as a function of s; 
      */

    double getRadiusSecDeriv(double s);
    /** @return The second derivative of radius function
      */

    double getScaleFactor(double x, double s);
    /** Returns the scale factor @f$ h_s = [(1 + x / \rho)^2 + 
     *  (d \rho / ds)^2]^0.5 @f$
      * rho -> radius
      * used in field expansion, comes from Laplacian when curvature is variable
      * if radius is fixed returns @f$  h_s = 1 + x / \rho @f$
      */

    double getScaleFactorDerivX(double x, double s);
    /** Helper function to bunch terms together
      * -> returns the partial deriv of the scale factor wrt. x
      * -> @f$ \partial_x h_s = (1 + x / \rho) / (\rho * h_s ^ 0.5) @f$
      */

    double getScaleFactorDerivS(double x, double s);
    /** Helper function to bunch terms together
      * -> returns @f$ \partial_s h_s @f$
      * -> @f$ = (h_s)^(-1/2) * \partial_s \rho * 
      * [-(1+x/\rho)*x/\rho^2 + {\partial_x}^2 \rho] @f$
      */ 
    
    double getFnDerivX(unsigned int n, double x, double s);
    /** Returns the partial derivative of f_n(x, s) wrt x
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */
    
    double getFnDerivS(unsigned int n, double x, double s);
    /** Returns the partial derivative of f_n(x, s) wrt s
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */

    double getFnSecDerivX(unsigned int n, double x, double s);
    /** Returns the second partial derivative of f_n(x, s) wrt x
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */

    double getFnSecDerivS(unsigned int n, double x, double s);
    /** Returns the second partial derivative of f_n(x, s) wrt s
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */

    double getFn(unsigned int n, double x, double s);
    /** Returns the function @f$ f_n(x, s) @f$ for curved geometry
      * -> based on recursion
      */
    
};
inline
    void MultipoleT::setVarStep(double step) {
    varStep_m = step;
}
inline
    double MultipoleT::getVarStep() const {
        return varStep_m;
}
inline
    void MultipoleT::setVarRadius() {
        variableRadius_m = true;
}
inline
    bool MultipoleT::getVarRadius() const {
        return variableRadius_m;
}
inline
    void MultipoleT::setEntranceAngle(double entranceAngle) {
        entranceAngle_m = entranceAngle;
}
inline
    double MultipoleT::getEntranceAngle() const {
        return entranceAngle_m;
}
inline 
    double MultipoleT::getTransProfile(int n) const {
        return transProfile_m[n];
}
inline
    std::vector<double> MultipoleT::getTransProfile() const {
        return transProfile_m;
}
inline 
    double MultipoleT::getDipoleConstant() const {
         return transProfile_m[0];
}
inline
    unsigned int MultipoleT::getMaxOrder() const {
         return maxOrder_m;
}
inline
    void MultipoleT::setMaxOrder(unsigned int maxOrder) {
         maxOrder_m = maxOrder;
}
inline
    unsigned int MultipoleT::getTransMaxOrder() const {
        return transMaxOrder_m;
}
inline 
    void MultipoleT::setTransMaxOrder(unsigned int transMaxOrder) {
        transMaxOrder_m = transMaxOrder;
	transProfile_m.resize(transMaxOrder + 1, 0.);
}
inline
    double MultipoleT::getRotation() const {
         return rotation_m;
}
inline 
    void MultipoleT::setRotation(double rot) {
         rotation_m = rot;
}
inline
    void MultipoleT::setBendAngle(double angle) {
        angle_m = angle;
}
inline
    double MultipoleT::getBendAngle() const {
        return angle_m;
}
inline
    void MultipoleT::setLength(double length) {
        length_m = std::abs(length);
}
inline
    double MultipoleT::getLength() const {
        return length_m;
}

#endif
