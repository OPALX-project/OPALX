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


#include "Algorithms/PartBunch.h"
#include "MultipoleT.h"

using namespace endfieldmodel;

MultipoleT::MultipoleT(const std::string &name):
    Component(name),   
    fringeField_l(endfieldmodel::Tanh()),
    fringeField_r(endfieldmodel::Tanh()),
    maxOrder_m(5),
    transMaxOrder_m(0),
    planarArcGeometry_m(1., 1.),
    varStep_m(0.1),
    length_m(1.0),
    angle_m(0.0),
    entranceAngle_m(0.0),
    rotation_m(0.0),
    variableRadius_m(false),
    verticalApert_m(0.5),
    horizApert_m(0.5) {
}

MultipoleT::MultipoleT(const MultipoleT &right):
    Component(right),
    fringeField_l(right.fringeField_l),
    fringeField_r(right.fringeField_r),
    maxOrder_m(right.maxOrder_m),
    transMaxOrder_m(right.transMaxOrder_m),
    transProfile_m(right.transProfile_m),
    planarArcGeometry_m(right.planarArcGeometry_m),
    varStep_m(right.varStep_m),
    length_m(right.length_m),
    angle_m(right.angle_m),
    entranceAngle_m(right.entranceAngle_m),
    rotation_m(right.rotation_m),
    variableRadius_m(right.variableRadius_m),
    verticalApert_m(right.verticalApert_m),
    horizApert_m(right.horizApert_m),
    dummy() {

    RefPartBunch_m = right.RefPartBunch_m;
}
    

MultipoleT::~MultipoleT() {
} 

ElementBase* MultipoleT::clone() const {
    return new MultipoleT(*this);
}

void MultipoleT::finalise() {
    RefPartBunch_m = NULL;
}

bool MultipoleT::apply(const Vector_t &R, const Vector_t &P,
               const double &t,Vector_t &E, Vector_t &B) {
    /** Rotate coordinates around the central axis of the magnet */
    Vector_t R_prime = rotateFrame(R);
    /** If magnet is not straight go to coords along the magnet */
    if(angle_m != 0.0) {
        Vector_t X = R_prime;
    R_prime = transformCoords(X);
    }
    if (insideAperture(R_prime)) {
        B[0] = getBx(R_prime);
    B[1] = getBz(R_prime);
    B[2] = getBs(R_prime);
    return false;
    }
    else {
        for(int i = 0; i < 3; i++) {
        B[i] = 0.0;
        }
        return true;
    }
}

bool MultipoleT::apply(const size_t &i, const double &t,
               Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

Vector_t MultipoleT::rotateFrame(const Vector_t &R) {
  Vector_t R_prime(3), R_pprime(3);
    /** Apply two 2D rotation matrices to coords vector
      * -> rotate around central axis => skew fields
      * -> rotate azymuthaly => entrance angle
      */
    // 1st rotation
    R_prime[0] = R[0] * cos(rotation_m) + R[1] * sin(rotation_m);
    R_prime[1] = - R[0] * sin(rotation_m) + R[1] * cos(rotation_m);
    R_prime[2] = R[2];
    // 2nd rotation
    R_pprime[0] = R_prime[2] * sin(entranceAngle_m) + 
                  R_prime[0] * cos(entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * cos(entranceAngle_m) - 
                  R_prime[0] * sin(entranceAngle_m);
    return R_pprime;

}

Vector_t MultipoleT::rotateFrameInverse(Vector_t &B) {
    /** This function represents the inverse of the rotation 
      * around the central axis performed by rotateFrame() method
      * -> used to rotate B field back to global coord system
      */
    Vector_t B_prime(3);
    B_prime[0] = B[0] * cos(rotation_m) + B[1] * sin(rotation_m);
    B_prime[1] = - B[0] * sin(rotation_m) + B[1] * cos(rotation_m);
    B_prime[2] = B[2];
    return B_prime;
}

bool MultipoleT::insideAperture(const Vector_t &R) {
    if (abs(R[1]) <= verticalApert_m / 2. && abs(R[0]) <= horizApert_m / 2.) {
        return true;
    }
    else {
        return false;
    }
}

Vector_t MultipoleT::transformCoords(const Vector_t &R) {
    Vector_t X(3), Y(3);
    // if radius not variable
    if (not variableRadius_m) {
        double radius = length_m / angle_m;
    // Transform from Cartesian to Frenet-Seret along the magnet
    double alpha = atan(R[2] / (R[0] + radius ));
    if (alpha != 0.0) {
            X[0] = R[2] / sin(alpha) - radius;
        X[1] = R[1];
        X[2] = radius * alpha;
    }
    else {
            X = R;
    }
    }
    // if radius is variable need to transform coordinates at each
    // point along the trajectory
    else {
        double deltaAlpha, S = 0.0, localRadius;
    double stepSize = varStep_m; // mm -> has a big effect on tracking time
    if (abs(R[2]) <= stepSize) {
        return R; // no transformation around origin
    }
    Y = R;
    Vector_t temp;
        while (abs(Y[2]) > stepSize and getRadius(S) != -1) {
        localRadius = getRadius(S);
        deltaAlpha = stepSize / localRadius;
        if (R[2] < 0) {
            deltaAlpha *= - 1.; // rotate in the other direction
        }
            temp = Y;
        // translation
        temp[0] += localRadius * (1 - cos(deltaAlpha));
        temp[2] -= localRadius * sin(deltaAlpha);
        // + rotation along the ideal trajectory
        Y[2] = temp[2] * cos(deltaAlpha) - temp[0] * sin(deltaAlpha);
        Y[0] = temp[2] * sin(deltaAlpha) + temp[0] * cos(deltaAlpha);  
        S += localRadius * deltaAlpha;
        // until we reach the actual point from Cartesian coordinates
    }
    X[0] = Y[0];
    X[1] = Y[1];
    X[2] = S;
    }     
    return X;
}

void MultipoleT::setTransProfile(unsigned int n, double dTn) {
    if (n > transMaxOrder_m) {
      transMaxOrder_m = n;
      transProfile_m.resize(n+1, 0.0);
    }
    transProfile_m[n] = dTn;
 }

bool MultipoleT::setFringeField(double s0, double lambda_l, double lambda_r) {
    fringeField_l.Tanh::setLambda(lambda_l);
    fringeField_r.Tanh::setLambda(lambda_r);
    fringeField_l.Tanh::setX0(s0);
    fringeField_l.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);
    fringeField_r.Tanh::setX0(s0);
    fringeField_r.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);
    return true;
}

double MultipoleT::getBz(const Vector_t &R) {
    /** Returns the vertical field component
      * -> sum_n  f_n * z^(2n) / (2n)!
      */
    double Bz = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> use corresponding field expansion 
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
        double f_n = 0.0;
        for(unsigned int i = 0; i <= n; i++)  {
            f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, R[0]) * 
                   getFringeDeriv(2 * n - 2 * i, R[2]);
        }
        f_n *= gsl_sf_pow_int(-1.0, n);
        Bz += gsl_sf_pow_int(R[1], 2 * n) / gsl_sf_fact(2 * n) * f_n;
    }
    }
    else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // return 0 if end of fringe field is reached
            // this is to avoid functions being called at infinite radius
            return 0.0; 
    }
        // Curved geometry -> use corresponding field expansion
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
            double f_n = getFn(n, R[0], R[2]);
        Bz += gsl_sf_pow_int(R[1], 2 * n) / gsl_sf_fact(2 * n) * f_n;
      }
    }
    return Bz;
}

double MultipoleT::getBx(const Vector_t &R) {
    /** Returns the radial component of the field
      * -> sum_n z^(2n+1) / (2n+1)! * \partial_x f_n 
      */
    double Bx = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> use corresponding field expansion
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
        double f_n = 0.0;
        for(unsigned int i = 0; i <= n; i++) {
            f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i + 1, R[0]) * 
               getFringeDeriv(2 * n - 2 * i, R[2]);
        }
        f_n *= gsl_sf_pow_int(-1.0, n);
        Bx += gsl_sf_pow_int(R[1], 2 * n + 1) / 
              gsl_sf_fact(2 * n + 1) * f_n;
    }
    }
    else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // return 0 if end of fringe field is reached
            // this is to avoid functions being called at infinite radius
            return 0.0; 
    }
        // Curved geometry -> use corresponding field expansion
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
            double partialX_fn = getFnDerivX(n, R[0], R[2]);
        Bx += gsl_sf_pow_int(R[1], 2 * n + 1) / gsl_sf_fact(2 * n + 1); 
        Bx *= partialX_fn;
    }
    }
    return Bx;
} 

double MultipoleT::getBs(const Vector_t &R) {
    /** Returns the component of the field along the central axis 
      * -> 1/h_s * sum_n z^(2n+1) / (2n+1)! \partial_s f_n
      */
    double Bs = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> use corresponding field expansion 
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
        double f_n = 0.0;
        for(unsigned int i = 0; i <= n; i++) {
            f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, R[0]) *
                   getFringeDeriv(2 * n - 2 * i + 1, R[2]);
        }
        f_n *= gsl_sf_pow_int(-1.0, n);
        Bs += gsl_sf_pow_int(R[1], 2 * n + 1) / 
                  gsl_sf_fact(2 * n + 1) * f_n;
    }
    }
    else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // return 0 if end of fringe field is reached
            // this is to avoid functions being called at infinite radius
            return 0.0; 
    }
        // Curved geometry -> use corresponding field expansion
        for(unsigned int n = 0; n <= maxOrder_m; n++) {
        double partialS_fn = getFnDerivS(n, R[0], R[2]);
        Bs += gsl_sf_pow_int(R[1], 2 * n + 1) / gsl_sf_fact(2 * n + 1);
        Bs *= partialS_fn;
    }
    Bs /= getScaleFactor(R[0], R[2]);
    }
    return Bs;
}

double MultipoleT::getFringeDeriv(int n, double s) {
    // Wraps around the getTanh method of endfield 
    double temp;
    temp = fringeField_l.getTanh(s, n) - fringeField_r.getNegTanh(s, n);
    temp /= 2;
    return temp;
}

double MultipoleT::getTransDeriv(unsigned int n, double x) {
    /** Sets a vector of the coefficients in the polynomial expansion
      * of transverse profile; shifts them to the left and multiply by
      * corresponding power each time to take derivative once; 
      * repeats until desired derivative is reached
      */   
    double func = 0;
    std::vector<double> temp = transProfile_m;
    if (n <= transMaxOrder_m) {
        if (n != 0) {
        for(unsigned int i = 1; i <= n; i++) {
            for(unsigned int j = 0; j <= transMaxOrder_m; j++) {
            if (j <= transMaxOrder_m - i) {
                // move terms to the left and multiply by power
                temp[j] = temp[j + 1] * (j + 1);
                }
            else {
                    // put 0 at the end for missing higher powers
                temp[j] = 0.0;
            }
        }
        }
    }
    // Now use the vector to calculate value of the function
        for(unsigned int k = 0; k <= transMaxOrder_m; k++) {
            func += temp[k] * gsl_sf_pow_int(x, k);
    }
    }
    return func;
}

void MultipoleT::setDipoleConstant(double B0) {
    if (transMaxOrder_m < 1) {
        transProfile_m.resize(1, 0.);
    }
    transProfile_m[0] = B0;
}

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    visitor.visitMultipoleT(*this);
}

void MultipoleT::getDimensions(double &zBegin, double &zEnd) const {
}

void MultipoleT::setAperture(double vertAp, double horizAp) {
    verticalApert_m = vertAp;
    horizApert_m = horizAp;
}

std::vector<double> MultipoleT::getAperture() const {
    std::vector<double> temp(2, 0.0);   
    temp[0] = verticalApert_m;
    temp[1] = horizApert_m;
    return temp;
}

std::vector<double> MultipoleT::getFringeLength() const {
    std::vector<double> temp(2, 0.0);
    temp[0] = fringeField_l.getLambda();
    temp[1] = fringeField_r.getLambda();
    return temp;
}

double MultipoleT::getRadius(double s) {
    /** The radius is calculated to be inverse proportional to the field on
      * the refrence trajectory (central axis of magnet)
      * @f$ \rho (s) = \rho(0) * S(0) / S(s) @f$
      */
    double centralRadius = length_m / angle_m;
    // at the centre of the magnet
    double propCoeff = centralRadius * getFringeDeriv(0, 0);
    // move to current position on central axis
    if (getFringeDeriv(0, s) != 0.0) {
       return propCoeff / getFringeDeriv(0, s);
    }
    else {
       return -1; // return -1 if radius is infinite 
    }
}

double MultipoleT::getRadiusFirstDeriv(double s) {
    /** @return The derivative of the function which describes the radius
      * as a function of s; 
      */
    return  -1.0 * getRadius(s) * getFringeDeriv(1, s) / getFringeDeriv(0, s);
}

double MultipoleT::getRadiusSecDeriv(double s) {
    /** @return The second derivative of radius function
      */
    double deriv = 0.0;
    deriv += gsl_sf_pow_int(getFringeDeriv(1, s), 2.0) / getFringeDeriv(0, s);
    deriv -= getFringeDeriv(2, s);
    deriv *= getRadius(s);
    return deriv;
}

double MultipoleT::getScaleFactor(double x, double s) {
    /** @return The scale factor h_s = [(1 + x / rho)^2 + (d rho / ds)^2]^0.5
      * rho -> radius
      * used in field expansion, comes from Laplacian when curvature is variable
      * if radius is fixed returns h_s = 1 + x / rho
      */
    if (not variableRadius_m) {
        double rho = length_m / angle_m;
        return 1 + x / rho;
    }
    else {
        double temp = 0.0;
    temp += gsl_sf_pow_int(getRadiusFirstDeriv(s), 2.0);
    temp += gsl_sf_pow_int(1 + x / getRadius(s), 2.0);
    return gsl_sf_pow_int(temp, 0.5);
    }
}

double MultipoleT::getScaleFactorDerivX(double x, double s) {
    /** Helper function to bunch terms together
      * -> returns the partial deriv of the scale factor wrt. x
      * -> partial h_s / partial x#include "Algorithms/PartBunch.h"
      * -> = (1 + x / rho) / (rho * h_s ^ 0.5)
      */
    double temp = 0.0;
    temp += (1 + x / getRadius(s)) / getRadius(s);
    temp /= gsl_sf_pow_int(getScaleFactor(x, s), 0.5);
    return temp;
}

double MultipoleT::getScaleFactorDerivS(double x, double s) {
    /** Helper function to bunch terms together
      * -> returns partial h_s / partial s
      * -> = (h_s)^(-1/2) * partial rho / partial s *
      * -> * [-(1+x/rho)*x/rho^2 + partial^2 rho / partial^x]
      */ 
    double temp = 0.0;
    temp -= (1 + x /getRadius(s)) * x / gsl_sf_pow_int(getRadius(s), 2.);
    temp += getRadiusSecDeriv(s);
    temp *= getRadiusFirstDeriv(s) / gsl_sf_pow_int(getScaleFactor(x, s), -0.5);
    return temp;
}

double MultipoleT::getFnDerivX(unsigned int n, double x, double s) {
    /** Returns the partial derivative of f_n(x, s) wrt x
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */
    if (n == 0) {
      return getTransDeriv(1, x) * getFringeDeriv(0, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += getFn(n, x - 2. * stepSize, s) - 8. * getFn(n, x - stepSize, s);
    deriv += 8. * getFn(n, x + stepSize, s) - getFn(n, x + 2. * stepSize, s);
    deriv /= 12.0 * stepSize;
    return deriv;
}

double MultipoleT::getFnDerivS(unsigned int n, double x, double s) {
    /** Returns the partial derivative of f_n(x, s) wrt s
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */
    if (n == 0) {
      return getTransDeriv(0, x) * getFringeDeriv(1, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += getFn(n, x, s - 2. * stepSize) - 8. * getFn(n, x, s - stepSize);
    deriv += 8. * getFn(n, x, s + stepSize) - getFn(n, x, s + 2. * stepSize);
    deriv /= 12.0 * stepSize;
    return deriv;
}

double MultipoleT::getFnSecDerivX(unsigned int n, double x, double s) {
    /** Returns the second partial derivative of f_n(x, s) wrt x
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */
    if (n == 0) {
      return getTransDeriv(2, x) * getFringeDeriv(0, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += getFn(n, x + 2. * stepSize, s) + getFn(n, x - 2. * stepSize, s);
    deriv -= getFn(n, x + stepSize, s) + getFn(n, x - stepSize, s);
    deriv /= 3. * gsl_sf_pow_int(stepSize, 2.);
    return deriv;
}

double MultipoleT::getFnSecDerivS(unsigned int n, double x, double s) {
    /** Returns the second partial derivative of f_n(x, s) wrt s
      * -> numerical differentiation
      * -> 5-points formula
      * -> error of order stepSize ^ 4
      */
    if (n == 0) {
      return getTransDeriv(0, x) * getFringeDeriv(2, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += getFn(n, x, s + 2. * stepSize) + getFn(n, x, s - 2. * stepSize);
    deriv -= getFn(n, x, s + stepSize) + getFn(n, x, s - stepSize);
    deriv /= 3. * gsl_sf_pow_int(stepSize, 2.);
    return deriv;
}

double MultipoleT::getFn(unsigned int n, double x, double s) {
    /** Returns the function f_n(x, s) for curved geometry
      * -> based on recursion
      */
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(0, s);
    }
    // If radius is constant use corresponding recursion
    if (not variableRadius_m) {
        double rho = length_m / angle_m;
        double h_s = 1 + x / rho;
    double temp = 0.0;
    temp += getFnDerivX(n - 1, x, s) / (h_s * rho);
    temp += getFnSecDerivX(n - 1, x, s);
    temp += getFnSecDerivS(n - 1, x, s) / gsl_sf_pow_int(h_s, 2.0);
    temp *= -1.0;
    return temp;
    }
    // If radius is variable use corresponding recursion 
    else {
        double temp = 0.0;
    double h_s = getScaleFactor(x, s);
        temp += getScaleFactorDerivX(x, s) * getFnDerivX(n - 1, x, s);
    temp -= getScaleFactorDerivS(x, s) * getFnDerivS(n - 1, x, s) / 
            gsl_sf_pow_int(h_s, 2.0);
    temp += getFnSecDerivX(n - 1, x, s) * h_s;
    temp += getFnSecDerivS(n - 1, x, s) / h_s;
    temp *= -1.0 / h_s; 
    return temp;
    }
}    
        

inline
void MultipoleT::initialise(PartBunchBase<double, 3>* bunch, 
                double &startField, 
                double &endField) {
    RefPartBunch_m = bunch;
}

inline
bool MultipoleT::bends() const {
    return (transProfile_m[0] != 0);
}

inline
PlanarArcGeometry& MultipoleT::getGeometry() {
    return planarArcGeometry_m;
}

inline
const PlanarArcGeometry& MultipoleT::getGeometry() const {
    return planarArcGeometry_m;
}

inline 
EMField &MultipoleT::getField() {
    return dummy;
}

inline
const EMField &MultipoleT::getField() const {
    return dummy;
}
