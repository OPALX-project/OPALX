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

#include "MultipoleTCurvedVarRadius.h"
#include <vector>

#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"
#include "GSLCompat.h"
#include "MultipoleT.h"
#include "tanhDeriv.h"

MultipoleTCurvedVarRadius::MultipoleTCurvedVarRadius(MultipoleT* element)
    : MultipoleTBase(element), varRadiusGeometry_m(1.0, 1.0, 1.0, 1.0, 1.0) {}

void MultipoleTCurvedVarRadius::initialise() {
    // Record geometry information
    varRadiusGeometry_m.setElementLength(
        element_m->getLength() + 2.0 * element_m->getEntryOffset());
    varRadiusGeometry_m.setRadius(element_m->getLength() / element_m->getBendAngle());
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    varRadiusGeometry_m.setS0(s0);
    varRadiusGeometry_m.setLambdaLeft(leftFringe);
    varRadiusGeometry_m.setLambdaRight(rightFringe);
    setMaxOrder(element_m->getMaxFOrder(), element_m->getMaxXOrder());
    // Now work out where the entry point will be in the local cartesian coordinate system
    // whose origin is at the center of the magnet.
    localCartesianEntryPoint_ = curvilinearToLocalCartesian(
        Vector_t<double, 3>{0.0, 0.0, element_m->getLength() / 2.0 + element_m->getEntryOffset()});
    // The tangent to the curve at this point forms the z axis of the coordinate system
    // opal addresses us in, so we can calculate the rotation required
    auto secondPoint = curvilinearToLocalCartesian(
        Vector_t<double, 3>{
            0.0, 0.0, element_m->getLength() / 2.0 + element_m->getEntryOffset() + TangentStep});
    localCartesianRotation_ = -atan2(
        secondPoint[0] - localCartesianEntryPoint_[0],
        secondPoint[2] - localCartesianEntryPoint_[2]);
}

void MultipoleTCurvedVarRadius::getField(
    const Kokkos::View<Vector_t<double, 3>*>& R, Kokkos::View<Vector_t<double, 3>*>& /*E*/,
    Kokkos::View<Vector_t<double, 3>*>& B, const double scaling) {
    // Local variables that are copied into the kernel
    const auto config           = element_m->getConfig();
    const auto tanhCoefficients = tanhCoefficientsDevice_m;
    // This version of the multipole magnet has not been converted to Kokkos.
    // See volume 2, section 4.4 of my notes.  Jon Thompson.
    for (size_t i = 0; i < R.extent(0); ++i) {
        computeBField(R(i), B(i), scaling, config);
    }
}

bool MultipoleTCurvedVarRadius::getField(
    const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B,
    const double scaling) {
    return computeBField(R, B, scaling, element_m->getConfig());
}

Vector_t<double, 3> MultipoleTCurvedVarRadius::toMagnetCoords(
    const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    /** Rotate coordinates around the central axis of the magnet */
    /* TODO:  I'm not sure this is the correct thing to do */
    auto result = rotateFrame(R, config);
    /** Go to local Frenet-Serret coordinates */
    result[2] *= -1;  // OPAL uses a different sign convention...
    transformCoords(result);
    return result;
}

Vector_t<double, 3> MultipoleTCurvedVarRadius::rotateFrame(
    const Vector_t<double, 3>& R, const MultipoleTConfig& config) const {
    Vector_t<double, 3> R_prime(3), R_pprime(3);
    /** Apply two 2D rotation matrices to coordinate vector
     * Rotate around central axis => skew fields
     * Rotate azimuthally => entrance angle
     */
    // 1st rotation
    R_prime[0] = R[0] * std::cos(config.rotation_m) + R[1] * std::sin(config.rotation_m);
    R_prime[1] = -R[0] * std::sin(config.rotation_m) + R[1] * std::cos(config.rotation_m);
    R_prime[2] = R[2];
    // 2nd rotation
    R_pprime[0] = R_prime[2] * std::sin(config.entranceAngle_m)
                  + R_prime[0] * std::cos(config.entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * std::cos(config.entranceAngle_m)
                  - R_prime[0] * std::sin(config.entranceAngle_m);
    return R_pprime;
}

void MultipoleTCurvedVarRadius::transformCoords(Vector_t<double, 3>& R) {
    // Rotate Opal supplied cartesian coordinates around its origin
    auto x_rotated = R[0] * cos(localCartesianRotation_) - R[2] * sin(localCartesianRotation_);
    auto z_rotated = R[0] * sin(localCartesianRotation_) + R[2] * cos(localCartesianRotation_);
    // Offset to the center of the magnet
    R = {x_rotated + localCartesianEntryPoint_[0], R[1], z_rotated + localCartesianEntryPoint_[2]};
    // And finally into curvilinear coordinates
    R = localCartesianToCurvilinear(R);
}

bool MultipoleTCurvedVarRadius::computeBField(
    const Vector_t<double, 3>& R, Vector_t<double, 3>& B, double scaling,
    const MultipoleTConfig& config) {
    const Vector_t<double, 3> RPrime = toMagnetCoords(R, config);
    const bool insideAperture        = Kokkos::abs(RPrime[1]) <= config.verticalAperture_m / 2.0
                                && Kokkos::abs(RPrime[0]) <= config.horizontalAperture_m / 2.0;
    const bool insideBoundingBox = config.boundingBoxLength_m == 0.0
                                   || Kokkos::abs(RPrime[2]) <= config.boundingBoxLength_m / 2.0;
    if (insideAperture && insideBoundingBox) {
        Vector_t<double, 3> result;
        /** Calculate B-field in the local Frenet-Serret frame */
        result[0] = getBx(RPrime, config);
        result[1] = getBz(RPrime, config);
        result[2] = getBs(RPrime, config);
        /** Transform B-field from local to lab coordinates */
        transformBField(result, RPrime);
        result[2] *= -1;  // OPAL uses a different sign convention...
        B += result * scaling;
    }
    return !insideAperture;
}

double MultipoleTCurvedVarRadius::getBs(
    const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    double Bs     = 0.0;
    std::size_t n = config.maxFOrder_m + 1;
    while (n != 0) {
        n--;
        Bs = Bs * R[1] * R[1] + getFnDerivS(n, R[0], R[2], config) / gsl_sf_fact(2 * n + 1);
    }
    Bs *= R[1] / getScaleFactor(R[0], R[2]);
    return Bs;
}

double MultipoleTCurvedVarRadius::getBx(
    const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    double Bx     = 0.0;
    std::size_t n = config.maxFOrder_m + 1;
    while (n != 0) {
        n--;
        Bx = Bx * R[1] * R[1] + getFnDerivX(n, R[0], R[2], config) / gsl_sf_fact(2 * n + 1);
    }
    Bx *= R[1];
    return Bx;
}

double MultipoleTCurvedVarRadius::getBz(
    const Vector_t<double, 3>& R, const MultipoleTConfig& config) {
    double Bz     = 0.0;
    std::size_t n = config.maxFOrder_m + 1;
    while (n != 0) {
        n--;
        Bz = Bz * R[1] * R[1] + getFn(n, R[0], R[2]) / gsl_sf_fact(2 * n);
    }
    return Bz;
}

double MultipoleTCurvedVarRadius::getFnDerivS(
    const std::size_t& n, const double& x, const double& s, const MultipoleTConfig& config) {
    if (n == 0) {
        return getTransDeriv(0, x, config) * getFringeDeriv(1, s);
    }
    double deriv    = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * getFn(n, x, s - 2. * stepSize);
    deriv += -8. * getFn(n, x, s - stepSize);
    deriv += 8. * getFn(n, x, s + stepSize);
    deriv += -1. * getFn(n, x, s + 2. * stepSize);
    deriv /= 12 * stepSize;
    return deriv;
}

double MultipoleTCurvedVarRadius::getFnDerivX(
    const std::size_t& n, const double& x, const double& s, const MultipoleTConfig& config) {
    if (n == 0) {
        return getTransDeriv(1, x, config) * getFringeDeriv(0, s);
    }
    double deriv    = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * getFn(n, x - 2. * stepSize, s);
    deriv += -8. * getFn(n, x - stepSize, s);
    deriv += 8. * getFn(n, x + stepSize, s);
    deriv += -1. * getFn(n, x + 2. * stepSize, s);
    deriv /= 12 * stepSize;
    return deriv;
}

double MultipoleTCurvedVarRadius::getTransDeriv(
    const std::size_t& n, const double& x, const MultipoleTConfig& config) const {
    double func               = 0.0;
    std::size_t transMaxOrder = config.transverseProfileMaxOrder_m;
    if (n > transMaxOrder) {
        return func;
    }
    auto temp = config.transverseProfile_m;
    for (std::size_t i = 1; i <= n; i++) {
        for (std::size_t j = 0; j <= transMaxOrder; j++) {
            if (j <= transMaxOrder - i) {
                temp[j] = temp[j + 1] * static_cast<double>(j + 1);
            } else {
                temp[j] = 0.0;
            }
        }
    }
    std::size_t k = transMaxOrder - n + 1;
    while (k != 0) {
        k--;
        func = func * x + temp[k];
    }
    return func;
}

double MultipoleTCurvedVarRadius::getFringeDeriv(const std::size_t& n, const double& s) {
    double result;
    if (n <= 10) {
        result = (fringeField_l.getTanh(s, static_cast<int>(n))
                  - fringeField_r.getNegTanh(s, static_cast<int>(n)))
                 / 2;
    } else {
        result = tanhderiv::integrate(
            s, fringeField_l.getX0(), fringeField_l.getLambda(), fringeField_r.getLambda(),
            static_cast<int>(n));
    }
    return result;
}

double MultipoleTCurvedVarRadius::getFn(size_t n, double x, double s) {
    double result{};
    if (n == 0) {
        result = element_m->getTransDeriv(0, x) * element_m->getFringeDeriv(0, s);
    } else {
        double rho = element_m->getLength() / element_m->getBendAngle();
        double S_0 = element_m->getFringeDeriv(0, 0);
        double y   = element_m->getFringeDeriv(0, s) / (S_0 * rho);
        std::vector<double> fringeDerivatives;
        for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
            fringeDerivatives.push_back(element_m->getFringeDeriv(j, s) / (S_0 * rho));
        }
        for (std::size_t i = 0; i <= recursion_m.at(n).getMaxXDerivatives(); i++) {
            double temp = 0.0;
            for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
                temp += recursion_m.at(n).evaluatePolynomial(x, y, i, j, fringeDerivatives)
                        * fringeDerivatives.at(j);
            }
            result += temp * element_m->getTransDeriv(i, x);
        }
        result *= gsl_sf_pow_int(-1.0, static_cast<int>(n)) * S_0 * rho;
    }
    return result;
}

Vector_t<double, 3> MultipoleTCurvedVarRadius::localCartesianToOpalCartesian(
    const Vector_t<double, 3>& r) {
    // Offset to the Opal origin
    auto x_offset = r[0] - localCartesianEntryPoint_[0];
    auto z_offset = r[2] - localCartesianEntryPoint_[2];
    // And rotate
    auto x_rotated =
        x_offset * cos(-localCartesianRotation_) - z_offset * sin(-localCartesianRotation_);
    auto z_rotated =
        x_offset * sin(-localCartesianRotation_) + z_offset * cos(-localCartesianRotation_);
    return {x_rotated, r[1], -z_rotated};
}

Vector_t<double, 3> MultipoleTCurvedVarRadius::localCartesianToCurvilinear(
    const Vector_t<double, 3>& r) {
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    double rho                         = element_m->getLength() / element_m->getBendAngle();
    coordinatetransform::CoordinateTransform t(r[0], r[1], r[2], s0, leftFringe, rightFringe, rho);
    std::vector<double> result = t.getTransformation();
    return {result[0], result[1], result[2]};
}

void MultipoleTCurvedVarRadius::transformBField(
    Vector_t<double, 3>& B, const Vector_t<double, 3>& R) {
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    double rho                         = element_m->getLength() / element_m->getBendAngle();
    double prefactor                   = rho * (tanh(s0 / leftFringe) + tanh(s0 / rightFringe));
    double theta                       = leftFringe * log(cosh((R[2] + s0) / leftFringe))
                   - rightFringe * log(cosh((R[2] - s0) / rightFringe));
    theta /= prefactor;
    double Bx = B[0], Bs = B[2];
    B[0] = Bx * cos(theta) - Bs * sin(theta);
    B[2] = Bx * sin(theta) + Bs * cos(theta);
}

void MultipoleTCurvedVarRadius::setMaxOrder(size_t orderZ, size_t orderX) {
    std::size_t N = recursion_m.size();
    while (orderZ >= N) {
        polynomial::RecursionRelationTwo r(N, 2 * (N + orderX + 1));
        r.resizeX(element_m->getTransMaxOrder());
        r.truncate(orderX);
        recursion_m.push_back(r);
        N = recursion_m.size();
    }
}

double MultipoleTCurvedVarRadius::getScaleFactor(double x, double s) {
    double result = 1.0;
    if (element_m->getFringeDeriv(0, s) > 1.0e-12) {
        double radius = element_m->getLength() * element_m->getFringeDeriv(0, 0)
                        / (element_m->getFringeDeriv(0, s) * element_m->getBendAngle());
        result += x / radius;
    }
    return result;
}

double MultipoleTCurvedVarRadius::getFn(unsigned int n, double x, double s) {
    double result{};
    if (n == 0) {
        result = element_m->getTransDeriv(0, x) * element_m->getFringeDeriv(0, s);
    } else {
        double rho = element_m->getLength() / element_m->getBendAngle();
        double S_0 = element_m->getFringeDeriv(0, 0);
        double y   = element_m->getFringeDeriv(0, s) / (S_0 * rho);
        std::vector<double> fringeDerivatives;
        for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
            fringeDerivatives.push_back(element_m->getFringeDeriv(j, s) / (S_0 * rho));
        }
        for (std::size_t i = 0; i <= recursion_m.at(n).getMaxXDerivatives(); i++) {
            double temp = 0.0;
            for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
                temp += recursion_m.at(n).evaluatePolynomial(x, y, i, j, fringeDerivatives)
                        * fringeDerivatives.at(j);
            }
            result += temp * element_m->getTransDeriv(i, x);
        }
        result *= powerInteger(-1.0, n) * S_0 * rho;
    }
    return result;
}

double MultipoleTCurvedVarRadius::reverseTransformResidual(
    const Vector_t<double, 3>& r, const Vector_t<double, 3>& target) {
    // Return the distance between the vector r and the target.
    // We only consider the first and last coordinates as the height coordinate
    // is invariant across these transforms.
    auto c    = localCartesianToCurvilinear(r);
    double dx = c[0] - target[0];
    double ds = c[2] - target[2];
    return sqrt(dx * dx + ds * ds);
}

Vector_t<double, 3> MultipoleTCurvedVarRadius::curvilinearToLocalCartesian(
    const Vector_t<double, 3>& r) {
    // This functions uses a minimize loop and coordinate descent with backtracking
    // to implement the inverse coordinate transform from the magnet's curvilinear
    // system to the local cartesian system whose origins are the centre of the magnet.
    // Note that this function is iterative and should therefore only be used occasionally.
    Vector_t<double, 3> result{r};
    double step     = 1.0;
    double best_res = reverseTransformResidual(result, r);
    for (size_t iter = 0; iter < ReverseTransformMaxIterations; ++iter) {
        bool improved = false;
        for (int dim = 0; dim < 2; ++dim) {
            for (int dir = -1; dir <= 1; dir += 2) {
                Vector_t<double, 3> trial = result;
                if (dim == 0) {
                    trial[0] += dir * step;
                } else {
                    trial[2] += dir * step;
                }
                double res = reverseTransformResidual(trial, r);
                if (res < best_res) {
                    result   = trial;
                    best_res = res;
                    improved = true;
                    break;
                }
            }
            if (improved) {
                break;
            }
        }
        if (!improved) {
            step *= 0.5;
        }
        if (step < ReverseTransformTolerance) {
            break;
        }
    }
    return result;
}
