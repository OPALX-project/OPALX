/*
 *  Copyright (c) 2017, Chris Rogers
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

#include <cmath>

#include "AbsBeamline/SpiralSector.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"

SpiralSector::SpiralSector(const std::string &name)
        : Component(name),
         planarArcGeometry_m(1., 1.), dummy() {
    setElType(isDrift);
}

SpiralSector::SpiralSector(const SpiralSector &right)
        : Component(right),
          planarArcGeometry_m(right.planarArcGeometry_m), dummy() {
    RefPartBunch_m = right.RefPartBunch_m;
    setElType(isDrift);
}

SpiralSector::~SpiralSector() {
}

ElementBase* SpiralSector::clone() const {
    return new SpiralSector(*this);
}

EMField &SpiralSector::getField() {
    return dummy;
}

const EMField &SpiralSector::getField() const {
    return dummy;
}

bool SpiralSector::apply(const size_t &i, const double &t,
                    Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

void SpiralSector::initialise(PartBunch *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    calculateDfCoefficients();
}

void SpiralSector::finalise() {
    RefPartBunch_m = NULL;
}

bool SpiralSector::bends() const {
    return true;
}

BGeometryBase& SpiralSector::getGeometry() {
    return planarArcGeometry_m;
}

const BGeometryBase& SpiralSector::getGeometry() const {
    return planarArcGeometry_m;
}

void SpiralSector::accept(BeamlineVisitor& visitor) const {
    visitor.visitSpiralSector(*this);
}


bool SpiralSector::apply(const Vector_t &R, const Vector_t &P,
                    const double &t, Vector_t &E, Vector_t &B) {
    Vector_t pos = R - centre_m;
    double r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    double phi = atan2(pos[2], pos[0]);
    Vector_t posCyl(r, pos[1], phi);
    Vector_t bCyl(0., 0., 0.);
    applyCylindrical(posCyl, P, t, E, bCyl);
    // this is cartesian coordinates
    B[1] = bCyl[1];
    B[0] = bCyl[0]*cos(phi) + bCyl[2]*sin(phi);
    B[2] = bCyl[2]*cos(phi) + bCyl[0]*sin(phi);
    return true;

}

void SpiralSector::calculateDfCoefficients() {
    dfCoefficients_m = std::vector<std::vector<double> >(maxOrder_m);
    dfCoefficients_m[0] = std::vector<double>(1, 1.); // f_0 = 1.*0th derivative
    for (size_t n = 0; n+1 < maxOrder_m; n += 2) { // n indexes the power in z
        dfCoefficients_m[n+1] = std::vector<double>(dfCoefficients_m[n].size()+1, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+1][i+1] = dfCoefficients_m[n][i]/(2*n+1);
        }
        if (n+2 == maxOrder_m) {
            break;
        }
        dfCoefficients_m[n+2] = std::vector<double>(dfCoefficients_m[n].size()+2, 0);
        for(size_t i = 0; i < dfCoefficients_m[n].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] = -(k_m-2*n)*(k_m-2*n)/(n+1)*dfCoefficients_m[n][i];
        }
        for(size_t i = 0; i < dfCoefficients_m[n+1].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] += 2*(k_m-2*n)*tanDelta_m*dfCoefficients_m[n+1][i];
            dfCoefficients_m[n+2][i+1] -= (1+tanDelta_m*tanDelta_m)*dfCoefficients_m[n+1][i];
        }
    }

}

BUG - IT COMPILES AND RUNS; I BELIEVE THE MATHS IS RIGHT; BUT MAXWELL DOES NOT IMPROVE!

STRATEGY... Try Testing M2a-c; validate each of the field elements in turn; then consider M1 again

bool SpiralSector::applyCylindrical(const Vector_t &pos, const Vector_t &P,
                    const double &t, Vector_t &E, Vector_t &B) {
    double r = pos[0];
    double phi = pos[2];
    double normRadius = r/r0_m;
    double g = tanDelta_m*log(normRadius);
    double phiSpiral = phi - g;
    double h = pow(normRadius, k_m)*Bz_m;
    std::vector<double> fringeDerivatives(maxOrder_m, 0.);
    std::cerr << "(r, y, phi) " << pos << " h " << h << " psi " << phiSpiral << std::endl;
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = endField_m->function(phiSpiral, i);
        std::cerr << fringeDerivatives[i] << " ";
    }
    std::cerr << std::endl;
    for (size_t n = 0; n < dfCoefficients_m.size(); n += 2) {
        double f2n = 0;
        Vector_t deltaB;
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            f2n += dfCoefficients_m[n][i]*fringeDerivatives[i];
        }
        double f2nplus1 = 0;
        for (size_t i = 0; i < dfCoefficients_m[n+1].size() && n+1 < dfCoefficients_m.size(); ++i) {
            f2nplus1 += dfCoefficients_m[n+1][i]*fringeDerivatives[i];
        }
        deltaB[1] = f2n*h*pow(pos[1]/r, n); // Bz = sum(f_2n * h * (z/r)^2n
        deltaB[2] = f2nplus1*h*pow(pos[1]/r, n+1); // Bphi = sum(f_2n+1 * h * (z/r)^2n+1
        deltaB[0] = (f2n*(k_m-n)/(n+1) - tanDelta_m*f2nplus1)*h*pow(pos[1]/r, n+1);
        std::cerr << "    " << n << " " << deltaB << std::endl;
        B += deltaB;
    }
    std::cerr << B << std::endl;
    return true;
}

void SpiralSector::setEndField(endfieldmodel::EndFieldModel* endField) {
    if (endField_m != NULL) {
        delete endField_m;
    }

    endField_m = endField;
}

