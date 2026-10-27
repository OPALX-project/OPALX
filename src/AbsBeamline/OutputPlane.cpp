// Copyright (c) 2023, Chris Rogers
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "AbsBeamline/OutputPlane.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/LossDataSink.h"

extern Inform* gmsg;

OutputPlane::OutputPlane() : OutputPlane("") {}

OutputPlane::OutputPlane(const std::string& name)
    : PluginElement(name),
      maxIterations_m(10),
      algorithm_m(algorithm::INTERPOLATION),
      verbose_m(0) {}

OutputPlane::OutputPlane(const OutputPlane& right) : PluginElement(right), maxIterations_m(10) {
    field_m            = right.field_m;
    normal_m           = right.normal_m;
    centre_m           = right.centre_m;
    radialExtent_m     = right.radialExtent_m;
    verticalExtent_m   = right.verticalExtent_m;
    horizontalExtent_m = right.horizontalExtent_m;
    maxIterations_m    = right.maxIterations_m;
    tolerance_m        = right.tolerance_m;
    geom_m             = right.geom_m;
    recentre_m         = right.recentre_m;
    algorithm_m        = right.algorithm_m;
    verbose_m          = right.verbose_m;
}

OutputPlane::~OutputPlane() {}

ElementBase* OutputPlane::clone() const {
    ElementBase* element = new OutputPlane(*this);
    return element;
}

void OutputPlane::accept(BeamlineVisitor& visitor) const { visitor.visitOutputPlane(*this); }

void OutputPlane::doInitialise(PartBunch_t* /*bunch*/) {
    *gmsg << "* Initialize OutputPlane at " << centre_m << " with normal " << normal_m << "\n*";
    if (radialExtent_m >= 0) {
        *gmsg << " Radial extent " << radialExtent_m;
    }
    if (horizontalExtent_m >= 0) {
        *gmsg << " Horizontal extent " << horizontalExtent_m;
    }
    if (verticalExtent_m >= 0) {
        *gmsg << " Vertical extent " << verticalExtent_m;
    }
    if (field_m != nullptr) {
        *gmsg << " Using field map " << field_m->getName();
    } else {
        *gmsg << " Using empty field map";
    }
    if (recentre_m >= 0) {
        *gmsg << " Recentre event " << recentre_m;
    }
    *gmsg << endl;
}

void OutputPlane::doGoOffline() { *gmsg << "* OutputPlane goes offline " << getName() << endl; }

KOKKOS_FUNCTION void OutputPlane::RK4Step(
        const double& tStep, const double& chargeToMass, const double& t, Vector_t<double>& R1,
        Vector_t<double>& P1) {
    const double thalfStep = tStep / 2.;
    const double tPlusHalf = t + thalfStep;
    const double tPlusStep = t + tStep;
    // f = dy/dt
    // f1 = f(x,t)
    double deriv1[6];

    getDerivatives(R1, P1, t, chargeToMass, deriv1);
    // f2 = f(x+dt*f1/2, t+dt/2 )
    double deriv2[6];
    const Vector_t<double> R2(
            {R1[0] + thalfStep * deriv1[0], R1[1] + thalfStep * deriv1[1],
             R1[2] + thalfStep * deriv1[2]});
    const Vector_t<double> P2(
            {P1[0] + thalfStep * deriv1[3], P1[1] + thalfStep * deriv1[4],
             P1[2] + thalfStep * deriv1[5]});
    getDerivatives(R2, P2, tPlusHalf, chargeToMass, deriv2);

    // f3 = f( x+dt*f2/2, t+dt/2 )
    double deriv3[6];
    const Vector_t<double> R3(
            {R1[0] + thalfStep * deriv2[0], R1[1] + thalfStep * deriv2[1],
             R1[2] + thalfStep * deriv2[2]});
    const Vector_t<double> P3(
            {P1[0] + thalfStep * deriv2[3], P1[1] + thalfStep * deriv2[4],
             P1[2] + thalfStep * deriv2[5]});
    getDerivatives(R3, P3, tPlusHalf, chargeToMass, deriv3);

    // f4 = f(x+dt*f3, t+dt ).
    double deriv4[6];
    const Vector_t<double> R4(
            {R1[0] + tStep * deriv3[0], R1[1] + tStep * deriv3[1], R1[2] + tStep * deriv3[2]});
    const Vector_t<double> P4(
            {P1[0] + tStep * deriv3[3], P1[1] + tStep * deriv3[4], P1[2] + tStep * deriv3[5]});
    getDerivatives(R4, P4, tPlusStep, chargeToMass, deriv4);

    // Return x(t+dt) computed from fourth-order R-K.
    R1[0] += (deriv1[0] + deriv4[0] + 2. * (deriv2[0] + deriv3[0])) * tStep / 6.;
    R1[1] += (deriv1[1] + deriv4[1] + 2. * (deriv2[1] + deriv3[1])) * tStep / 6.;
    R1[2] += (deriv1[2] + deriv4[2] + 2. * (deriv2[2] + deriv3[2])) * tStep / 6.;
    P1[0] += (deriv1[3] + deriv4[3] + 2. * (deriv2[3] + deriv3[3])) * tStep / 6.;
    P1[1] += (deriv1[4] + deriv4[4] + 2. * (deriv2[4] + deriv3[4])) * tStep / 6.;
    P1[2] += (deriv1[5] + deriv4[5] + 2. * (deriv2[5] + deriv3[5])) * tStep / 6.;
}

KOKKOS_FUNCTION void OutputPlane::getDerivatives(
        const Vector_t<double>& R, const Vector_t<double>& P, const double& t,
        const double& chargeToMass, double* yp) {
    const double gamma = std::sqrt(1 + (P[0] * P[0] + P[1] * P[1] + P[2] * P[2]));

    Vector_t<double> beta = P / gamma;
    const double betaX    = beta[0];
    const double betaY    = beta[1];
    const double betaZ    = beta[2];

    yp[0] = Physics::c * betaX;
    yp[1] = Physics::c * betaY;
    yp[2] = Physics::c * betaZ;

    if (field_m == nullptr) {
        throw OpalException("OutputPlane::getDerivatives", "Field was null");
    }
    Vector_t<double> externalB, externalE;
    field_m->apply(R, P, t, externalE, externalB);
    externalB *= Units::kG2T;
    externalE *= Units::kV2V / Units::mm2m / Physics::c;

    yp[3] = chargeToMass * (externalB(2) * betaY - externalB(1) * betaZ + externalE(0));
    yp[4] = chargeToMass * (externalB(0) * betaZ - externalB(2) * betaX + externalE(1));
    yp[5] = chargeToMass * (externalB(1) * betaX - externalB(0) * betaY + externalE(2));
}

KOKKOS_FUNCTION bool OutputPlane::checkOne(
        const int index, const double tStep, const double chargeToMass, double& t,
        Vector_t<double>& R, Vector_t<double>& P, Vector_t<double> centre,
                  Vector_t<double> normal) {
    // distance from particle to the output plane
    // time units are ns
    Vector_t<double> delta = R - centre;
    const double distance  = normal.dot(delta);

    // maximum step, assuming no curvature, rough guess for relativistic beta
    double betaEstimate = euclidean_norm(P);
    if (betaEstimate > 1) {
        betaEstimate = 1.0;
    }
    const double sStep = tStep * betaEstimate * Physics::c;
    if (std::abs(distance) > sStep) {
        // we can't cross the plane
        return false;
    }

    Vector_t rTest(R);
    Vector_t pTest(P);
    RK4Step(tStep, chargeToMass, t, rTest, pTest);
    const double distanceTest = normal.dot(rTest - centre);

    if (distance != 0 && distanceTest / distance > 0) {
        // step does not cross the plane - give up
        // note that particle could cross and cross back the plane in a single
        // time-step; in this case the particle is not registered
        return false;
    }

    if (algorithm_m == algorithm::RK4STEP) {
        rk4Test(tStep, chargeToMass, t, R, P);
    } else if (algorithm_m == algorithm::INTERPOLATION) {
        interpolation(t, R, P);
    }

    delta = R - centre_m;
    if (horizontalExtent_m > 0
        && delta[0] * delta[0] + delta[1] * delta[1] > horizontalExtent_m * horizontalExtent_m) {
        // out of horizontal extent. Again, defined in global coordinates
        return false;
    }
    if (verticalExtent_m > 0 && abs(delta[2]) > verticalExtent_m) {
        // out of vertical extent; note that this is defined in global coordinates
        // not in the coordinates of the plane
        return false;
    }
    if (radialExtent_m > 0
        && delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]
                   > radialExtent_m * radialExtent_m) {
        // out of radial extent
        return false;
    }
    if (recentre_m == index) {
        recentre(R, P);
    }
    return true;
}

void OutputPlane::rk4Test(
        double tStep, const double chargeToMass, double& t, Vector_t<double>& R,
        Vector_t<double>& P) const {
    double preStepDistance  = 0.0;
    double postStepDistance = normal_m.dot(R - centre_m);
    size_t iteration        = 0;
    while (std::abs(postStepDistance) > tolerance_m && iteration < maxIterations_m) {
        preStepDistance = postStepDistance;
        RK4Step(tStep, chargeToMass, t, R, P);  // this updates R and P
        if (verbose_m > 2) {
            *gmsg << "    RK4 iteration " << iteration << " step distance " << preStepDistance
                  << " R " << R << " P " << P << " centre " << centre_m
                  << " d: " << Vector_t<double>{R - centre_m} << " t: " << t << " dt: " << tStep
                  << endl;
            Vector_t<double> externalB, externalE;
            field_m->apply(R, P, t, externalE, externalB);
            *gmsg << "        B " << externalB << " [kG] E " << externalE << " [MV/m] " << endl;
        }
        t += tStep;
        postStepDistance = normal_m.dot(R - centre_m);
        if (postStepDistance / preStepDistance < 0) {  // straddling the plane
            // we stepped too far; step in opposite direction
            // step length in ratio of distance to plane
            tStep *= -abs(postStepDistance) / abs(postStepDistance - preStepDistance);
        } else {
            // we didn't step far enough; step in same direction
            tStep *= abs(postStepDistance) / abs(postStepDistance - preStepDistance);
        }
        iteration++;
    }
}

void OutputPlane::interpolation(double& t, Vector_t<double>& R, Vector_t<double>& P) {
    // trajectory R = R0 + V dt
    // plane = (X-X0).n
    // intersection time t0 = (X0-R0).N / (V.N)
    // intersection position = R0 + V dt
    // relativistic gamma:
    const double gamma              = std::sqrt(1 + P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
    const Vector_t<double> velocity = P / gamma * Physics::c;  // m/ns
    const double dt = Vector_t<double>{centre_m - R}.dot(normal_m) / velocity.dot(normal_m);
    R += velocity * dt;
    t += dt;
}

bool OutputPlane::doCheck(
        PartBunch_t* bunch, const int turnnumber, const double t, const double tstep) {
    for (auto& pc : bunch->getParticleContainers()) {
        const auto rView       = pc->R.getView();
        const auto pView       = pc->P.getView();
        const auto qView       = pc->getQView();
        const auto invalidView = pc->InvalidMask.getView();
        const auto mView       = pc->getMView();
        const bool singleQM    = pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::SingleValue;
        Kokkos::View<bool*> hitView("hits", pc->getLocalNum());
        Kokkos::parallel_for("OutputPlane::doCheck", pc->getLocalNum(),
            KOKKOS_LAMBDA(const size_t i) {
            doOne(i, rView, pView, qView, mView, invalidView, t, tstep, singleQM, hitView);
        });
    }
    return false;
}

KOKKOS_FUNCTION void OutputPlane::doOne(size_t index,
        Kokkos::View<ippl::Vector<double, 3>*> rView,
        Kokkos::View<ippl::Vector<double, 3>*> pView,
        Kokkos::View<double*> qView, Kokkos::View<double*> mView, Kokkos::View<bool*> invalidView,
        double t, double tStep, bool singleQM, Kokkos::View<bool*> hitView) {
    if (!invalidView(index)) {
        double t0{t};
        double q;
        double m;
        if (singleQM) {
            q = qView(0);
            m = mView(0);
        } else {
            q = qView(index);
            m = mView(index);
        }
        double chargeToMass = q / Physics::q_e * Physics::c * Physics::c / m
                              / Units::GeV2eV;  // electron charge cancels
        hitView(index) = checkOne(index, tStep, chargeToMass, t0, rView(index), pView(index));
    }
}


void OutputPlane::recentre(Vector_t R, Vector_t P) {
    setCentre(R);
    setNormal(P);
    recentre_m = -1;  // don't recentre again
}

ElementType OutputPlane::getType() const { return ElementType::OUTPUTPLANE; }
