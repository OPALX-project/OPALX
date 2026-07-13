#include "AbsBeamline/RBend.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/Geometry.h"
#include "PartBunch/PartBunch.h"

RBend::RBend() : RBend("") {}

RBend::RBend(const RBend& right)
    : ElementBase(right),
      normalComponents_m(right.normalComponents_m),
      skewComponents_m(right.skewComponents_m),
      normalComponentsHost_m(right.normalComponentsHost_m),
      skewComponentsHost_m(right.skewComponentsHost_m),
      maxNormal_m(right.maxNormal_m),
      maxSkew_m(right.maxSkew_m),
      gap_m(right.gap_m),
      fringeIntegral_m(right.fringeIntegral_m),
      designEnergy_m(right.designEnergy_m),
      designEnergyChangeable_m(true) {}

RBend::RBend(const std::string& name)
    : ElementBase(name),
      gap_m(0.0),
      fringeIntegral_m(0.5),
      designEnergy_m(0.0),
      designEnergyChangeable_m(true) {}

RBend::~RBend() = default;

void RBend::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    online_m       = true;
}

void RBend::finalise() { online_m = false; }

bool RBend::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    // Calc. field extent.
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);

    // Field info capture by value for kernel.
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();

    Kokkos::parallel_for(
            "RBend::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const Vector_t<double, 3>& point = Rview(i);
                if (point(2) < zBegin || point(2) > zEnd) {
                    return;  // exits the lambda for current particle
                }

                const Vector_t<double, 3> Bf = BendFieldModel::bendField(point, inputs);

                Bview(i)(0) += Bf(0);
                Bview(i)(1) += Bf(1);
                Bview(i)(2) += Bf(2);
            });

    return false;
}

bool RBend::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    double zBegin = 0.0, zEnd = 0.0;
    getFieldExtent(zBegin, zEnd);
    if (R(2) < zBegin || R(2) >= zEnd) {
        return false;
    }
    if (!isInsideTransverse(R)) {
        return getFlagDeleteOnTransverseExit();
    }

    computeFieldHost(R, B);
    (void)E;
    return false;
}

bool RBend::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    double zBegin = 0.0, zEnd = 0.0;
    getFieldExtent(zBegin, zEnd);
    if (R(2) < zBegin || R(2) >= zEnd) {
        return false;
    }
    if (!isInsideTransverse(R)) {
        return true;
    }

    computeFieldHost(R, B);
    (void)E;
    return false;
}

void RBend::getFieldExtent(double& zBegin, double& zEnd) const {
    // Element length plus one Enge fringe half width past each pole face
    const double half = BendFieldModel::fringeHalfWidth(gap_m);
    zBegin            = -half;
    zEnd              = getGeometry().getElementLength() + half;
}

void RBend::computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();
    const Vector_t<double, 3> Bf             = BendFieldModel::bendField(R, inputs);
    for (unsigned d = 0; d < 3; ++d) {
        B(d) += Bf(d);
    }
}

double RBend::referenceCurvature() const {
    // Design-orbit curvature 1/rho = angle / arc length = 2 sin(angle/2) / L. Used only to
    // scale the pole-face edge-focusing kick; the field itself is uniform in the box.
    const double arc = getGeometry().getArcLength();
    return (arc > 1.0e-12) ? getGeometry().getBendAngle() / arc : 0.0;
}

double RBend::edgeAngleEntrance() const { return 0.5 * getGeometry().getBendAngle(); }

double RBend::edgeAngleExit() const { return 0.5 * getGeometry().getBendAngle(); }

bool RBend::isInside(const Vector_t<double, 3>& r) const {
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);
    return r(2) >= zBegin && r(2) < zEnd && isInsideTransverse(r);
}

BendFieldModel::FieldInputs RBend::makeFieldInputs() const {
    BendFieldModel::FieldInputs in{};

    // Multipole coefficients (dipole + quadrupole), read once from the device views.
    // Divide by the reference charge so the physical field is B = (p/q)·k (see SBend);
    // the species is only known once a bunch is attached. Defaults to q = 1 (unit tests).
    double charge = 1.0;
    if (RefPartBunch_m != nullptr) {
        const double q = RefPartBunch_m->getParticleContainer()->getReference()->getQ();
        if (std::abs(q) > 1.0e-15) {
            charge = q;
        }
    }
    // Read the pre-built host mirrors directly; no per-apply device->host copy.
    in.dipoleNormal = ((maxNormal_m > 0) ? normalComponentsHost_m(0) : 0.0) / charge;
    in.quadNormal   = ((maxNormal_m > 1) ? normalComponentsHost_m(1) : 0.0) / charge;
    in.dipoleSkew   = ((maxSkew_m > 0) ? skewComponentsHost_m(0) : 0.0) / charge;
    in.quadSkew     = ((maxSkew_m > 1) ? skewComponentsHost_m(1) : 0.0) / charge;

    // RBend geometry
    in.bodyLength = getGeometry().getElementLength();
    in.curvature  = 0.0;
    in.profileGap = gap_m;

    // Compute Enge ramp kick coefficient.
    in.entryEdgeCoefficient = 0.0;
    in.exitEdgeCoefficient  = 0.0;
    if (in.profileGap > 0.0) {
        const double h    = referenceCurvature();
        const double half = BendFieldModel::fringeHalfWidth(in.profileGap);
        const double span = std::abs(
                BendFieldModel::engeProfile(-half, in.profileGap).value
                - BendFieldModel::engeProfile(half, in.profileGap).value);
        if (std::abs(h) > 1.0e-15 && span > 1.0e-15) {
            const double rigidity   = in.dipoleNormal / h;
            in.entryEdgeCoefficient = rigidity
                                      * BendFieldModel::edgeVerticalKickCoefficient(
                                              h, 0.5 * gap_m, fringeIntegral_m, edgeAngleEntrance())
                                      / span;
            in.exitEdgeCoefficient = rigidity
                                     * BendFieldModel::edgeVerticalKickCoefficient(
                                             h, 0.5 * gap_m, fringeIntegral_m, edgeAngleExit())
                                     / span;
        }
    }

    return in;
}

void RBend::setFieldComponents(const std::vector<double>& normal, const std::vector<double>& skew) {
    maxNormal_m = static_cast<int>(normal.size());
    maxSkew_m   = static_cast<int>(skew.size());

    normalComponents_m = Kokkos::View<double*>("RBend::normal", maxNormal_m);
    skewComponents_m   = Kokkos::View<double*>("RBend::skew", maxSkew_m);

    normalComponentsHost_m = Kokkos::create_mirror_view(normalComponents_m);
    skewComponentsHost_m   = Kokkos::create_mirror_view(skewComponents_m);
    for (int i = 0; i < maxNormal_m; ++i) {
        normalComponentsHost_m(i) = normal[i];
    }
    for (int i = 0; i < maxSkew_m; ++i) {
        skewComponentsHost_m(i) = skew[i];
    }
    Kokkos::deep_copy(normalComponents_m, normalComponentsHost_m);
    Kokkos::deep_copy(skewComponents_m, skewComponentsHost_m);
}

double RBend::getB() const {
    if (maxNormal_m < 1) {
        return 0.0;
    }
    return normalComponentsHost_m(0);
}

void RBend::setB(double B) {
    if (maxNormal_m < 1) {
        maxNormal_m            = 1;
        normalComponents_m     = Kokkos::View<double*>("RBend::normal", 1);
        normalComponentsHost_m = Kokkos::create_mirror_view(normalComponents_m);
    }
    normalComponentsHost_m(0) = B;
    Kokkos::deep_copy(Kokkos::subview(normalComponents_m, 0), B);
}

void RBend::accept(BeamlineVisitor& visitor) const { visitor.visitRBend(*this); }

ElementType RBend::getType() const { return ElementType::RBEND; }
