#include "AbsBeamline/RBend.h"

#include "AbsBeamline/BeamlineVisitor.h"

#include "BeamlineGeometry/Geometry.h"

#include "PartBunch/PartBunch.h"

#include <Kokkos_Core.hpp>

#include <cmath>

RBend::RBend() : RBend("") {}

RBend::RBend(const RBend& right)
    : ElementBase(right),
      normalComponents_m(right.normalComponents_m),
      skewComponents_m(right.skewComponents_m),
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

    // Field-support extent (single source, shared with isInside selection), captured on the
    // host before the kernel launch (getFieldExtent is not device-callable).
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);

    // Coefficients, fringe geometry and edge-focusing built once on the host and captured by value.
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();

    Kokkos::parallel_for(
            "RBend::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                // The local frame is the straight box (+z along the box axis), so the field
                // is a uniform vertical dipole with an Enge ramp in the box z — gate on the
                // box z, no arc conversion. The pusher curves the orbit through it.
                const Vector_t<double, 3>& box = Rview(i);
                if (box(2) < zBegin || box(2) > zEnd) {
                    return;
                }
                const Vector_t<double, 3> Bf = BendFieldModel::bendField(box, inputs);
                for (unsigned d = 0; d < 3; ++d) {
                    Bview(i)(d) += Bf(d);
                }
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
    const double half = BendFieldModel::fringeHalfWidth(
            gap_m);
    zBegin = -half;
    zEnd   = getGeometry().getElementLength() + half;
}

void RBend::computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    // R is already in the straight box frame: a uniform vertical dipole with an Enge ramp
    // in the box z, no arc conversion (mirrors the device kernel in apply(pc)).
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

double RBend::edgeAngleEntrance() const {
    // Angle between the design orbit and the entrance pole face, for the edge focusing:
    // the orbit meets the box face at half the bend angle, plus any explicit rotation E1.
    return 0.5 * getGeometry().getBendAngle() + getGeometry().getEntranceAngle();
}

double RBend::edgeAngleExit() const {
    return 0.5 * getGeometry().getBendAngle() + getGeometry().getExitAngle();
}

bool RBend::isInside(const Vector_t<double, 3>& r) const {
    // Straight box: gate on the box z and the box transverse aperture directly.
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
    auto normalHost = Kokkos::create_mirror_view(normalComponents_m);
    auto skewHost   = Kokkos::create_mirror_view(skewComponents_m);
    Kokkos::deep_copy(normalHost, normalComponents_m);
    Kokkos::deep_copy(skewHost, skewComponents_m);
    in.dipoleNormal = ((maxNormal_m > 0) ? normalHost(0) : 0.0) / charge;
    in.quadNormal   = ((maxNormal_m > 1) ? normalHost(1) : 0.0) / charge;
    in.dipoleSkew   = ((maxSkew_m > 0) ? skewHost(0) : 0.0) / charge;
    in.quadSkew     = ((maxSkew_m > 1) ? skewHost(1) : 0.0) / charge;

    // Straight box frame: the field is a uniform vertical dipole gated on the box z, so no
    // curvature or frame de-tilt is applied (curvature/faceAngle stay 0). The fringe runs
    // over the box length with the faces perpendicular to the box axis (unit projections);
    // pole-face angles enter only the edge-focusing kick below.
    in.bodyLength  = getGeometry().getElementLength();
    in.curvature   = 0.0;
    in.faceAngle   = 0.0;
    in.profileGap  = gap_m;
    in.cosEntrance = 1.0;
    in.cosExit     = 1.0;

    // Distributed pole-face vertical edge focusing, active only with a fringe. The kick
    // coefficient is spread over the Enge ramp so its integral matches the hard-edge kick.
    in.entryEdgeCoefficient = 0.0;
    in.exitEdgeCoefficient  = 0.0;
    if (in.profileGap > 0.0) {
        const double h    = referenceCurvature();
        const double half = BendFieldModel::fringeHalfWidth(in.profileGap);
        const double span = std::abs(
                BendFieldModel::engeProfile(-half, in.profileGap).value
                - BendFieldModel::engeProfile(half, in.profileGap).value);
        if (std::abs(h) > 1.0e-15 && span > 1.0e-15) {
            const double rigidity = in.dipoleNormal / h;
            in.entryEdgeCoefficient =
                    rigidity
                    * BendFieldModel::edgeVerticalKickCoefficient(
                            h, 0.5 * gap_m, fringeIntegral_m, edgeAngleEntrance())
                    / span;
            in.exitEdgeCoefficient =
                    rigidity
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

    auto normalHost = Kokkos::create_mirror_view(normalComponents_m);
    auto skewHost   = Kokkos::create_mirror_view(skewComponents_m);
    for (int i = 0; i < maxNormal_m; ++i) {
        normalHost(i) = normal[i];
    }
    for (int i = 0; i < maxSkew_m; ++i) {
        skewHost(i) = skew[i];
    }
    Kokkos::deep_copy(normalComponents_m, normalHost);
    Kokkos::deep_copy(skewComponents_m, skewHost);
}

double RBend::getB() const {
    if (maxNormal_m < 1) {
        return 0.0;
    }
    double val;
    Kokkos::deep_copy(val, Kokkos::subview(normalComponents_m, 0));
    return val;
}

void RBend::setB(double B) {
    if (maxNormal_m < 1) {
        maxNormal_m = 1;
        Kokkos::resize(normalComponents_m, 1);
    }
    Kokkos::deep_copy(Kokkos::subview(normalComponents_m, 0), B);
}

void RBend::accept(BeamlineVisitor& visitor) const { visitor.visitRBend(*this); }

ElementType RBend::getType() const { return ElementType::RBEND; }
