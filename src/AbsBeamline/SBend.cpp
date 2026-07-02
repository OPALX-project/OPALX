#include "AbsBeamline/SBend.h"

#include "AbsBeamline/BeamlineVisitor.h"

#include "BeamlineGeometry/Geometry.h"

#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>

SBend::SBend() : SBend("") {}

SBend::SBend(const SBend& right)
    : ElementBase(right),
      normalComponents_m(right.normalComponents_m),
      skewComponents_m(right.skewComponents_m),
      maxNormal_m(right.maxNormal_m),
      maxSkew_m(right.maxSkew_m),
      gap_m(right.gap_m),
      fringeHalfGap_m(right.fringeHalfGap_m),
      fringeIntegral_m(right.fringeIntegral_m),
      designEnergy_m(right.designEnergy_m),
      designEnergyChangeable_m(true) {}

SBend::SBend(const std::string& name)
    : ElementBase(name),
      gap_m(0.0),
      fringeHalfGap_m(0.0),
      fringeIntegral_m(0.5),
      designEnergy_m(0.0),
      designEnergyChangeable_m(true) {}

SBend::~SBend() = default;

void SBend::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    online_m       = true;
}

void SBend::finalise() { online_m = false; }

bool SBend::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
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
            "SBend::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                // Convert to bend coordinates (radial x, y, arc-length s) so the fringe and
                // gate act on the arc length, exact at any bend angle.
                const Vector_t<double, 3> aligned =
                        GeometryHelper::rotateAboutY(Rview(i), -inputs.faceAngle);
                const Vector_t<double, 3> arc =
                        GeometryHelper::toBendArcCoords(aligned, inputs.curvature, inputs.bodyLength);
                if (arc(2) < zBegin || arc(2) > zEnd) {
                    return;
                }
                Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
                        BendFieldModel::bendField(arc, inputs), arc(2), inputs.curvature,
                        inputs.bodyLength);
                Bf = GeometryHelper::rotateAboutY(Bf, inputs.faceAngle);
                for (unsigned d = 0; d < 3; ++d) {
                    Bview(i)(d) += Bf(d);
                }
            });

    return false;
}

bool SBend::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const Vector_t<double, 3> arc = bendCoords(R);
    if (!isInsideArc(arc)) {
        return false;
    }
    if (!isInsideTransverse(arc)) {
        return getFlagDeleteOnTransverseExit();
    }

    computeFieldHost(R, B);
    (void)E;
    return false;
}

bool SBend::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const Vector_t<double, 3> arc = bendCoords(R);
    if (!isInsideArc(arc)) {
        return false;
    }
    if (!isInsideTransverse(arc)) {
        return true;
    }

    computeFieldHost(R, B);
    (void)E;
    return false;
}

void SBend::getFieldExtent(double& zBegin, double& zEnd) const {
    // Single source of the longitudinal field extent: the body plus one Enge fringe
    // half width past each pole face, projected by the face angle. With no gap the
    // half width is zero, so this is the plain body extent [0, L].
    const double half = BendFieldModel::fringeHalfWidth(
            BendFieldModel::profileGap(gap_m, fringeHalfGap_m));
    zBegin = -half / GeometryHelper::safeAbsCos(getGeometry().getEntranceAngle());
    zEnd   = getGeometry().getElementLength() + half / GeometryHelper::safeAbsCos(getGeometry().getExitAngle());
}

// isInside() is inherited from ElementBase (field extent + transverse aperture).

void SBend::computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();
    const Vector_t<double, 3> aligned = GeometryHelper::rotateAboutY(R, -inputs.faceAngle);
    const Vector_t<double, 3> arc =
            GeometryHelper::toBendArcCoords(aligned, inputs.curvature, inputs.bodyLength);
    Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
            BendFieldModel::bendField(arc, inputs), arc(2), inputs.curvature, inputs.bodyLength);
    Bf = GeometryHelper::rotateAboutY(Bf, inputs.faceAngle);
    for (unsigned d = 0; d < 3; ++d) {
        B(d) += Bf(d);
    }
}

Vector_t<double, 3> SBend::bendCoords(const Vector_t<double, 3>& r) const {
    // Remove the entrance pole-face tilt so the arc coordinate is measured relative to the
    // design orbit (the placed frame is face-aligned, tilted by E1).
    const Vector_t<double, 3> aligned = GeometryHelper::rotateAboutY(r, -getGeometry().getEntranceAngle());
    return GeometryHelper::toBendArcCoords(aligned, getGeometry().getCurvature(), getGeometry().getElementLength());
}

bool SBend::isInsideArc(const Vector_t<double, 3>& arc) const {
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);
    return arc(2) >= zBegin && arc(2) < zEnd;
}

bool SBend::isInside(const Vector_t<double, 3>& r) const {
    // Selection/containment uses the arc-length s and the radial offset (not the
    // straight-frame z/x), so the bend stays selected as the orbit curves through it
    // and the aperture is measured relative to the design orbit, not the entry frame.
    const Vector_t<double, 3> arc = bendCoords(r);
    return isInsideArc(arc) && isInsideTransverse(arc);
}

BendFieldModel::FieldInputs SBend::makeFieldInputs() const {
    BendFieldModel::FieldInputs in{};

    // Multipole coefficients (dipole + quadrupole), read once from the device views.
    // The parse side scales them by the reference momentum (p/c) but not the charge; the
    // physical field is B = (p/q)·k, so divide by the reference charge here (the species is
    // only known once a bunch is attached). This gives the correct bend direction for either
    // charge sign. Defaults to q = 1 when no bunch is attached (unit tests set the field
    // directly and are unaffected).
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

    in.bodyLength  = getGeometry().getElementLength();
    in.curvature   = getGeometry().getCurvature();
    in.faceAngle   = getGeometry().getEntranceAngle();
    in.profileGap  = BendFieldModel::profileGap(gap_m, fringeHalfGap_m);
    in.cosEntrance = GeometryHelper::safeAbsCos(getGeometry().getEntranceAngle());
    in.cosExit     = GeometryHelper::safeAbsCos(getGeometry().getExitAngle());

    // Distributed pole-face vertical edge focusing, active only with a fringe. The kick
    // coefficient is spread over the Enge ramp so its integral matches the hard-edge kick.
    in.entryEdgeCoefficient = 0.0;
    in.exitEdgeCoefficient  = 0.0;
    if (in.profileGap > 0.0) {
        const double arcLength = getGeometry().getArcLength();
        const double h         = (arcLength > 0.0) ? getGeometry().getBendAngle() / arcLength : 0.0;
        const double half      = BendFieldModel::fringeHalfWidth(in.profileGap);
        const double span      = std::abs(
                BendFieldModel::engeProfile(-half, in.profileGap).value
                - BendFieldModel::engeProfile(half, in.profileGap).value);
        if (std::abs(h) > 1.0e-15 && span > 1.0e-15) {
            const double rigidity = in.dipoleNormal / h;
            in.entryEdgeCoefficient =
                    rigidity
                    * BendFieldModel::edgeVerticalKickCoefficient(
                            h, fringeHalfGap_m, fringeIntegral_m, getGeometry().getEntranceAngle())
                    / span;
            in.exitEdgeCoefficient =
                    rigidity
                    * BendFieldModel::edgeVerticalKickCoefficient(
                            h, fringeHalfGap_m, fringeIntegral_m, getGeometry().getExitAngle())
                    / span;
        }
    }

    return in;
}

void SBend::setFieldComponents(const std::vector<double>& normal, const std::vector<double>& skew) {
    maxNormal_m = static_cast<int>(normal.size());
    maxSkew_m   = static_cast<int>(skew.size());

    normalComponents_m = Kokkos::View<double*>("SBend::normal", maxNormal_m);
    skewComponents_m   = Kokkos::View<double*>("SBend::skew", maxSkew_m);

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

double SBend::getB() const {
    if (maxNormal_m < 1) {
        return 0.0;
    }
    double val;
    Kokkos::deep_copy(val, Kokkos::subview(normalComponents_m, 0));
    return val;
}

void SBend::setB(double B) {
    if (maxNormal_m < 1) {
        maxNormal_m = 1;
        Kokkos::resize(normalComponents_m, 1);
    }
    Kokkos::deep_copy(Kokkos::subview(normalComponents_m, 0), B);
}

void SBend::accept(BeamlineVisitor& visitor) const { visitor.visitSBend(*this); }

ElementType SBend::getType() const { return ElementType::SBEND; }

