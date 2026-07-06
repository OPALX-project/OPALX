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
      normalComponentsHost_m(right.normalComponentsHost_m),
      skewComponentsHost_m(right.skewComponentsHost_m),
      maxNormal_m(right.maxNormal_m),
      maxSkew_m(right.maxSkew_m),
      gap_m(right.gap_m),
      fringeIntegral_m(right.fringeIntegral_m),
      designEnergy_m(right.designEnergy_m),
      designEnergyChangeable_m(true) {}

SBend::SBend(const std::string& name)
    : ElementBase(name),
      gap_m(0.0),
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
                // Convert (x,y,z) -> (x,y,arc s) in the element's frame.
                const Vector_t<double, 3> arc = GeometryHelper::toBendArcCoords(
                        Rview(i), inputs.curvature, inputs.bodyLength);

                if (arc(2) < zBegin || arc(2) > zEnd) {
                    return; // return this particle's lambda 
                }

                // Rotate the element's field to the entrance frame.
                const Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
                        BendFieldModel::bendField(arc, inputs), arc(2), inputs.curvature,
                        inputs.bodyLength);

                // Apply field.
                Bview(i)(0) += Bf(0);
                Bview(i)(1) += Bf(1);
                Bview(i)(2) += Bf(2);
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
    // Single source of the longitudinal field extent: the body plus one Enge fringe half width
    // past each face. With no gap the half width is zero, so this is the plain body extent [0, L].
    const double half = BendFieldModel::fringeHalfWidth(gap_m);
    zBegin            = -half;
    zEnd              = getGeometry().getElementLength() + half;
}

// isInside() is inherited from ElementBase (field extent + transverse aperture).

void SBend::computeFieldHost(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();
    const Vector_t<double, 3> arc =
            GeometryHelper::toBendArcCoords(R, inputs.curvature, inputs.bodyLength);
    const Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
            BendFieldModel::bendField(arc, inputs), arc(2), inputs.curvature, inputs.bodyLength);
    for (unsigned d = 0; d < 3; ++d) {
        B(d) += Bf(d);
    }
}

Vector_t<double, 3> SBend::bendCoords(const Vector_t<double, 3>& r) const {
    // The stored frame is the design-orbit entrance tangent, so the arc coordinate is measured
    // directly (no pole-face de-tilt).
    return GeometryHelper::toBendArcCoords(
            r, getGeometry().getCurvature(), getGeometry().getElementLength());
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
    // Read the pre-built host mirrors directly; no per-apply device->host copy.
    in.dipoleNormal = ((maxNormal_m > 0) ? normalComponentsHost_m(0) : 0.0) / charge;
    in.quadNormal   = ((maxNormal_m > 1) ? normalComponentsHost_m(1) : 0.0) / charge;
    in.dipoleSkew   = ((maxSkew_m > 0) ? skewComponentsHost_m(0) : 0.0) / charge;
    in.quadSkew     = ((maxSkew_m > 1) ? skewComponentsHost_m(1) : 0.0) / charge;

    in.bodyLength = getGeometry().getElementLength();
    in.curvature  = getGeometry().getCurvature();
    in.profileGap = gap_m;

    // Vertical edge focusing, active only with a fringe. A sector bend's faces are perpendicular
    // to the design orbit (edge angle 0), so only the fringe-field (FINT) term remains; the kick
    // is spread over the Enge ramp so its integral matches the hard-edge kick.
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
            const double rigidity    = in.dipoleNormal / h;
            const double coefficient = rigidity
                                       * BendFieldModel::edgeVerticalKickCoefficient(
                                               h, 0.5 * gap_m, fringeIntegral_m, 0.0)
                                       / span;
            in.entryEdgeCoefficient = coefficient;
            in.exitEdgeCoefficient  = coefficient;
        }
    }

    return in;
}

void SBend::setFieldComponents(const std::vector<double>& normal, const std::vector<double>& skew) {
    maxNormal_m = static_cast<int>(normal.size());
    maxSkew_m   = static_cast<int>(skew.size());

    normalComponents_m = Kokkos::View<double*>("SBend::normal", maxNormal_m);
    skewComponents_m   = Kokkos::View<double*>("SBend::skew", maxSkew_m);

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

double SBend::getB() const {
    if (maxNormal_m < 1) {
        return 0.0;
    }
    return normalComponentsHost_m(0);
}

void SBend::setB(double B) {
    if (maxNormal_m < 1) {
        maxNormal_m            = 1;
        normalComponents_m     = Kokkos::View<double*>("SBend::normal", 1);
        normalComponentsHost_m = Kokkos::create_mirror_view(normalComponents_m);
    }
    normalComponentsHost_m(0) = B;
    Kokkos::deep_copy(Kokkos::subview(normalComponents_m, 0), B);
}

void SBend::accept(BeamlineVisitor& visitor) const { visitor.visitSBend(*this); }

ElementType SBend::getType() const { return ElementType::SBEND; }

