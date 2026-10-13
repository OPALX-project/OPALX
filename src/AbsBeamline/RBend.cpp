#include "AbsBeamline/RBend.h"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/Geometry.h"
#include "PartBunch/PartBunch.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/OpalException.h"

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
      designEnergyChangeable_m(right.designEnergyChangeable_m),
      legacy_m(right.legacy_m),
      legacyFilename_m(right.legacyFilename_m),
      legacyPlacementLength_m(right.legacyPlacementLength_m) {}

RBend::RBend(const std::string& name)
    : ElementBase(name),
      gap_m(0.0),
      fringeIntegral_m(0.5),
      designEnergy_m(0.0),
      designEnergyChangeable_m(true) {}

RBend::~RBend() = default;

void RBend::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    if (legacy_m.enabled) {
        const auto* reference = bunch->getParticleContainer()->getReference();
        initialiseLegacyField(reference->getM(), reference->getQ());
        initialiseLegacyReferencePath(reference->getM(), reference->getQ(), bunch->getdT());
    }
    online_m = true;
}

void RBend::finalise() { online_m = false; }

void RBend::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    // Calc. field extent.
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);
    const ApertureType apertureType = aperture_m.first;
    const double apertureX = aperture_m.second[0], apertureY = aperture_m.second[1];

    // Field info capture by value for kernel.
    const BendFieldModel::FieldInputs inputs = makeFieldInputs();

    Kokkos::parallel_for(
            "RBend::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const Vector_t<double, 3>& point = Rview(i);
                // Match isInside(): entrance belongs to this support, exit does not.
                if (point(2) < zBegin || point(2) >= zEnd
                    || !ApertureHelper::isInsideAperture(
                            point(0), point(1), apertureType, apertureX, apertureY)) {
                    return;  // exits the lambda for current particle
                }

                const Vector_t<double, 3> Bf = BendFieldModel::bendField(point, inputs);

                Bview(i)(0) += Bf(0);
                Bview(i)(1) += Bf(1);
                Bview(i)(2) += Bf(2);
            });
}

void RBend::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    double zBegin = 0.0, zEnd = 0.0;
    getFieldExtent(zBegin, zEnd);
    if (R(2) < zBegin || R(2) >= zEnd) {
        return;
    }
    if (!ApertureHelper::isInsideAperture(R, aperture_m)) {
        return;
    }

    computeFieldHost(R, B);
    (void)E;
}

bool RBend::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    double zBegin = 0.0, zEnd = 0.0;
    getFieldExtent(zBegin, zEnd);
    if (R(2) < zBegin || R(2) >= zEnd) {
        return false;
    }
    if (!ApertureHelper::isInsideAperture(R, aperture_m)) {
        return true;
    }

    computeFieldHost(R, B);
    (void)E;
    return false;
}

void RBend::getFieldExtent(double& zBegin, double& zEnd) const {
    if (legacy_m.enabled) {
        zBegin = legacy_m.entrance.begin;
        zEnd   = getGeometry().getElementLength() + legacy_m.exit.end;
        return;
    }
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
    return r(2) >= zBegin && r(2) < zEnd && ApertureHelper::isInsideAperture(r, aperture_m);
}

BendFieldModel::FieldInputs RBend::makeFieldInputs() const {
    BendFieldModel::FieldInputs in{};
    if (legacy_m.enabled) {
        in.legacy     = legacy_m;
        in.bodyLength = getGeometry().getElementLength();
        return in;
    }

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

void RBend::setLegacyFieldMap(const std::string& filename) {
    legacyPlacementLength_m = -1.0;
    legacy_m                = LegacyRBendField::Inputs{};
    legacyFilename_m        = filename;
    legacy_m.enabled        = !filename.empty();
    getGeometry().setLegacyRBend(legacy_m.enabled);
    if (!legacy_m.enabled) return;
    legacy_m.length               = getGeometry().getElementLength();
    legacy_m.exit.coefficients[5] = 0.31811;  // OPAL's default exit coefficient.
    if (filename != "1DPROFILE1-DEFAULT") {
        std::ifstream file(filename);
        if (!file) throw OpalException("RBend::setLegacyFieldMap", "Cannot open " + filename);
        // Accept blank lines and comments without weakening token validation.
        std::ostringstream cleaned;
        std::string line;
        while (std::getline(file, line)) {
            auto comment = line.find('#');
            auto slash   = line.find("//");
            line.erase(
                    std::min(comment, slash) == std::string::npos ? line.size()
                                                                  : std::min(comment, slash));
            cleaned << line << '\n';
        }
        std::istringstream data(cleaned.str());
        std::string type;
        int entryFlag, exitFlag;
        auto& a = legacy_m.entrance;
        auto& b = legacy_m.exit;
        if (!(data >> type >> a.order >> b.order >> legacy_m.gap) || type != "1DProfile1"
            || a.order < 0 || b.order < 0 || a.order > LegacyRBendField::maximumOrder
            || b.order > LegacyRBendField::maximumOrder
            || !(data >> a.begin >> a.origin >> a.end >> entryFlag >> b.begin >> b.origin >> b.end
                 >> exitFlag))
            throw OpalException(
                    "RBend::setLegacyFieldMap", "Invalid 1DProfile1 header in " + filename);
        for (auto* p : {&a, &b}) {
            for (int i = 0; i <= p->order; ++i)
                if (!(data >> p->coefficients[i]) || !std::isfinite(p->coefficients[i]))
                    throw OpalException(
                            "RBend::setLegacyFieldMap", "Invalid Enge coefficient in " + filename);
            p->begin *= 0.01;
            p->origin *= 0.01;
            p->end *= 0.01;
            if (!std::isfinite(p->begin) || !std::isfinite(p->origin) || !std::isfinite(p->end)
                || !(p->begin < p->origin && p->origin < p->end))
                throw OpalException(
                        "RBend::setLegacyFieldMap", "Invalid fringe bounds in " + filename);
        }
        std::string extra;
        if (data >> extra)
            throw OpalException("RBend::setLegacyFieldMap", "Trailing data in " + filename);
        legacy_m.gap *= 0.01;
        if (!std::isfinite(legacy_m.gap) || legacy_m.gap <= 0.0)
            throw OpalException(
                    "RBend::setLegacyFieldMap", "Profile gap must be positive in " + filename);
        // Match OPAL's orientation convention: each polynomial falls outward.
        for (auto* p : {&a, &b})
            if (LegacyRBendField::profile(p->begin, legacy_m.gap, *p).f
                < LegacyRBendField::profile(p->end, legacy_m.gap, *p).f)
                for (int i = 1; i <= p->order; i += 2)
                    p->coefficients[i] *= -1.0;
    }
    if (gap_m > 0.0) {
        const double ratio = gap_m / legacy_m.gap;
        for (auto* p : {&legacy_m.entrance, &legacy_m.exit}) {
            p->begin = p->origin + (p->begin - p->origin) * ratio;
            p->end   = p->origin + (p->end - p->origin) * ratio;
        }
        legacy_m.gap = gap_m;
    } else {
        gap_m = legacy_m.gap;
    }
}

namespace {
    double integratedProfile(const LegacyRBendField::Inputs& in) {
        // Split at every change of branch; midpoint quadrature avoids counting
        // a hard support edge twice. Richardson extrapolation gives O(h^4).
        std::vector<double> cuts = {
                in.entrance.begin, in.entrance.end, in.length + in.exit.begin,
                in.length + in.exit.end, 0.5 * (in.length + in.entrance.origin + in.exit.origin)};
        std::sort(cuts.begin(), cuts.end());
        const double lo = in.entrance.begin, hi = in.length + in.exit.end;
        auto sum = [&](int n) {
            double integral = 0.0;
            for (size_t j = 1; j < cuts.size(); ++j) {
                const double a = std::max(lo, cuts[j - 1]), b = std::min(hi, cuts[j]);
                if (b <= a) continue;
                const double step = (b - a) / n;
                for (int i = 0; i < n; ++i)
                    integral += step * LegacyRBendField::envelope(a + (i + 0.5) * step, in).f;
            }
            return integral;
        };
        return (4.0 * sum(2048) - sum(1024)) / 3.0;
    }
}  // namespace

void RBend::initialiseLegacyField(double mass, double charge) {
    if (!legacy_m.enabled) return;
    if (!std::isfinite(mass) || mass <= 0.0 || !std::isfinite(charge) || charge == 0.0
        || !std::isfinite(designEnergy_m) || designEnergy_m <= 0.0)
        throw OpalException(
                "RBend::initialiseLegacyField",
                "Positive DESIGNENERGY and a charged design particle are required.");
    // Reload so repeated initialisation never rescales bounds or shifts twice.
    setLegacyFieldMap(legacyFilename_m);
    const double length = getGeometry().getElementLength();
    if (legacyFilename_m == "1DPROFILE1-DEFAULT") {
        const auto original = legacy_m;
        double lo = -2.0 * gap_m, hi = 2.0 * gap_m;
        auto shift = [&](double delta) {
            legacy_m = original;
            legacy_m.entrance.begin += delta;
            legacy_m.entrance.origin += delta;
            legacy_m.entrance.end += delta;
            legacy_m.exit.begin -= delta;
            legacy_m.exit.origin -= delta;
            legacy_m.exit.end -= delta;
        };
        shift(lo);
        const double lowerIntegral = integratedProfile(legacy_m);
        shift(hi);
        if (lowerIntegral < length || integratedProfile(legacy_m) > length)
            throw OpalException(
                    "RBend::initialiseLegacyField",
                    "Cannot match default fringe effective length; supply an explicit 1DProfile1 "
                    "map.");
        for (int i = 0; i < 50; ++i) {
            const double mid = 0.5 * (lo + hi);
            shift(mid);
            if (integratedProfile(legacy_m) > length)
                lo = mid;
            else
                hi = mid;
        }
        shift(0.5 * (lo + hi));
    }
    const double integral = integratedProfile(legacy_m);
    if (!std::isfinite(integral) || integral <= 0.0)
        throw OpalException(
                "RBend::initialiseLegacyField", "Fringe profile has no positive integrated field.");
    const double angle = getGeometry().getBendAngle();
    const double e1    = getGeometry().getEntranceAngle();
    const double rigidity =
            std::sqrt(designEnergy_m * (designEnergy_m + 2.0 * mass)) / (Physics::c * charge);
    // For parallel faces dp_x/dz = -q B_y. This integral matches the design
    // entrance/exit directions without a step-size-dependent tracking search.
    legacy_m.amplitude = rigidity * (std::sin(e1) + std::sin(angle - e1)) / integral;
}

void RBend::initialiseLegacyReferencePath(double mass, double charge, double timeStep) {
    if (!legacy_m.enabled) return;
    legacyPlacementLength_m = -1.0;
    if (!std::isfinite(timeStep) || timeStep <= 0.0 || !std::isfinite(mass) || mass <= 0.0
        || !std::isfinite(charge) || charge == 0.0 || !std::isfinite(designEnergy_m)
        || designEnergy_m <= 0.0)
        throw OpalException(
                "RBend::initialiseLegacyReferencePath", "Invalid design particle or time step.");
    const auto& geometry = getGeometry();
    const double e1      = geometry.getEntranceAngle();
    const double e2      = geometry.getBendAngle() - e1;
    if (std::cos(e1) <= 0.0 || std::cos(e2) <= 0.0)
        throw OpalException(
                "RBend::initialiseLegacyReferencePath",
                "The design orbit must cross both pole faces forwards.");
    const double entryWidth = legacy_m.entrance.origin - legacy_m.entrance.begin;
    const double exitWidth  = legacy_m.exit.end - legacy_m.exit.origin;
    const double span =
            geometry.getArcLength() + entryWidth / std::cos(e1) + exitWidth / std::cos(e2);
    const double bg    = std::sqrt(designEnergy_m * (designEnergy_m + 2.0 * mass)) / mass;
    const double gamma = std::sqrt(1.0 + bg * bg);
    const double step  = bg / gamma * Physics::c * timeStep;
    if (!std::isfinite(span) || !std::isfinite(step) || span + step <= span)
        throw OpalException(
                "RBend::initialiseLegacyReferencePath", "Nonfinite path or unresolved time step.");
    const double halfDrift = 0.5 * Physics::c * timeStep / gamma;
    Vector_t<double, 3> r(-entryWidth * std::tan(e1), 0.0, -entryWidth);
    Vector_t<double, 3> p      = bg * Vector_t<double, 3>(std::sin(e1), 0.0, std::cos(e1));
    const double entryDistance = std::sqrt(dot(r, r));
    const auto exitPoint       = geometry.getEdgeToEnd().getOrigin();
    const BorisPusher pusher;
    for (double s = 0.0; s < span; s += step) {
        r += halfDrift * p;
        const auto b = LegacyRBendField::field(r, legacy_m);
        pusher.kick(r, p, Vector_t<double, 3>(0.0), b, timeStep, mass, charge);
        r += halfDrift * p;
    }
    // OPAL-T subtracts the exterior straight distances from the nominal field
    // span, using the endpoint of the last FULL integration step. In particular,
    // do not replace span by the rounded-up tracked distance: that would change
    // the legacy timing convention by up to one step per bend.
    const Vector_t<double, 3> tail = r - exitPoint;
    const double length            = span - entryDistance - std::sqrt(dot(tail, tail));
    if (!std::isfinite(length) || length <= 0.0)
        throw OpalException(
                "RBend::initialiseLegacyReferencePath", "Invalid effective reference path length.");
    legacyPlacementLength_m = length;
}

double RBend::getPlacementPathLength() const {
    if (!legacy_m.enabled) return getGeometry().getArcLength();
    if (legacyPlacementLength_m < 0.0)
        throw OpalException(
                "RBend::getPlacementPathLength",
                "Initialize the legacy reference path before resolving ELEMEDGE placement.");
    return legacyPlacementLength_m;
}
