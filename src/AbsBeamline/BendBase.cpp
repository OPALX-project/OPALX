#include "AbsBeamline/BendBase.h"

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Rotation3D.h"
#include "BeamlineGeometry/Vector3D.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>

namespace {
    CoordinateSystemTrafo toCoordinateSystemTrafo(const Euclid3D& frame) {
        matrix3x3_t rotation;
        const Rotation3D& euclidRotation = frame.getRotation();
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                rotation(row, col) = euclidRotation(row, col);
            }
        }

        const Vector3D& displacement = frame.getVector();
        const Vector_t<double, 3> origin(
                displacement.getX(), displacement.getY(), displacement.getZ());

        return CoordinateSystemTrafo(origin, Quaternion(rotation).conjugate());
    }
}  // namespace

BendBase::BendBase() : BendBase("") {}

BendBase::BendBase(const BendBase& right)
    : Component(right),
      startField_m(right.startField_m),
      endField_m(right.endField_m),
      angle_m(right.angle_m),
      entranceAngle_m(right.entranceAngle_m),
      exitAngle_m(right.exitAngle_m),
      gap_m(right.gap_m),
      designEnergy_m(right.designEnergy_m),
      designEnergyChangeable_m(true),
      fieldAmplitudeX_m(right.fieldAmplitudeX_m),
      fieldAmplitudeY_m(right.fieldAmplitudeY_m),
      fieldAmplitude_m(right.fieldAmplitude_m),
      fileName_m(right.fileName_m),
      entryFaceRotation_m(right.entryFaceRotation_m),
      exitFaceRotation_m(right.exitFaceRotation_m),
      entryFaceCurvature_m(right.entryFaceCurvature_m),
      exitFaceCurvature_m(right.exitFaceCurvature_m),
      slices_m(right.slices_m),
      stepSize_m(right.stepSize_m),
      nSlices_m(right.nSlices_m),
      k1_m(right.k1_m) {}

BendBase::BendBase(const std::string& name)
    : Component(name),
      startField_m(0.0),
      endField_m(0.0),
      angle_m(0.0),
      entranceAngle_m(0.0),
      exitAngle_m(0.0),
      gap_m(0.0),
      designEnergy_m(0.0),
      designEnergyChangeable_m(true),
      fieldAmplitudeX_m(0.0),
      fieldAmplitudeY_m(0.0),
      fieldAmplitude_m(0.0),
      fileName_m(),
      entryFaceRotation_m(0.0),
      exitFaceRotation_m(0.0),
      entryFaceCurvature_m(0.0),
      exitFaceCurvature_m(0.0),
      slices_m(1.0),
      stepSize_m(0.0),
      nSlices_m(1),
      k1_m(0.0) {}

BendBase::~BendBase() = default;

void BendBase::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    startField_m   = startField;
    endField_m     = startField + getElementLength();
    endField       = endField_m;
    online_m       = true;
}

void BendBase::finalise() { online_m = false; }

bool BendBase::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    auto Rview          = pc->R.getView();
    auto Eview          = pc->E.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    const BMultipoleField& field = getField();
    const int order              = field.order();
    Kokkos::View<double*> normal("BendBase::normal", order);
    Kokkos::View<double*> skew("BendBase::skew", order);
    auto normalHost = Kokkos::create_mirror_view(normal);
    auto skewHost   = Kokkos::create_mirror_view(skew);
    for (int i = 0; i < order; ++i) {
        normalHost(i) = field.getNormalComponent(i);
        skewHost(i)   = field.getSkewComponent(i);
    }
    Kokkos::deep_copy(normal, normalHost);
    Kokkos::deep_copy(skew, skewHost);

    const double elemLength = getElementLength();

    Kokkos::parallel_for(
            "BendBase::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                if (Rview(i)(2) < 0.0 || Rview(i)(2) > elemLength) {
                    return;
                }

                Vector_t<double, 3> Bf(0.0);
                const double x = Rview(i)(0);
                const double y = Rview(i)(1);

                if (normal.extent(0) > 0) {
                    Bf(1) += normal(0);
                }
                if (skew.extent(0) > 0) {
                    Bf(0) -= skew(0);
                }
                if (normal.extent(0) > 1) {
                    Bf(0) += normal(1) * y;
                    Bf(1) += normal(1) * x;
                }
                if (skew.extent(0) > 1) {
                    Bf(0) -= skew(1) * x;
                    Bf(1) += skew(1) * y;
                }

                for (unsigned d = 0; d < 3; ++d) {
                    Eview(i)(d) += 0.0;
                    Bview(i)(d) += Bf(d);
                }
            });

    return false;
}

bool BendBase::apply(
        const size_t& i, const double&, Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    const Vector_t<double, 3> R             = pc->R.getView()(i);

    if (!isInside(R)) {
        return false;
    }
    if (!isInsideTransverse(R)) {
        return getFlagDeleteOnTransverseExit();
    }

    computeFieldHost(R, getField(), B);
    (void)E;
    return false;
}

bool BendBase::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    if (!isInside(R)) {
        return false;
    }
    if (!isInsideTransverse(R)) {
        return getFlagDeleteOnTransverseExit();
    }

    computeFieldHost(R, getField(), B);
    (void)E;
    return false;
}

bool BendBase::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    if (!isInside(R)) {
        return false;
    }
    if (!isInsideTransverse(R)) {
        return true;
    }

    computeFieldHost(R, getField(), B);
    (void)E;
    return false;
}

void BendBase::getFieldExtend(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd   = getElementLength();
}

bool BendBase::isInside(const Vector_t<double, 3>& r) const {
    return r(2) >= 0.0 && r(2) < getElementLength() && isInsideTransverse(r);
}

CoordinateSystemTrafo BendBase::getEdgeToBegin() const {
    return toCoordinateSystemTrafo(getEntranceFrame());
}

CoordinateSystemTrafo BendBase::getEdgeToEnd() const {
    return toCoordinateSystemTrafo(getExitFrame());
}

double BendBase::getChordLength() const {
    const Euclid3D entrance = getEntranceFrame();
    const Euclid3D exit     = getExitFrame();
    const Vector3D delta    = exit.getVector() - entrance.getVector();
    return std::sqrt(
            delta.getX() * delta.getX() + delta.getY() * delta.getY()
            + delta.getZ() * delta.getZ());
}

std::vector<Vector_t<double, 3>> BendBase::getDesignPath(std::size_t minSamples) const {
    const double sBegin = getEntrance();
    const double sEnd   = getExit();
    const double span   = std::abs(sEnd - sBegin);
    const std::size_t samples =
            std::max<std::size_t>(minSamples, static_cast<std::size_t>(std::ceil(span / 0.01)) + 1);

    std::vector<Vector_t<double, 3>> path;
    path.reserve(samples);
    for (std::size_t i = 0; i < samples; ++i) {
        const double alpha =
                (samples > 1) ? static_cast<double>(i) / static_cast<double>(samples - 1) : 0.0;
        const double s      = sBegin + alpha * (sEnd - sBegin);
        const Euclid3D pose = getTransform(s);
        path.emplace_back(pose.getX(), pose.getY(), pose.getZ());
    }

    return path;
}

double BendBase::calcDesignRadius(double fieldAmplitude) const {
    const auto& reference  = *RefPartBunch_m->getParticleContainer()->getReference();
    const double mass      = reference.getM();
    const double betaGamma = calcBetaGamma();
    const double charge    = reference.getQ();
    return std::abs(betaGamma * mass / (Physics::c * fieldAmplitude * charge));
}

double BendBase::calcFieldAmplitude(double radius) const {
    const auto& reference  = *RefPartBunch_m->getParticleContainer()->getReference();
    const double mass      = reference.getM();
    const double betaGamma = calcBetaGamma();
    const double charge    = reference.getQ();
    return betaGamma * mass / (Physics::c * radius * charge);
}

double BendBase::calcBendAngle(double chordLength, double radius) const {
    return 2.0 * std::asin(chordLength / (2.0 * radius));
}

double BendBase::calcDesignRadius(double chordLength, double angle) const {
    return chordLength / (2.0 * std::sin(angle / 2.0));
}

double BendBase::calcGamma() const {
    const auto& reference = *RefPartBunch_m->getParticleContainer()->getReference();
    const double mass     = reference.getM();
    return designEnergy_m / mass + 1.0;
}

double BendBase::calcBetaGamma() const {
    const double gamma = calcGamma();
    return std::sqrt(gamma * gamma - 1.0);
}

void BendBase::computeFieldHost(
        const Vector_t<double, 3>& R, const BMultipoleField& field, Vector_t<double, 3>& B) {
    const int order = field.order();
    if (order > 0) {
        B(1) += field.getNormalComponent(0);
        B(0) -= field.getSkewComponent(0);
    }
    if (order > 1) {
        B(0) += field.getNormalComponent(1) * R(1);
        B(1) += field.getNormalComponent(1) * R(0);
        B(0) -= field.getSkewComponent(1) * R(0);
        B(1) += field.getSkewComponent(1) * R(1);
    }
}
