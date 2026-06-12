#include "AbsBeamline/BendBase.h"

#include <Kokkos_Core.hpp>
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Rotation3D.h"
#include "BeamlineGeometry/Vector3D.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"

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

    /**
     * Recover the signed reference-path phase from entry-frame Cartesian bend coordinates.
     *
     * For the circular reference-path parameterization
     * \f[
     * x + \rho = \rho \cos\phi, \qquad z = \rho \sin\phi,
     * \f]
     * the geometric phase must satisfy \f$\phi = h s\f$ for both positive and
     * negative curvature \f$h = 1/\rho\f$. When \f$\rho < 0\f$, applying
     * `atan2(z, x + \rho)` directly shifts the phase by \f$\pi\f$ because both
     * arguments change sign simultaneously. Multiplying both arguments by the
     * sign of \f$\rho\f$ keeps the branch consistent with the signed bend
     * curvature and therefore recovers the correct local path coordinate.
     */
    KOKKOS_INLINE_FUNCTION double recoverReferencePhaseFromEntryCartesian(
            const Vector_t<double, 3>& entryCartesian, const double curvature) {
        const double radius = 1.0 / curvature;
        const double sign   = (radius >= 0.0) ? 1.0 : -1.0;
        return Kokkos::atan2(sign * entryCartesian(2), sign * (entryCartesian(0) + radius));
    }

    /**
     * @brief Device-compatible entry-frame to bend-field chart conversion.
     *
     * The many-particle slice kernel reconstructs particle positions in the
     * rigid entrance frame. The analytic SBEND field chart is piecewise:
     * straight upstream fringe, circular body, and straight downstream fringe
     * tangent to the body exit. Continuing the circular phase through the exit
     * fringe shortens the sampled downstream support for explicitly placed
     * bends.
     */
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> convertEntryCartesianToFieldLocalDevice(
            const Vector_t<double, 3>& entryCartesian, const double bodyLength,
            const double curvature) {
        if (Kokkos::abs(curvature) <= 1.0e-15) {
            return entryCartesian;
        }

        if (entryCartesian(2) <= 0.0) {
            return entryCartesian;
        }

        const double exitPhi    = curvature * bodyLength;
        const double cosExit    = Kokkos::cos(exitPhi);
        const double sinExit    = Kokkos::sin(exitPhi);
        const double exitX      = (cosExit - 1.0) / curvature;
        const double exitZ      = sinExit / curvature;
        const double dx         = entryCartesian(0) - exitX;
        const double dz         = entryCartesian(2) - exitZ;
        const double exitLocalX = cosExit * dx + sinExit * dz;
        const double exitLocalZ = -sinExit * dx + cosExit * dz;
        if (exitLocalZ >= 0.0) {
            return Vector_t<double, 3>(exitLocalX, entryCartesian(1), bodyLength + exitLocalZ);
        }

        const double radius = 1.0 / curvature;
        const double phi    = recoverReferencePhaseFromEntryCartesian(entryCartesian, curvature);
        const double s      = phi / curvature;
        const double radialDistance =
                Kokkos::hypot(entryCartesian(0) + radius, entryCartesian(2)) - Kokkos::abs(radius);
        return Vector_t<double, 3>(radialDistance, entryCartesian(1), s);
    }

    /**
     * @brief Convert an RBEND entry-frame point into the OPAL rectangular field chart.
     *
     * Old OPAL's `Bend2D` treats a rectangular bend as a circular central
     * region with entrance and exit fringe regions measured perpendicular to
     * the physical pole faces.  The circular reference orbit exits upstream of
     * the rectangular exit face for nonzero bend angle, so using the circular
     * exit tangent inside the exit fringe makes the field decay too early and
     * introduces a branch discontinuity at the rectangular exit plane.  This
     * helper switches to the rectangular exit-face coordinate throughout the
     * exit fringe support.
     */
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> convertRBendEntryCartesianToFieldLocalDevice(
            const Vector_t<double, 3>& entryCartesian, const double referenceLength,
            const double curvature, const double chordLength, const double bendAngle,
            const double exitFringe) {
        if (Kokkos::abs(curvature) <= 1.0e-15) {
            return entryCartesian;
        }
        if (entryCartesian(2) <= 0.0) {
            return entryCartesian;
        }

        const double halfAngle     = 0.5 * bendAngle;
        const double exitOriginX   = -chordLength * Kokkos::tan(halfAngle);
        const double exitOriginZ   = chordLength;
        const double dx            = entryCartesian(0) - exitOriginX;
        const double dz            = entryCartesian(2) - exitOriginZ;
        const double cosAngle      = Kokkos::cos(bendAngle);
        const double sinAngle      = Kokkos::sin(bendAngle);
        const double exitLocalX    = cosAngle * dx + sinAngle * dz;
        const double exitLocalZ    = -sinAngle * dx + cosAngle * dz;
        const double exitThreshold = (exitFringe > 0.0) ? -exitFringe : 0.0;
        if (exitLocalZ >= exitThreshold) {
            return Vector_t<double, 3>(exitLocalX, entryCartesian(1), referenceLength + exitLocalZ);
        }

        return convertEntryCartesianToFieldLocalDevice(entryCartesian, referenceLength, curvature);
    }

    KOKKOS_INLINE_FUNCTION double clampUnitInterval(const double value) {
        if (value < 0.0) {
            return 0.0;
        }
        if (value > 1.0) {
            return 1.0;
        }
        return value;
    }

    struct FringeFieldScale {
        double scale;
        double derivative;
        double secondDerivative;
    };

    KOKKOS_INLINE_FUNCTION double safeAbsCos(const double angle) {
        const double cosAngle = Kokkos::abs(Kokkos::cos(angle));
        return (cosAngle > 1.0e-6) ? cosAngle : 1.0e-6;
    }

    KOKKOS_INLINE_FUNCTION double resolveFringeProfileGap(
            const double fullGap, const double fringeHalfGap) {
        if (fullGap > 0.0) {
            return fullGap;
        }
        return 2.0 * fringeHalfGap;
    }

    KOKKOS_INLINE_FUNCTION double defaultOpalFringeHalfWidth(const double profileGap) {
        // OPAL's FM1DProfile1 default map uses p1=-0.1 m, p2=0, p3=0.1 m
        // for a 0.02 m full gap. Bend2D scales these distances by GAP/0.02.
        return 5.0 * profileGap;
    }

    KOKKOS_INLINE_FUNCTION FringeFieldScale
    evaluateDefaultOpalEngeProfile(const double longitudinalCoordinate, const double profileGap) {
        if (profileGap <= 0.0) {
            return FringeFieldScale{1.0, 0.0, 0.0};
        }

        const double u  = longitudinalCoordinate / profileGap;
        const double u2 = u * u;
        const double u3 = u2 * u;
        const double u4 = u2 * u2;
        const double u5 = u4 * u;

        const double a0 = 0.478959;
        const double a1 = 1.911289;
        const double a2 = -1.185953;
        const double a3 = 1.630554;
        const double a4 = -1.082657;
        const double a5 = 0.318111;

        const double exponent = a0 + a1 * u + a2 * u2 + a3 * u3 + a4 * u4 + a5 * u5;
        if (exponent > 80.0) {
            return FringeFieldScale{0.0, 0.0, 0.0};
        }
        if (exponent < -80.0) {
            return FringeFieldScale{1.0, 0.0, 0.0};
        }

        const double invGap = 1.0 / profileGap;
        const double polynomialPrime =
                a1 + 2.0 * a2 * u + 3.0 * a3 * u2 + 4.0 * a4 * u3 + 5.0 * a5 * u4;
        const double polynomialSecond   = 2.0 * a2 + 6.0 * a3 * u + 12.0 * a4 * u2 + 20.0 * a5 * u3;
        const double exponentDerivative = polynomialPrime * invGap;
        const double exponentSecond     = polynomialSecond * invGap * invGap;
        const double engeExp            = Kokkos::exp(exponent);
        const double scale              = 1.0 / (1.0 + engeExp);
        const double scaleSquared       = scale * scale;
        const double expDerivative      = exponentDerivative * engeExp;
        const double expSecond =
                (exponentSecond + exponentDerivative * exponentDerivative) * engeExp;
        const double derivative = -expDerivative * scaleSquared;
        const double secondDerivative =
                -expSecond * scaleSquared
                + 2.0 * expDerivative * expDerivative * scaleSquared * scale;

        return FringeFieldScale{scale, derivative, secondDerivative};
    }

    KOKKOS_INLINE_FUNCTION FringeFieldScale evaluateFieldScale(
            const double z, const double fieldBegin, const double bodyLength, const double fieldEnd,
            const double profileGap) {
        const FringeFieldScale zero{0.0, 0.0, 0.0};
        if (z < fieldBegin || z > fieldEnd) {
            return zero;
        }

        const double entryFringe = -fieldBegin;
        const double exitFringe  = fieldEnd - bodyLength;
        if (profileGap <= 0.0 || (entryFringe <= 0.0 && exitFringe <= 0.0)) {
            return FringeFieldScale{1.0, 0.0, 0.0};
        }

        const double fringeHalfWidth = defaultOpalFringeHalfWidth(profileGap);
        FringeFieldScale field{1.0, 0.0, 0.0};

        if (entryFringe > 0.0 && z < entryFringe) {
            const double entryProjection =
                    (fringeHalfWidth > 0.0) ? fringeHalfWidth / entryFringe : 1.0;
            const double coordinate      = -z * entryProjection;
            const FringeFieldScale entry = evaluateDefaultOpalEngeProfile(coordinate, profileGap);
            field.scale                  = clampUnitInterval(entry.scale);
            field.derivative             = -entryProjection * entry.derivative;
            field.secondDerivative = entryProjection * entryProjection * entry.secondDerivative;
        }

        if (exitFringe > 0.0 && z > bodyLength - exitFringe) {
            const double exitProjection =
                    (fringeHalfWidth > 0.0) ? fringeHalfWidth / exitFringe : 1.0;
            const double coordinate     = (z - bodyLength) * exitProjection;
            const FringeFieldScale exit = evaluateDefaultOpalEngeProfile(coordinate, profileGap);
            const double exitScale      = clampUnitInterval(exit.scale);
            if (exitScale < field.scale) {
                field.scale            = exitScale;
                field.derivative       = exitProjection * exit.derivative;
                field.secondDerivative = exitProjection * exitProjection * exit.secondDerivative;
            }
        }

        return field;
    }

    double integrateFieldScale(
            const double fieldBegin, const double bodyLength, const double fieldEnd,
            const double profileGap) {
        if (fieldEnd <= fieldBegin) {
            return 0.0;
        }

        const double span     = fieldEnd - fieldBegin;
        const double stepGoal = (profileGap > 0.0) ? 0.02 * profileGap : 1.0e-3;
        int intervals         = std::max(256, static_cast<int>(std::ceil(span / stepGoal)));
        intervals             = std::min(intervals, 16384);
        intervals += intervals % 2;
        const double dz = span / static_cast<double>(intervals);

        double integral =
                evaluateFieldScale(fieldBegin, fieldBegin, bodyLength, fieldEnd, profileGap).scale
                + evaluateFieldScale(fieldEnd, fieldBegin, bodyLength, fieldEnd, profileGap).scale;
        for (int i = 1; i < intervals; ++i) {
            const double z      = fieldBegin + static_cast<double>(i) * dz;
            const double weight = (i % 2 == 0) ? 2.0 : 4.0;
            integral += weight
                        * evaluateFieldScale(z, fieldBegin, bodyLength, fieldEnd, profileGap).scale;
        }

        return integral * dz / 3.0;
    }

    Euclid3D makeReferencePathTransformFromEntry(
            const double s, const double bodyLength, const double curvature) {
        Euclid3D transform;
        if (std::abs(curvature) <= 1.0e-15) {
            transform.setZ(s);
            return transform;
        }

        if (s <= 0.0) {
            transform.setZ(s);
            return transform;
        }

        if (s >= bodyLength) {
            const double exitPhi = curvature * bodyLength;
            transform            = Euclid3D::YRotation(-exitPhi);
            transform.setX(
                    (std::cos(exitPhi) - 1.0) / curvature - std::sin(exitPhi) * (s - bodyLength));
            transform.setZ(std::sin(exitPhi) / curvature + std::cos(exitPhi) * (s - bodyLength));
            return transform;
        }

        const double phi = curvature * s;
        transform        = Euclid3D::YRotation(-phi);
        transform.setX((std::cos(phi) - 1.0) / curvature);
        transform.setZ(std::sin(phi) / curvature);
        return transform;
    }

}  // namespace

BendBase::BendBase() : BendBase("") {}

BendBase::BendBase(const BendBase& right)
    : Component(right),
      startField_m(right.startField_m),
      endField_m(right.endField_m),
      fieldBegin_m(right.fieldBegin_m),
      fieldEnd_m(right.fieldEnd_m),
      effectiveFieldLength_m(right.effectiveFieldLength_m),
      angle_m(right.angle_m),
      entranceAngle_m(right.entranceAngle_m),
      exitAngle_m(right.exitAngle_m),
      gap_m(right.gap_m),
      fringeHalfGap_m(right.fringeHalfGap_m),
      fringeIntegral_m(right.fringeIntegral_m),
      designEnergy_m(right.designEnergy_m),
      designEnergyChangeable_m(true),
      fieldAmplitudeX_m(right.fieldAmplitudeX_m),
      fieldAmplitudeY_m(right.fieldAmplitudeY_m),
      fieldAmplitude_m(right.fieldAmplitude_m),
      normalizedField_m(right.normalizedField_m),
      hasNormalizedField_m(right.hasNormalizedField_m),
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
      fieldBegin_m(0.0),
      fieldEnd_m(0.0),
      effectiveFieldLength_m(0.0),
      angle_m(0.0),
      entranceAngle_m(0.0),
      exitAngle_m(0.0),
      gap_m(0.0),
      fringeHalfGap_m(0.0),
      fringeIntegral_m(0.5),
      designEnergy_m(0.0),
      designEnergyChangeable_m(true),
      fieldAmplitudeX_m(0.0),
      fieldAmplitudeY_m(0.0),
      fieldAmplitude_m(0.0),
      normalizedField_m(),
      hasNormalizedField_m(false),
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
    updateFieldSupportExtent();
    updatePhysicalFieldFromReference();
    startField_m = startField + fieldBegin_m;
    endField_m   = startField + fieldEnd_m;
    endField     = endField_m;
    online_m     = true;
}

void BendBase::finalise() { online_m = false; }

CoordinateSystemTrafo BendBase::getFieldCSTrafoLab2Local(const PlacedElement& placed) const {
    return placed.getNominalEntryTransform();
}

Vector_t<double, 3> BendBase::transformFieldFrameToLocal(const Vector_t<double, 3>& r) const {
    return convertEntryCartesianToFieldLocal(r);
}

Vector_t<double, 3> BendBase::rotateFieldFrameToLocal(
        const Vector_t<double, 3>& v, const Vector_t<double, 3>& fieldLocalPosition) const {
    return rotateEntryCartesianVectorToFieldLocal(v, fieldLocalPosition(2));
}

Vector_t<double, 3> BendBase::rotateFieldLocalToFieldFrame(
        const Vector_t<double, 3>& v, const Vector_t<double, 3>& fieldLocalPosition) const {
    return rotateFieldLocalVectorToEntryCartesian(v, fieldLocalPosition(2));
}

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

    const double fieldBegin                   = fieldBegin_m;
    const double elemLength                   = getReferencePathLength();
    const double fieldEnd                     = fieldEnd_m;
    const double entryFringe                  = -fieldBegin;
    const double exitFringe                   = fieldEnd - elemLength;
    const double profileGap                   = resolveFringeProfileGap(gap_m, fringeHalfGap_m);
    const double entryEdgeVerticalCoefficient = getEntryEdgeVerticalFieldCoefficient();
    const double exitEdgeVerticalCoefficient  = getExitEdgeVerticalFieldCoefficient();

    Kokkos::parallel_for(
            "BendBase::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const double z = Rview(i)(2);
                if (z < fieldBegin || z > fieldEnd) {
                    return;
                }

                Vector_t<double, 3> Bf(0.0);
                const double x = Rview(i)(0);
                const double y = Rview(i)(1);
                const FringeFieldScale fringe =
                        evaluateFieldScale(z, fieldBegin, elemLength, fieldEnd, profileGap);
                const double scale = fringe.scale;

                if (normal.extent(0) > 0) {
                    Bf(1) += normal(0) * (scale - 0.5 * fringe.secondDerivative * y * y);
                    Bf(2) += normal(0) * fringe.derivative * y;
                }
                if (skew.extent(0) > 0) {
                    Bf(0) -= scale * skew(0);
                }
                if (normal.extent(0) > 1) {
                    Bf(0) += scale * normal(1) * y;
                    Bf(1) += scale * normal(1) * x;
                }
                if (skew.extent(0) > 1) {
                    Bf(0) -= scale * skew(1) * x;
                    Bf(1) += scale * skew(1) * y;
                }

                double verticalGradient = 0.0;
                if (z <= entryFringe && entryFringe > 0.0) {
                    // Weight the vertical edge force with the same Enge derivative as OPAL's
                    // fringe map; this preserves the integrated kick and its longitudinal centroid.
                    verticalGradient = entryEdgeVerticalCoefficient * fringe.derivative;
                } else if (z >= elemLength - exitFringe && exitFringe > 0.0) {
                    verticalGradient = -exitEdgeVerticalCoefficient * fringe.derivative;
                }
                if (verticalGradient != 0.0) {
                    Bf(0) += verticalGradient * y;
                }

                for (unsigned d = 0; d < 3; ++d) {
                    Eview(i)(d) += 0.0;
                    Bview(i)(d) += Bf(d);
                }
            });

    return false;
}

bool BendBase::applyToBunch(
        const std::shared_ptr<ParticleContainer_t>& pc,
        const CoordinateSystemTrafo& refToFieldCSTrafo) {
    const auto slices = buildTrackingSlices();
    for (const auto& slice : slices) {
        const CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * refToFieldCSTrafo;
        const CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
        pc->transformBunch(refToSliceLocal);
        applySlice(pc, slice);
        pc->transformBunch(sliceLocalToRef);
    }
    return false;
}

std::vector<BendBase::TrackingSlice> BendBase::buildTrackingSlices() const {
    const double bodyLength  = getReferencePathLength();
    const double sBegin      = -getEntryFringeSupportLength();
    const double sEnd        = bodyLength + getExitFringeSupportLength();
    const double totalLength = std::max(0.0, sEnd - sBegin);
    std::size_t nSlices      = 1u;
    if (getStepsize() > 1.0e-15) {
        nSlices = std::max<std::size_t>(
                1u, static_cast<std::size_t>(std::ceil(totalLength / getStepsize())));
    } else if (getNSlices() > 1) {
        nSlices = getNSlices();
    } else if (getSlices() > 1.0) {
        nSlices = std::max<std::size_t>(1u, static_cast<std::size_t>(std::ceil(getSlices())));
    } else {
        const std::size_t lengthDriven =
                std::max<std::size_t>(1u, static_cast<std::size_t>(std::ceil(totalLength / 0.05)));
        const std::size_t angleDriven = std::max<std::size_t>(
                1u, static_cast<std::size_t>(
                            std::ceil(std::abs(getBendAngle()) / (Physics::pi / 32.0))));
        nSlices = std::max<std::size_t>(8u, std::max(lengthDriven, angleDriven));
    }
    const double ds        = (nSlices > 0) ? (sEnd - sBegin) / static_cast<double>(nSlices) : 0.0;
    const double curvature = getReferencePathCurvature();
    const CoordinateSystemTrafo entryToBodyLocal = toCoordinateSystemTrafo(
            makeReferencePathTransformFromEntry(0.5 * bodyLength, bodyLength, curvature));
    const CoordinateSystemTrafo bodyToEntryLocal = entryToBodyLocal.inverted();

    std::vector<TrackingSlice> slices;
    slices.reserve(nSlices);
    for (std::size_t i = 0; i < nSlices; ++i) {
        const double sliceBegin                       = sBegin + static_cast<double>(i) * ds;
        const double sliceEnd                         = sBegin + static_cast<double>(i + 1) * ds;
        const double sliceCenter                      = 0.5 * (sliceBegin + sliceEnd);
        const CoordinateSystemTrafo entryToSliceLocal = toCoordinateSystemTrafo(
                makeReferencePathTransformFromEntry(sliceCenter, bodyLength, curvature));
        const CoordinateSystemTrafo bodyToSliceLocal = entryToSliceLocal * bodyToEntryLocal;
        slices.push_back(
                TrackingSlice{
                        entryToSliceLocal, bodyToSliceLocal, entryToSliceLocal.getRotationMatrix(),
                        entryToSliceLocal.getOrigin(), sliceBegin, sliceCenter, sliceEnd});
    }
    return slices;
}

Vector_t<double, 3> BendBase::convertEntryCartesianToFieldLocal(
        const Vector_t<double, 3>& entryCartesian) const {
    const double bodyLength = getReferencePathLength();
    const double curvature  = getReferencePathCurvature();
    if (std::abs(curvature) <= 1.0e-15) {
        return entryCartesian;
    }

    if (entryCartesian(2) <= 0.0) {
        return entryCartesian;
    }

    if (getType() == ElementType::RBEND) {
        return convertRBendEntryCartesianToFieldLocalDevice(
                entryCartesian, bodyLength, curvature, getElementLength(), getBendAngle(),
                getExitFringeSupportLength());
    }

    const CoordinateSystemTrafo entryToExitLocal = toCoordinateSystemTrafo(
            makeReferencePathTransformFromEntry(bodyLength, bodyLength, curvature));
    const Vector_t<double, 3> exitLocal = entryToExitLocal.transformTo(entryCartesian);
    if (exitLocal(2) >= 0.0) {
        return Vector_t<double, 3>(exitLocal(0), exitLocal(1), bodyLength + exitLocal(2));
    }

    const double radius = 1.0 / curvature;
    const double phi    = recoverReferencePhaseFromEntryCartesian(entryCartesian, curvature);
    const double s      = phi / curvature;
    const double radialDistance =
            std::hypot(entryCartesian(0) + radius, entryCartesian(2)) - std::abs(radius);
    return Vector_t<double, 3>(radialDistance, entryCartesian(1), s);
}

Vector_t<double, 3> BendBase::convertBodyCartesianToFieldLocal(
        const Vector_t<double, 3>& bodyCartesian) const {
    const double bodyLength = getReferencePathLength();
    const double curvature  = getReferencePathCurvature();

    if (getType() == ElementType::RBEND) {
        const CoordinateSystemTrafo entryLocal   = getEdgeToBegin();
        const Vector_t<double, 3> entryCartesian = entryLocal.transformTo(bodyCartesian);
        if (std::abs(curvature) <= 1.0e-15 || entryCartesian(2) <= 0.0) {
            return entryCartesian;
        }

        const CoordinateSystemTrafo exitLocal   = getEdgeToEnd();
        const Vector_t<double, 3> exitCartesian = exitLocal.transformTo(bodyCartesian);
        const double exitFringe                 = getExitFringeSupportLength();
        const double exitThreshold              = (exitFringe > 0.0) ? -exitFringe : 0.0;
        if (exitCartesian(2) >= exitThreshold) {
            return Vector_t<double, 3>(
                    exitCartesian(0), exitCartesian(1), bodyLength + exitCartesian(2));
        }

        /**
         * For an RBEND the hardware body remains straight while the design
         * reference path is curved.  The entrance fringe and central region use
         * the entry-based circular chart, while the exit fringe uses the
         * rectangular exit-face chart.  The physical pole face sits downstream
         * of the circular reference-orbit exit, so using the exit-face chart for
         * the whole exit-fringe support avoids an artificial branch switch at
         * the rectangular exit plane:
         * \f[
         *   \mathbf r_\mathrm{entry}
         *     = T_{b\rightarrow e}\,\mathbf r_\mathrm{body}, \qquad
         *   \mathbf r_\mathrm{exit}
         *     = T_{b\rightarrow x}\,\mathbf r_\mathrm{body}, \qquad
         *   s = L_\mathrm{ref} + z_\mathrm{exit}
         * \f]
         * for \f$z_\mathrm{exit}\ge -\ell_\mathrm{fringe,exit}\f$.  This keeps
         * upstream and downstream drift extensions outside the finite field
         * support instead of projecting them through the straight body midpoint.
         */
        return convertEntryCartesianToFieldLocal(entryCartesian);
    }

    const CoordinateSystemTrafo entryLocal   = toCoordinateSystemTrafo(getEntranceFrame());
    const Vector_t<double, 3> entryCartesian = entryLocal.transformTo(bodyCartesian);
    if (std::abs(curvature) <= 1.0e-15 || entryCartesian(2) <= 0.0) {
        return entryCartesian;
    }

    const CoordinateSystemTrafo exitLocal   = toCoordinateSystemTrafo(getExitFrame());
    const Vector_t<double, 3> exitCartesian = exitLocal.transformTo(bodyCartesian);
    if (exitCartesian(2) >= 0.0) {
        return Vector_t<double, 3>(
                exitCartesian(0), exitCartesian(1), bodyLength + exitCartesian(2));
    }

    return convertEntryCartesianToFieldLocal(entryCartesian);
}

Vector_t<double, 3> BendBase::rotateBodyCartesianVectorToFieldLocal(
        const Vector_t<double, 3>& bodyVector, const double s) const {
    const double bodyLength                = getReferencePathLength();
    const double curvature                 = getReferencePathCurvature();
    const CoordinateSystemTrafo entryLocal = (getType() == ElementType::RBEND)
                                                     ? getEdgeToBegin()
                                                     : toCoordinateSystemTrafo(getEntranceFrame());
    const Vector_t<double, 3> entryVector  = entryLocal.rotateTo(bodyVector);
    if (std::abs(curvature) <= 1.0e-15 || s <= 0.0) {
        return entryVector;
    }

    const CoordinateSystemTrafo exitLocal = (getType() == ElementType::RBEND)
                                                    ? getEdgeToEnd()
                                                    : toCoordinateSystemTrafo(getExitFrame());
    const double exitThreshold =
            (getType() == ElementType::RBEND && getExitFringeSupportLength() > 0.0)
                    ? bodyLength - getExitFringeSupportLength()
                    : bodyLength;
    if (s >= exitThreshold) {
        return exitLocal.rotateTo(bodyVector);
    }

    return rotateEntryCartesianVectorToFieldLocal(entryVector, s);
}

Vector_t<double, 3> BendBase::rotateEntryCartesianVectorToFieldLocal(
        const Vector_t<double, 3>& entryVector, const double s) const {
    const double bodyLength = getReferencePathLength();
    const double curvature  = getReferencePathCurvature();
    if (std::abs(curvature) <= 1.0e-15) {
        return entryVector;
    }

    if (s <= 0.0) {
        return entryVector;
    }

    if (s >= bodyLength) {
        const CoordinateSystemTrafo entryToExitLocal = toCoordinateSystemTrafo(
                makeReferencePathTransformFromEntry(bodyLength, bodyLength, curvature));
        return entryToExitLocal.rotateTo(entryVector);
    }

    const double phi  = curvature * s;
    const double cphi = std::cos(phi);
    const double sphi = std::sin(phi);
    return Vector_t<double, 3>(
            cphi * entryVector(0) + sphi * entryVector(2), entryVector(1),
            -sphi * entryVector(0) + cphi * entryVector(2));
}

Vector_t<double, 3> BendBase::rotateFieldLocalVectorToEntryCartesian(
        const Vector_t<double, 3>& fieldVector, const double s) const {
    const double bodyLength = getReferencePathLength();
    const double curvature  = getReferencePathCurvature();
    if (std::abs(curvature) <= 1.0e-15) {
        return fieldVector;
    }

    if (s <= 0.0) {
        return fieldVector;
    }

    if (s >= bodyLength) {
        const CoordinateSystemTrafo entryToExitLocal = toCoordinateSystemTrafo(
                makeReferencePathTransformFromEntry(bodyLength, bodyLength, curvature));
        return entryToExitLocal.rotateFrom(fieldVector);
    }

    const double phi  = curvature * s;
    const double cphi = std::cos(phi);
    const double sphi = std::sin(phi);
    return Vector_t<double, 3>(
            cphi * fieldVector(0) - sphi * fieldVector(2), fieldVector(1),
            sphi * fieldVector(0) + cphi * fieldVector(2));
}

bool BendBase::applySlice(
        const std::shared_ptr<ParticleContainer_t>& pc, const TrackingSlice& slice) const {
    auto Rview          = pc->R.getView();
    auto Eview          = pc->E.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    const BMultipoleField& field = getField();
    const int order              = field.order();
    Kokkos::View<double*> normal("BendBase::sliceNormal", order);
    Kokkos::View<double*> skew("BendBase::sliceSkew", order);
    auto normalHost = Kokkos::create_mirror_view(normal);
    auto skewHost   = Kokkos::create_mirror_view(skew);
    for (int i = 0; i < order; ++i) {
        normalHost(i) = field.getNormalComponent(i);
        skewHost(i)   = field.getSkewComponent(i);
    }
    Kokkos::deep_copy(normal, normalHost);
    Kokkos::deep_copy(skew, skewHost);

    const double bodyLength                   = getReferencePathLength();
    const double fieldBegin                   = -getEntryFringeSupportLength();
    const double fieldEnd                     = bodyLength + getExitFringeSupportLength();
    const double entryFringe                  = getEntryFringeSupportLength();
    const double exitFringe                   = getExitFringeSupportLength();
    const double profileGap                   = resolveFringeProfileGap(gap_m, fringeHalfGap_m);
    const double referenceCurvature           = getReferencePathCurvature();
    const double chordLength                  = getElementLength();
    const double bendAngle                    = getBendAngle();
    const bool isRectangularBend              = getType() == ElementType::RBEND;
    const double entryEdgeVerticalCoefficient = getEntryEdgeVerticalFieldCoefficient();
    const double exitEdgeVerticalCoefficient  = getExitEdgeVerticalFieldCoefficient();
    const matrix3x3_t entryToSliceRotation    = slice.entryToSliceRotation;
    const Vector_t<double, 3> entryOrigin     = slice.entryOrigin;

    Kokkos::parallel_for(
            "BendBase::applySlice", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const Vector_t<double, 3> entryCartesian =
                        prod_vector_transpose(entryToSliceRotation, Rview(i)) + entryOrigin;
                const Vector_t<double, 3> entryChartPosition =
                        isRectangularBend ? convertRBendEntryCartesianToFieldLocalDevice(
                                                    entryCartesian, bodyLength, referenceCurvature,
                                                    chordLength, bendAngle, exitFringe)
                                          : convertEntryCartesianToFieldLocalDevice(
                                                    entryCartesian, bodyLength, referenceCurvature);

                const double zEntry = entryChartPosition(2);
                if (zEntry < slice.sBegin || zEntry > slice.sEnd) {
                    return;
                }
                if (zEntry < fieldBegin || zEntry > fieldEnd) {
                    return;
                }

                Vector_t<double, 3> Bf(0.0);
                const double x = entryChartPosition(0);
                const double y = entryChartPosition(1);
                const FringeFieldScale fringe =
                        evaluateFieldScale(zEntry, fieldBegin, bodyLength, fieldEnd, profileGap);
                const double scale = fringe.scale;

                if (normal.extent(0) > 0) {
                    Bf(1) += normal(0) * (scale - 0.5 * fringe.secondDerivative * y * y);
                    Bf(2) += normal(0) * fringe.derivative * y;
                }
                if (skew.extent(0) > 0) {
                    Bf(0) -= scale * skew(0);
                }
                if (normal.extent(0) > 1) {
                    Bf(0) += scale * normal(1) * y;
                    Bf(1) += scale * normal(1) * x;
                }
                if (skew.extent(0) > 1) {
                    Bf(0) -= scale * skew(1) * x;
                    Bf(1) += scale * skew(1) * y;
                }

                double verticalGradient = 0.0;
                if (zEntry <= entryFringe && entryFringe > 0.0) {
                    // Weight the vertical edge force with the same Enge derivative as OPAL's
                    // fringe map; this preserves the integrated kick and its longitudinal centroid.
                    verticalGradient = entryEdgeVerticalCoefficient * fringe.derivative;
                } else if (zEntry >= bodyLength - exitFringe && exitFringe > 0.0) {
                    verticalGradient = -exitEdgeVerticalCoefficient * fringe.derivative;
                }
                if (verticalGradient != 0.0) {
                    Bf(0) += verticalGradient * y;
                }

                Vector_t<double, 3> Bentry = Bf;
                if (isRectangularBend && zEntry >= bodyLength - exitFringe && exitFringe > 0.0) {
                    const double cosAngle = Kokkos::cos(bendAngle);
                    const double sinAngle = Kokkos::sin(bendAngle);
                    Bentry(0)             = cosAngle * Bf(0) - sinAngle * Bf(2);
                    Bentry(1)             = Bf(1);
                    Bentry(2)             = sinAngle * Bf(0) + cosAngle * Bf(2);
                } else if (Kokkos::abs(referenceCurvature) > 1.0e-15 && zEntry > 0.0) {
                    const double rotationS = (zEntry >= bodyLength) ? bodyLength : zEntry;
                    const double phi       = referenceCurvature * rotationS;
                    const double cphi      = Kokkos::cos(phi);
                    const double sphi      = Kokkos::sin(phi);
                    Bentry(0)              = cphi * Bf(0) - sphi * Bf(2);
                    Bentry(1)              = Bf(1);
                    Bentry(2)              = sphi * Bf(0) + cphi * Bf(2);
                }

                const Vector_t<double, 3> Bslice = prod_vector(entryToSliceRotation, Bentry);

                for (unsigned d = 0; d < 3; ++d) {
                    Eview(i)(d) += 0.0;
                    Bview(i)(d) += Bslice(d);
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

    const FringeFieldScale fringe = evaluateFieldScale(
            R(2), fieldBegin_m, getReferencePathLength(), fieldEnd_m,
            resolveFringeProfileGap(gap_m, fringeHalfGap_m));
    computeFieldHost(R, getField(), B);
    B *= fringe.scale;
    if (getField().order() > 0) {
        const double dipole = getField().getNormalComponent(0);
        B(1) += -0.5 * dipole * fringe.secondDerivative * R(1) * R(1);
        B(2) += dipole * fringe.derivative * R(1);
    }
    const double entryFringe     = getEntryFringeSupportLength();
    const double exitFringe      = getExitFringeSupportLength();
    const double referenceLength = getReferencePathLength();
    if (R(2) <= entryFringe && entryFringe > 0.0) {
        B(0) += getEntryEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    } else if (R(2) >= referenceLength - exitFringe && exitFringe > 0.0) {
        B(0) += -getExitEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    }
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

    const FringeFieldScale fringe = evaluateFieldScale(
            R(2), fieldBegin_m, getReferencePathLength(), fieldEnd_m,
            resolveFringeProfileGap(gap_m, fringeHalfGap_m));
    computeFieldHost(R, getField(), B);
    B *= fringe.scale;
    if (getField().order() > 0) {
        const double dipole = getField().getNormalComponent(0);
        B(1) += -0.5 * dipole * fringe.secondDerivative * R(1) * R(1);
        B(2) += dipole * fringe.derivative * R(1);
    }
    const double entryFringe     = getEntryFringeSupportLength();
    const double exitFringe      = getExitFringeSupportLength();
    const double referenceLength = getReferencePathLength();
    if (R(2) <= entryFringe && entryFringe > 0.0) {
        B(0) += getEntryEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    } else if (R(2) >= referenceLength - exitFringe && exitFringe > 0.0) {
        B(0) += -getExitEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    }
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

    const FringeFieldScale fringe = evaluateFieldScale(
            R(2), fieldBegin_m, getReferencePathLength(), fieldEnd_m,
            resolveFringeProfileGap(gap_m, fringeHalfGap_m));
    computeFieldHost(R, getField(), B);
    B *= fringe.scale;
    if (getField().order() > 0) {
        const double dipole = getField().getNormalComponent(0);
        B(1) += -0.5 * dipole * fringe.secondDerivative * R(1) * R(1);
        B(2) += dipole * fringe.derivative * R(1);
    }
    const double entryFringe     = getEntryFringeSupportLength();
    const double exitFringe      = getExitFringeSupportLength();
    const double referenceLength = getReferencePathLength();
    if (R(2) <= entryFringe && entryFringe > 0.0) {
        B(0) += getEntryEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    } else if (R(2) >= referenceLength - exitFringe && exitFringe > 0.0) {
        B(0) += -getExitEdgeVerticalFieldCoefficient() * fringe.derivative * R(1);
    }
    (void)E;
    return false;
}

void BendBase::getFieldExtend(double& zBegin, double& zEnd) const {
    zBegin = fieldBegin_m;
    zEnd   = fieldEnd_m;
}

bool BendBase::isInside(const Vector_t<double, 3>& r) const {
    return r(2) >= fieldBegin_m && r(2) < fieldEnd_m && isInsideTransverse(r);
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

Vector_t<double, 3> BendBase::getDesignPathPoint(const double z) const {
    const double sBegin = getEntrance();
    const double sEnd   = getExit();
    const double s      = sBegin + z;

    const auto pointFromPose = [](const Euclid3D& pose) {
        return Vector_t<double, 3>(pose.getX(), pose.getY(), pose.getZ());
    };

    if (s < sBegin) {
        const Euclid3D entryPose = getTransform(sBegin);
        const Vector_t<double, 3> entry =
                Vector_t<double, 3>(entryPose.getX(), entryPose.getY(), entryPose.getZ());
        const Vector3D entryTangent = entryPose.getRotation() * Vector3D(0.0, 0.0, 1.0);
        const Vector_t<double, 3> tangent =
                Vector_t<double, 3>(entryTangent.getX(), entryTangent.getY(), entryTangent.getZ());
        return entry + (s - sBegin) * tangent;
    }

    if (s > sEnd) {
        const Euclid3D exitPose = getTransform(sEnd);
        const Vector_t<double, 3> exit =
                Vector_t<double, 3>(exitPose.getX(), exitPose.getY(), exitPose.getZ());
        const Vector3D exitTangent = exitPose.getRotation() * Vector3D(0.0, 0.0, 1.0);
        const Vector_t<double, 3> tangent =
                Vector_t<double, 3>(exitTangent.getX(), exitTangent.getY(), exitTangent.getZ());
        return exit + (s - sEnd) * tangent;
    }

    return pointFromPose(getTransform(s));
}

double BendBase::calcDesignRadius(double fieldAmplitude) const {
    const auto& reference  = *RefPartBunch_m->getParticleContainer()->getReference();
    const double mass      = reference.getM();
    const double momentum  = getReferenceMomentumEV();
    const double betaGamma = momentum / mass;
    const double charge    = reference.getQ();
    return std::abs(betaGamma * mass / (Physics::c * fieldAmplitude * charge));
}

double BendBase::calcFieldAmplitude(double radius) const {
    const auto& reference  = *RefPartBunch_m->getParticleContainer()->getReference();
    const double mass      = reference.getM();
    const double momentum  = getReferenceMomentumEV();
    const double betaGamma = momentum / mass;
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
    if (designEnergy_m > 0.0) {
        return designEnergy_m / mass + 1.0;
    }

    const double betaGamma = reference.getP() / mass;
    return std::sqrt(1.0 + betaGamma * betaGamma);
}

double BendBase::calcBetaGamma() const {
    const double gamma = calcGamma();
    return std::sqrt(gamma * gamma - 1.0);
}

double BendBase::getReferenceMomentumEV() const {
    const auto container  = RefPartBunch_m ? RefPartBunch_m->getParticleContainer()
                                           : std::shared_ptr<ParticleContainer_t>();
    const auto* reference = container ? container->getReference() : nullptr;
    if (reference == nullptr) {
        return 0.0;
    }

    if (designEnergy_m <= 0.0 && container->getTotalNum() > 0) {
        const Vector_t<double, 3> meanP = container->getMeanP();
        return std::sqrt(dot(meanP, meanP)) * reference->getM();
    }

    return resolveReferenceMomentumEV(reference->getP(), reference->getM());
}

void BendBase::updatePhysicalFieldFromMomentumEV(
        const double referenceMomentumEV, const double charge) {
    if (!hasNormalizedField_m || referenceMomentumEV <= 0.0 || std::abs(charge) <= 1.0e-15) {
        return;
    }

    // OPALX's bend reference geometry is defined so that a positive horizontal
    // bend angle deflects the design trajectory toward negative entry-frame x.
    // The required dipole sign therefore satisfies
    // \f[
    // q v_z B_y = \frac{p_0 c}{\rho},
    // \qquad
    // B_y = \frac{p_0}{q c} h,
    // \f]
    // with signed curvature \f$h = 1/\rho\f$. This sign convention is shared
    // by both `SBEND` and `RBEND` runtime normalization and is verified by the
    // reference-particle bend regression tests expressed in the entry/lab
    // geometry chart for positive and negative bend angles.
    const double factor = referenceMomentumEV / (charge * Physics::c);
    BMultipoleField field;
    const int order = normalizedField_m.order();
    for (int i = 0; i < order; ++i) {
        field.setNormalComponent(i, factor * normalizedField_m.getNormalComponent(i));
        field.setSkewComponent(i, factor * normalizedField_m.getSkewComponent(i));
    }
    getField() = field;
}

void BendBase::updatePhysicalFieldFromReference() {
    const auto container  = RefPartBunch_m ? RefPartBunch_m->getParticleContainer()
                                           : std::shared_ptr<ParticleContainer_t>();
    const auto* reference = container ? container->getReference() : nullptr;
    if (reference == nullptr) {
        return;
    }
    updatePhysicalFieldFromMomentumEV(getReferenceMomentumEV(), reference->getQ());
}

double BendBase::getFieldScale(const double z) const {
    return evaluateFieldScale(
                   z, fieldBegin_m, getReferencePathLength(), fieldEnd_m,
                   resolveFringeProfileGap(gap_m, fringeHalfGap_m))
            .scale;
}

double BendBase::getReferencePathCurvature() const {
    const double referenceLength = getReferencePathLength();
    if (std::abs(referenceLength) <= 1.0e-15) {
        return 0.0;
    }
    return getBendAngle() / referenceLength;
}

double BendBase::getSignedCurvature() const {
    const double effectiveLength = getEffectiveFieldLength();
    if (std::abs(effectiveLength) <= 1.0e-15) {
        return 0.0;
    }
    return getBendAngle() / effectiveLength;
}

double BendBase::getEntryEdgeHorizontalStrength() const {
    return getSignedCurvature() * std::tan(getEntranceAngle());
}

double BendBase::getExitEdgeHorizontalStrength() const {
    return getSignedCurvature() * std::tan(getExitAngle());
}

double BendBase::getEntryEdgeVerticalStrength() const {
    const double edgeAngle = getEntranceAngle();
    const double cosEdge   = std::cos(edgeAngle);
    const double safeCos = (std::abs(cosEdge) > 1.0e-6) ? cosEdge : std::copysign(1.0e-6, cosEdge);
    const double sinEdge = std::sin(edgeAngle);
    const double psi     = getSignedCurvature() * getFringeHalfGap() * getFringeIntegral()
                       * (1.0 + sinEdge * sinEdge) / safeCos;
    return -getSignedCurvature() * std::tan(edgeAngle - psi);
}

double BendBase::getExitEdgeVerticalStrength() const {
    const double edgeAngle = getExitAngle();
    const double cosEdge   = std::cos(edgeAngle);
    const double safeCos = (std::abs(cosEdge) > 1.0e-6) ? cosEdge : std::copysign(1.0e-6, cosEdge);
    const double sinEdge = std::sin(edgeAngle);
    const double psi     = getSignedCurvature() * getFringeHalfGap() * getFringeIntegral()
                       * (1.0 + sinEdge * sinEdge) / safeCos;
    return -getSignedCurvature() * std::tan(edgeAngle - psi);
}

double BendBase::getEntryEdgeVerticalFieldCoefficient() const {
    const double entryFringe = getEntryFringeSupportLength();
    if (entryFringe <= 0.0 || std::abs(getSignedCurvature()) <= 1.0e-15) {
        return 0.0;
    }

    const double fieldScaleSpan =
            std::abs(getFieldScale(entryFringe) - getFieldScale(-entryFringe));
    if (fieldScaleSpan <= 1.0e-15) {
        return 0.0;
    }

    const double rigidityScale = getField().getNormalComponent(0) / getSignedCurvature();
    return rigidityScale * getEntryEdgeVerticalStrength() / fieldScaleSpan;
}

double BendBase::getExitEdgeVerticalFieldCoefficient() const {
    const double exitFringe = getExitFringeSupportLength();
    if (exitFringe <= 0.0 || std::abs(getSignedCurvature()) <= 1.0e-15) {
        return 0.0;
    }

    const double bodyLength     = getReferencePathLength();
    const double fieldScaleSpan = std::abs(
            getFieldScale(bodyLength - exitFringe) - getFieldScale(bodyLength + exitFringe));
    if (fieldScaleSpan <= 1.0e-15) {
        return 0.0;
    }

    const double rigidityScale = getField().getNormalComponent(0) / getSignedCurvature();
    return rigidityScale * getExitEdgeVerticalStrength() / fieldScaleSpan;
}

void BendBase::updateFieldSupportExtent() {
    const double bodyLength       = getReferencePathLength();
    const double profileGap       = resolveFringeProfileGap(gap_m, fringeHalfGap_m);
    const double fringeHalfWidth  = defaultOpalFringeHalfWidth(profileGap);
    const double entryDenominator = safeAbsCos(getEntranceAngle());
    const double exitDenominator  = safeAbsCos(getExitAngle());
    const double entryFringe = (fringeHalfWidth > 0.0) ? fringeHalfWidth / entryDenominator : 0.0;
    const double exitFringe  = (fringeHalfWidth > 0.0) ? fringeHalfWidth / exitDenominator : 0.0;
    fieldBegin_m             = -entryFringe;
    fieldEnd_m               = bodyLength + exitFringe;
    effectiveFieldLength_m = integrateFieldScale(fieldBegin_m, bodyLength, fieldEnd_m, profileGap);
    if (effectiveFieldLength_m <= 0.0) {
        effectiveFieldLength_m = bodyLength;
    }
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
