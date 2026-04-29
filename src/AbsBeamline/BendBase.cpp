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

    double evaluateFieldScale(
            const double z, const double fieldBegin, const double bodyLength,
            const double fieldEnd) {
        if (z < fieldBegin || z > fieldEnd) {
            return 0.0;
        }
        if (z <= 0.0) {
            const double entryFringe = -fieldBegin;
            if (entryFringe <= 0.0) {
                return 1.0;
            }
            return std::clamp((z - fieldBegin) / entryFringe, 0.0, 1.0);
        }
        if (z >= bodyLength) {
            const double exitFringe = fieldEnd - bodyLength;
            if (exitFringe <= 0.0) {
                return 1.0;
            }
            return std::clamp((fieldEnd - z) / exitFringe, 0.0, 1.0);
        }
        return 1.0;
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

    const double fieldBegin       = fieldBegin_m;
    const double elemLength       = getElementLength();
    const double fieldEnd         = fieldEnd_m;
    const double entryFringe      = -fieldBegin;
    const double exitFringe       = fieldEnd - elemLength;
    const double signedCurvature  = getSignedCurvature();
    const double dipoleFieldTesla = getField().getNormalComponent(0);
    const double bendRigidityScale =
            (std::abs(signedCurvature) > 1.0e-15) ? dipoleFieldTesla / signedCurvature : 0.0;
    const double entryEdgeHorizontalGradient =
            (entryFringe > 0.0) ? bendRigidityScale * getEntryEdgeHorizontalStrength() / entryFringe
                                : 0.0;
    const double exitEdgeHorizontalGradient =
            (exitFringe > 0.0) ? bendRigidityScale * getExitEdgeHorizontalStrength() / exitFringe
                               : 0.0;
    const double entryEdgeVerticalGradient =
            (entryFringe > 0.0) ? -bendRigidityScale * getEntryEdgeVerticalStrength() / entryFringe
                                : 0.0;
    const double exitEdgeVerticalGradient =
            (exitFringe > 0.0) ? -bendRigidityScale * getExitEdgeVerticalStrength() / exitFringe
                               : 0.0;

    Kokkos::parallel_for(
            "BendBase::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const double z = Rview(i)(2);
                if (z < fieldBegin || z > fieldEnd) {
                    return;
                }

                Vector_t<double, 3> Bf(0.0);
                const double x = Rview(i)(0);
                const double y = Rview(i)(1);
                double scale   = 1.0;
                if (z <= 0.0) {
                    const double entryFringe = -fieldBegin;
                    if (entryFringe > 0.0) {
                        scale = (z - fieldBegin) / entryFringe;
                        if (scale < 0.0) {
                            scale = 0.0;
                        } else if (scale > 1.0) {
                            scale = 1.0;
                        }
                    }
                } else if (z >= elemLength) {
                    const double exitFringe = fieldEnd - elemLength;
                    if (exitFringe > 0.0) {
                        scale = (fieldEnd - z) / exitFringe;
                        if (scale < 0.0) {
                            scale = 0.0;
                        } else if (scale > 1.0) {
                            scale = 1.0;
                        }
                    }
                }

                if (normal.extent(0) > 0) {
                    Bf(1) += scale * normal(0);
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

                double horizontalGradient = 0.0;
                double verticalGradient   = 0.0;
                if (z <= 0.0 && entryFringe > 0.0) {
                    horizontalGradient = entryEdgeHorizontalGradient;
                    verticalGradient   = entryEdgeVerticalGradient;
                } else if (z >= elemLength && exitFringe > 0.0) {
                    horizontalGradient = exitEdgeHorizontalGradient;
                    verticalGradient   = exitEdgeVerticalGradient;
                }
                if (horizontalGradient != 0.0 || verticalGradient != 0.0) {
                    Bf(0) += verticalGradient * y;
                    Bf(1) += horizontalGradient * x;
                }

                for (unsigned d = 0; d < 3; ++d) {
                    Eview(i)(d) += 0.0;
                    Bview(i)(d) += Bf(d);
                }
            });

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

    std::vector<TrackingSlice> slices;
    slices.reserve(nSlices);
    for (std::size_t i = 0; i < nSlices; ++i) {
        const double sliceBegin                       = sBegin + static_cast<double>(i) * ds;
        const double sliceEnd                         = sBegin + static_cast<double>(i + 1) * ds;
        const double sliceCenter                      = 0.5 * (sliceBegin + sliceEnd);
        const CoordinateSystemTrafo entryToSliceLocal = toCoordinateSystemTrafo(
                makeReferencePathTransformFromEntry(sliceCenter, bodyLength, curvature));
        slices.push_back(
                TrackingSlice{
                        entryToSliceLocal, entryToSliceLocal.getRotationMatrix(),
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

    const CoordinateSystemTrafo entryToExitLocal = toCoordinateSystemTrafo(
            makeReferencePathTransformFromEntry(bodyLength, bodyLength, curvature));
    const Vector_t<double, 3> exitLocal = entryToExitLocal.transformTo(entryCartesian);
    if (exitLocal(2) >= 0.0) {
        return Vector_t<double, 3>(exitLocal(0), exitLocal(1), bodyLength + exitLocal(2));
    }

    const double radius = 1.0 / curvature;
    const double phi    = std::atan2(entryCartesian(2), entryCartesian(0) + radius);
    const double s      = phi / curvature;
    const double radialDistance =
            std::hypot(entryCartesian(0) + radius, entryCartesian(2)) - std::abs(radius);
    return Vector_t<double, 3>(radialDistance, entryCartesian(1), s);
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

    const double bodyLength         = getReferencePathLength();
    const double fieldBegin         = -getEntryFringeSupportLength();
    const double fieldEnd           = bodyLength + getExitFringeSupportLength();
    const double entryFringe        = getEntryFringeSupportLength();
    const double exitFringe         = getExitFringeSupportLength();
    const double referenceCurvature = getReferencePathCurvature();
    const double signedCurvature    = getSignedCurvature();
    const double dipoleFieldTesla   = getField().getNormalComponent(0);
    const double bendRigidityScale =
            (std::abs(signedCurvature) > 1.0e-15) ? dipoleFieldTesla / signedCurvature : 0.0;
    const double entryEdgeHorizontalGradient =
            (entryFringe > 0.0) ? bendRigidityScale * getEntryEdgeHorizontalStrength() / entryFringe
                                : 0.0;
    const double exitEdgeHorizontalGradient =
            (exitFringe > 0.0) ? bendRigidityScale * getExitEdgeHorizontalStrength() / exitFringe
                               : 0.0;
    const double entryEdgeVerticalGradient =
            (entryFringe > 0.0) ? -bendRigidityScale * getEntryEdgeVerticalStrength() / entryFringe
                                : 0.0;
    const double exitEdgeVerticalGradient =
            (exitFringe > 0.0) ? -bendRigidityScale * getExitEdgeVerticalStrength() / exitFringe
                               : 0.0;
    const matrix3x3_t entryToSliceRotation = slice.entryToSliceRotation;
    const Vector_t<double, 3> entryOrigin  = slice.entryOrigin;

    Kokkos::parallel_for(
            "BendBase::applySlice", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const Vector_t<double, 3> entryCartesian =
                        prod_vector_transpose(entryToSliceRotation, Rview(i)) + entryOrigin;
                Vector_t<double, 3> entryChartPosition = entryCartesian;
                if (std::abs(referenceCurvature) > 1.0e-15) {
                    const double radius = 1.0 / referenceCurvature;
                    const double phi    = std::atan2(entryCartesian(2), entryCartesian(0) + radius);
                    const double s      = phi / referenceCurvature;
                    const double radialDistance =
                            std::hypot(entryCartesian(0) + radius, entryCartesian(2))
                            - std::abs(radius);
                    entryChartPosition(0) = radialDistance;
                    entryChartPosition(1) = entryCartesian(1);
                    entryChartPosition(2) = s;
                }

                const double zEntry = entryChartPosition(2);
                if (zEntry < slice.sBegin || zEntry > slice.sEnd) {
                    return;
                }
                if (zEntry < fieldBegin || zEntry > fieldEnd) {
                    return;
                }

                Vector_t<double, 3> Bf(0.0);
                const double x     = entryChartPosition(0);
                const double y     = entryChartPosition(1);
                const double scale = evaluateFieldScale(zEntry, fieldBegin, bodyLength, fieldEnd);

                if (normal.extent(0) > 0) {
                    Bf(1) += scale * normal(0);
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

                double horizontalGradient = 0.0;
                double verticalGradient   = 0.0;
                if (zEntry <= 0.0 && entryFringe > 0.0) {
                    horizontalGradient = entryEdgeHorizontalGradient;
                    verticalGradient   = entryEdgeVerticalGradient;
                } else if (zEntry >= bodyLength && exitFringe > 0.0) {
                    horizontalGradient = exitEdgeHorizontalGradient;
                    verticalGradient   = exitEdgeVerticalGradient;
                }
                if (horizontalGradient != 0.0 || verticalGradient != 0.0) {
                    Bf(0) += verticalGradient * y;
                    Bf(1) += horizontalGradient * x;
                }

                Vector_t<double, 3> Bentry = Bf;
                if (std::abs(referenceCurvature) > 1.0e-15) {
                    const double phi  = referenceCurvature * zEntry;
                    const double cphi = std::cos(phi);
                    const double sphi = std::sin(phi);
                    Bentry(0)         = cphi * Bf(0) - sphi * Bf(2);
                    Bentry(1)         = Bf(1);
                    Bentry(2)         = sphi * Bf(0) + cphi * Bf(2);
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

    computeFieldHost(R, getField(), B);
    B *= getFieldScale(R(2));
    const double entryFringe   = getEntryFringeSupportLength();
    const double exitFringe    = getExitFringeSupportLength();
    const double rigidityScale = (std::abs(getSignedCurvature()) > 1.0e-15)
                                         ? getField().getNormalComponent(0) / getSignedCurvature()
                                         : 0.0;
    if (R(2) <= 0.0 && entryFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getEntryEdgeHorizontalStrength() / entryFringe);
        const double verticalGradient =
                -rigidityScale * (getEntryEdgeVerticalStrength() / entryFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
    } else if (R(2) >= getElementLength() && exitFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getExitEdgeHorizontalStrength() / exitFringe);
        const double verticalGradient =
                -rigidityScale * (getExitEdgeVerticalStrength() / exitFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
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

    computeFieldHost(R, getField(), B);
    B *= getFieldScale(R(2));
    const double entryFringe   = getEntryFringeSupportLength();
    const double exitFringe    = getExitFringeSupportLength();
    const double rigidityScale = (std::abs(getSignedCurvature()) > 1.0e-15)
                                         ? getField().getNormalComponent(0) / getSignedCurvature()
                                         : 0.0;
    if (R(2) <= 0.0 && entryFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getEntryEdgeHorizontalStrength() / entryFringe);
        const double verticalGradient =
                -rigidityScale * (getEntryEdgeVerticalStrength() / entryFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
    } else if (R(2) >= getElementLength() && exitFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getExitEdgeHorizontalStrength() / exitFringe);
        const double verticalGradient =
                -rigidityScale * (getExitEdgeVerticalStrength() / exitFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
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

    computeFieldHost(R, getField(), B);
    B *= getFieldScale(R(2));
    const double entryFringe   = getEntryFringeSupportLength();
    const double exitFringe    = getExitFringeSupportLength();
    const double rigidityScale = (std::abs(getSignedCurvature()) > 1.0e-15)
                                         ? getField().getNormalComponent(0) / getSignedCurvature()
                                         : 0.0;
    if (R(2) <= 0.0 && entryFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getEntryEdgeHorizontalStrength() / entryFringe);
        const double verticalGradient =
                -rigidityScale * (getEntryEdgeVerticalStrength() / entryFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
    } else if (R(2) >= getElementLength() && exitFringe > 0.0) {
        const double horizontalGradient =
                rigidityScale * (getExitEdgeHorizontalStrength() / exitFringe);
        const double verticalGradient =
                -rigidityScale * (getExitEdgeVerticalStrength() / exitFringe);
        B(0) += verticalGradient * R(1);
        B(1) += horizontalGradient * R(0);
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

    // For a positive bend angle the placed design path curves toward negative
    // entry-frame x. With OPALX's local field convention this requires
    // q v_z B_y to point toward the reference-path center of curvature, giving
    // the usual rigidity scaling B_y = p h / (q c).
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
    return evaluateFieldScale(z, fieldBegin_m, getElementLength(), fieldEnd_m);
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

void BendBase::updateFieldSupportExtent() {
    const double cosFloor         = 1.0e-6;
    const double entryDenominator = std::max(std::abs(std::cos(getEntranceAngle())), cosFloor);
    const double exitDenominator  = std::max(std::abs(std::cos(getExitAngle())), cosFloor);
    const double entryFringe      = fringeHalfGap_m * fringeIntegral_m / entryDenominator;
    const double exitFringe       = fringeHalfGap_m * fringeIntegral_m / exitDenominator;
    fieldBegin_m                  = -entryFringe;
    fieldEnd_m                    = getElementLength() + exitFringe;
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
