// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_LINEAR_TRANSFER_MAP_BUILDER_H
#define OPAL_LINEAR_TRANSFER_MAP_BUILDER_H

#include <optional>
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/IndexMap.h"
#include "Structure/LinearTransferMap.h"

class OpalBeamline;
class PartData;

/**
 * @brief Construct external-field linear maps from a sampled reference orbit.
 *
 * Owns Bishop-frame transport, boundary refinement, the twelve perturbed rays,
 * the pivoted finite-difference solve and ordered composition. The caller owns
 * the reference pass and attaches the returned maps to their runtime elements.
 * The builder neither modifies the beamline nor writes diagnostics to stdout.
 * See LinearTransferMap for coordinates, units and the finite-difference equations.
 */
class LinearTransferMapBuilder {
public:
    struct ReferenceSample {
        LinearTransferMapReference state;
    };
    struct Segment {
        LinearTransferMap map;
        IndexMap::value_t owners;
    };
    struct Result {
        /// Includes field-free segments (owners empty); each overlap occurs only once.
        std::vector<Segment> segments;
        std::optional<matrix6x6_t> combined;
        std::optional<double> determinantResidual;
        std::optional<double> symplecticResidual;
    };

    LinearTransferMapBuilder(OpalBeamline& beamline, const PartData& reference, double dt);
    /// Samples are copied so pre-roll clipping cannot mutate the caller's reference orbit.
    Result build(std::vector<ReferenceSample> samples, double startPathLength);

    /// Initialize the transverse axes from the beamline frame, without a roll singularity.
    static LinearTransferMapReference initialFrame(
            OpalBeamline& beamline, const Vector_t<double, 3>& momentum);
    /// Minimum-rotation transport of an existing orthonormal frame to a new tangent.
    static LinearTransferMapReference transportFrame(
            const LinearTransferMapReference& frame, const Vector_t<double, 3>& momentum);

private:
    using RayState = ExternalFieldRayTracker::State;
    OpalBeamline& itsOpalBeamline_m;
    ExternalFieldRayTracker tracker_m;
    double dt_m;

    LinearTransferMapReference refineBoundary(
            const std::shared_ptr<ElementBase>& element, const LinearTransferMapReference& before,
            const LinearTransferMapReference& after, bool entering);
    RayState trackRayToExit(
            const RayState& initial, const LinearTransferMapReference& exit,
            double referenceFlightTime);
    LinearTransferMap makeLinearTransferMap(
            const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
            std::size_t pass);
    static std::array<double, 6> coordinates(
            const RayState& ray, const LinearTransferMapReference& reference);
    static RayState rayFromCoordinates(
            const std::array<double, 6>& coordinates, const LinearTransferMapReference& reference);
};
#endif
