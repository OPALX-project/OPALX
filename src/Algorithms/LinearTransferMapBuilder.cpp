#include "Algorithms/LinearTransferMapBuilder.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

namespace {
    constexpr std::array<double, 6> transferMapSteps{1.0e-3, 1.0e-3, 1.0e-3,
                                                     1.0e-3, 1.0e-3, 1.0e-3};
    constexpr double boundaryTolerance = 1.0e-12;

    Vector_t<double, 3> crossProduct(const Vector_t<double, 3>& a, const Vector_t<double, 3>& b) {
        return Vector_t<double, 3>(
                a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0));
    }

    Vector_t<double, 3> normalized(const Vector_t<double, 3>& value) {
        const double norm = euclidean_norm(value);
        if (!(norm > 0.0)) {
            throw OpalException(
                    "LinearTransferMapBuilder::normalized", "Cannot normalize a zero vector.");
        }
        return value / norm;
    }

    Vector_t<double, 3> rotateRodrigues(
            const Vector_t<double, 3>& value, const Vector_t<double, 3>& axis, const double angle) {
        return value * std::cos(angle) + crossProduct(axis, value) * std::sin(angle)
               + axis * dot(axis, value) * (1.0 - std::cos(angle));
    }

    double invert6x6(const matrix6x6_t& input, matrix6x6_t& inverse) {
        double augmented[6][12]{};
        double normInput = 0.0;
        for (int row = 0; row < 6; ++row) {
            double rowSum = 0.0;
            for (int col = 0; col < 6; ++col) {
                augmented[row][col]     = input(row, col);
                augmented[row][col + 6] = row == col ? 1.0 : 0.0;
                rowSum += std::abs(input(row, col));
            }
            normInput = std::max(normInput, rowSum);
        }

        for (int col = 0; col < 6; ++col) {
            int pivot = col;
            for (int row = col + 1; row < 6; ++row) {
                if (std::abs(augmented[row][col]) > std::abs(augmented[pivot][col])) pivot = row;
            }
            if (std::abs(augmented[pivot][col]) < 1.0e-18) {
                throw OpalException(
                        "LinearTransferMapBuilder::makeLinearTransferMap",
                        "The finite-difference input matrix is singular.");
            }
            if (pivot != col) {
                for (int j = 0; j < 12; ++j)
                    std::swap(augmented[pivot][j], augmented[col][j]);
            }
            const double diagonal = augmented[col][col];
            for (int j = 0; j < 12; ++j)
                augmented[col][j] /= diagonal;
            for (int row = 0; row < 6; ++row) {
                if (row == col) continue;
                const double factor = augmented[row][col];
                for (int j = 0; j < 12; ++j)
                    augmented[row][j] -= factor * augmented[col][j];
            }
        }

        double normInverse = 0.0;
        for (int row = 0; row < 6; ++row) {
            double rowSum = 0.0;
            for (int col = 0; col < 6; ++col) {
                inverse(row, col) = augmented[row][col + 6];
                rowSum += std::abs(inverse(row, col));
            }
            normInverse = std::max(normInverse, rowSum);
        }
        return normInput * normInverse;
    }

    double symplecticResidual(const matrix6x6_t& map) {
        matrix6x6_t j(0.0);
        for (int pair = 0; pair < 3; ++pair) {
            j(2 * pair, 2 * pair + 1) = 1.0;
            j(2 * pair + 1, 2 * pair) = -1.0;
        }
        const matrix6x6_t residual = prod(get_transpose(map), prod(j, map));
        double maximum             = 0.0;
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                maximum = std::max(maximum, std::abs(residual(row, col) - j(row, col)));
            }
        }
        return maximum;
    }

    double determinant(const matrix6x6_t& map) {
        double work[6][6]{};
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                work[row][col] = map(row, col);
            }
        }

        double result = 1.0;
        for (int col = 0; col < 6; ++col) {
            int pivot = col;
            for (int row = col + 1; row < 6; ++row) {
                if (std::abs(work[row][col]) > std::abs(work[pivot][col])) pivot = row;
            }
            if (work[pivot][col] == 0.0) return 0.0;
            if (pivot != col) {
                for (int j = col; j < 6; ++j) {
                    std::swap(work[pivot][j], work[col][j]);
                }
                result = -result;
            }
            const double diagonal = work[col][col];
            result *= diagonal;
            for (int row = col + 1; row < 6; ++row) {
                const double factor = work[row][col] / diagonal;
                for (int j = col + 1; j < 6; ++j) {
                    work[row][j] -= factor * work[col][j];
                }
            }
        }
        return result;
    }

    double determinantResidual(const matrix6x6_t& map) { return std::abs(determinant(map) - 1.0); }
}  // namespace

LinearTransferMapBuilder::LinearTransferMapBuilder(
        OpalBeamline& beamline, const PartData& reference, const double dt)
    : itsOpalBeamline_m(beamline), tracker_m(beamline, reference), dt_m(dt) {}

LinearTransferMapReference LinearTransferMapBuilder::initialFrame(
        OpalBeamline& beamline, const Vector_t<double, 3>& momentum) {
    LinearTransferMapReference state;
    state.sAxis = normalized(momentum);
    Vector_t<double, 3> xAxis =
            beamline.getCSTrafoLab2Local().rotateFrom(Vector_t<double, 3>(1.0, 0.0, 0.0));
    xAxis -= state.sAxis * dot(xAxis, state.sAxis);
    if (euclidean_norm(xAxis) < 1.0e-12) {
        xAxis = beamline.getCSTrafoLab2Local().rotateFrom(Vector_t<double, 3>(0.0, 1.0, 0.0));
        xAxis -= state.sAxis * dot(xAxis, state.sAxis);
    }
    state.xAxis = normalized(xAxis);
    state.yAxis = normalized(crossProduct(state.sAxis, state.xAxis));

    return state;
}

LinearTransferMapReference LinearTransferMapBuilder::transportFrame(
        const LinearTransferMapReference& frame, const Vector_t<double, 3>& momentum) {
    LinearTransferMapReference result  = frame;
    const Vector_t<double, 3> newS     = normalized(momentum);
    const Vector_t<double, 3> rotation = crossProduct(frame.sAxis, newS);
    const double sine                  = euclidean_norm(rotation);
    const double cosine                = std::max(-1.0, std::min(1.0, dot(frame.sAxis, newS)));

    Vector_t<double, 3> newX = frame.xAxis;
    if (sine > 1.0e-14) {
        const Vector_t<double, 3> axis = rotation / sine;
        newX = rotateRodrigues(frame.xAxis, axis, std::atan2(sine, cosine));
    } else if (cosine < 0.0) {
        // A 180-degree reversal has no unique minimum-rotation axis.  Preserve the old
        // transverse x direction, which is already perpendicular to both tangents.
        newX = frame.xAxis;
    }

    newX         = normalized(newX - newS * dot(newX, newS));
    result.sAxis = newS;
    result.xAxis = newX;
    result.yAxis = normalized(crossProduct(newS, newX));
    return result;
}

LinearTransferMapReference LinearTransferMapBuilder::refineBoundary(
        const std::shared_ptr<ElementBase>& element, const LinearTransferMapReference& before,
        const LinearTransferMapReference& after, const bool entering) {
    const auto isInside = [&](const RayState& ray) {
        return element->isInside(itsOpalBeamline_m.transformToLocalCS(element, ray.position));
    };

    const RayState start{before.position, before.momentum, before.time};
    double lower   = 0.0;
    double upper   = after.time - before.time;
    RayState trial = start;
    for (int iteration = 0; iteration < 60; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        trial               = tracker_m.advance(start, middle);
        if (isInside(trial) == entering) {
            upper = middle;
        } else {
            lower = middle;
        }
        if (std::abs(upper - lower) <= boundaryTolerance * std::abs(dt_m)) break;
    }
    const double fraction               = (after.time != before.time)
                                                  ? (trial.time - before.time) / (after.time - before.time)
                                                  : 0.0;
    LinearTransferMapReference boundary = transportFrame(before, trial.momentum);
    boundary.position                   = trial.position;
    boundary.momentum                   = trial.momentum;
    boundary.time                       = trial.time;
    boundary.pathLength = before.pathLength + fraction * (after.pathLength - before.pathLength);
    return boundary;
}

std::array<double, 6> LinearTransferMapBuilder::coordinates(
        const RayState& ray, const LinearTransferMapReference& reference) {
    const Vector_t<double, 3> displacement = ray.position - reference.position;
    const double longitudinalMomentum      = dot(ray.momentum, reference.sAxis);
    if (std::abs(longitudinalMomentum) < 1.0e-14) {
        throw OpalException(
                "LinearTransferMapBuilder::coordinates",
                "A shadow ray has zero longitudinal momentum.");
    }
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    return {dot(displacement, reference.xAxis),
            dot(ray.momentum, reference.xAxis) / longitudinalMomentum,
            dot(displacement, reference.yAxis),
            dot(ray.momentum, reference.yAxis) / longitudinalMomentum,
            -beta * Physics::c * (ray.time - reference.time),
            euclidean_norm(ray.momentum) / referenceMomentum - 1.0};
}

LinearTransferMapBuilder::RayState LinearTransferMapBuilder::rayFromCoordinates(
        const std::array<double, 6>& coordinate, const LinearTransferMapReference& reference) {
    RayState ray;
    ray.position =
            reference.position + coordinate[0] * reference.xAxis + coordinate[2] * reference.yAxis;
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double momentum          = referenceMomentum * (1.0 + coordinate[5]);
    Vector_t<double, 3> direction =
            coordinate[1] * reference.xAxis + coordinate[3] * reference.yAxis + reference.sAxis;
    ray.momentum      = momentum * normalized(direction);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    ray.time          = reference.time - coordinate[4] / (beta * Physics::c);
    return ray;
}

LinearTransferMapBuilder::RayState LinearTransferMapBuilder::trackRayToExit(
        const RayState& initial, const LinearTransferMapReference& exit,
        const double referenceFlightTime) {
    RayState previous                        = initial;
    Vector_t<double, 3> previousDisplacement = previous.position - exit.position;
    double previousDistance                  = dot(previousDisplacement, exit.sAxis);
    const double expectedDirection           = dot(exit.momentum, exit.sAxis) >= 0.0 ? 1.0 : -1.0;
    const auto crossed = [&](const double oldDistance, const double newDistance) {
        return expectedDirection * oldDistance <= 0.0 && expectedDirection * newDistance >= 0.0;
    };

    const std::size_t referenceSteps =
            static_cast<std::size_t>(std::ceil(std::abs(referenceFlightTime / dt_m)));
    const std::size_t maximumSteps = std::max<std::size_t>(100, 4 * referenceSteps + 100);
    for (std::size_t step = 0; step < maximumSteps; ++step) {
        RayState next                              = tracker_m.advance(previous, dt_m);
        const Vector_t<double, 3> nextDisplacement = next.position - exit.position;
        const double nextDistance                  = dot(nextDisplacement, exit.sAxis);
        if (crossed(previousDistance, nextDistance)) {
            double lower   = 0.0;
            double upper   = dt_m;
            RayState trial = next;
            for (int iteration = 0; iteration < 60; ++iteration) {
                const double middle                    = 0.5 * (lower + upper);
                trial                                  = tracker_m.advance(previous, middle);
                const Vector_t<double, 3> displacement = trial.position - exit.position;
                const double distance                  = dot(displacement, exit.sAxis);
                if (crossed(previousDistance, distance)) {
                    upper = middle;
                } else {
                    lower = middle;
                }
                if (std::abs(upper - lower) <= boundaryTolerance * std::abs(dt_m)) break;
            }
            return trial;
        }
        previous         = next;
        previousDistance = nextDistance;
    }
    throw OpalException(
            "LinearTransferMapBuilder::trackRayToExit",
            "A transfer-map shadow ray did not reach the recorded exit plane.");
}

LinearTransferMap LinearTransferMapBuilder::makeLinearTransferMap(
        const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
        const std::size_t pass) {
    matrix6x6_t inputDifferences(0.0);
    matrix6x6_t outputDifferences(0.0);
    for (int column = 0; column < 6; ++column) {
        std::array<double, 6> plus{};
        std::array<double, 6> minus{};
        plus[column]                 = transferMapSteps[column];
        minus[column]                = -transferMapSteps[column];
        const RayState plusEntrance  = rayFromCoordinates(plus, entrance);
        const RayState minusEntrance = rayFromCoordinates(minus, entrance);
        const RayState plusExit  = trackRayToExit(plusEntrance, exit, exit.time - entrance.time);
        const RayState minusExit = trackRayToExit(minusEntrance, exit, exit.time - entrance.time);
        const auto encodedPlusEntrance  = coordinates(plusEntrance, entrance);
        const auto encodedMinusEntrance = coordinates(minusEntrance, entrance);
        const auto encodedPlusExit      = coordinates(plusExit, exit);
        const auto encodedMinusExit     = coordinates(minusExit, exit);
        for (int row = 0; row < 6; ++row) {
            inputDifferences(row, column) =
                    0.5 * (encodedPlusEntrance[row] - encodedMinusEntrance[row]);
            outputDifferences(row, column) = 0.5 * (encodedPlusExit[row] - encodedMinusExit[row]);
        }
    }

    matrix6x6_t inverseInput(0.0);
    const double condition = invert6x6(inputDifferences, inverseInput);
    if (!std::isfinite(condition) || condition > 1.0e12) {
        throw OpalException(
                "LinearTransferMapBuilder::makeLinearTransferMap",
                "The finite-difference input matrix is ill-conditioned (condition number "
                        + std::to_string(condition) + ").");
    }

    LinearTransferMap result;
    result.matrix                    = prod(outputDifferences, inverseInput);
    result.finiteDifferenceSteps     = transferMapSteps;
    result.entrance                  = entrance;
    result.exit                      = exit;
    result.pass                      = pass;
    result.inputConditionNumber      = condition;
    result.determinantResidual       = determinantResidual(result.matrix);
    result.symplecticResidual        = symplecticResidual(result.matrix);
    result.includesOverlappingFields = false;
    return result;
}

LinearTransferMapBuilder::Result LinearTransferMapBuilder::build(
        std::vector<ReferenceSample> referenceSamples_m, double transferMapStartPathLength_m) {
    Result result;
    if (referenceSamples_m.size() < 2) return result;

    // OrbitThreader first tracks backwards by the bunch extent.  That path is needed to build
    // the IndexMap, but a LINE/one-turn transfer map must start at the requested reference s.
    const auto firstAtOrAfterStart = std::lower_bound(
            referenceSamples_m.begin(), referenceSamples_m.end(), transferMapStartPathLength_m,
            [](const ReferenceSample& sample, const double pathLength) {
                return sample.state.pathLength < pathLength;
            });
    if (firstAtOrAfterStart == referenceSamples_m.end()) return result;
    if (firstAtOrAfterStart != referenceSamples_m.begin()
        && firstAtOrAfterStart->state.pathLength > transferMapStartPathLength_m) {
        const auto& before    = std::prev(firstAtOrAfterStart)->state;
        const auto& after     = firstAtOrAfterStart->state;
        const double fraction = (transferMapStartPathLength_m - before.pathLength)
                                / (after.pathLength - before.pathLength);
        const RayState clipped = tracker_m.advance(
                {before.position, before.momentum, before.time},
                fraction * (after.time - before.time));
        LinearTransferMapReference start = transportFrame(before, clipped.momentum);
        start.position                   = clipped.position;
        start.momentum                   = clipped.momentum;
        start.time                       = clipped.time;
        start.pathLength                 = transferMapStartPathLength_m;

        std::vector<ReferenceSample> clippedSamples;
        clippedSamples.reserve(
                1 + static_cast<std::size_t>(referenceSamples_m.end() - firstAtOrAfterStart));
        clippedSamples.push_back({start});
        clippedSamples.insert(clippedSamples.end(), firstAtOrAfterStart, referenceSamples_m.end());
        referenceSamples_m = std::move(clippedSamples);
    } else {
        referenceSamples_m.erase(referenceSamples_m.begin(), firstAtOrAfterStart);
    }
    if (referenceSamples_m.size() < 2) return result;

    const auto validateActiveSet = [](const IndexMap::value_t& active) {
        for (const auto& element : active) {
            if (element->getType() == ElementType::RFCAVITY
                || element->getType() == ElementType::TRAVELINGWAVE) {
                throw OpalException(
                        "LinearTransferMapBuilder::calculateLinearTransferMaps",
                        "RF elements are not supported by the first linear-transfer-map "
                        "implementation (element '"
                                + element->getName() + "').");
            }
        }
    };

    using ElementPtr = std::shared_ptr<ElementBase>;
    struct BoundaryEvent {
        LinearTransferMapReference state;
        ElementPtr element;
        bool entering;
    };
    struct Segment {
        LinearTransferMapReference entrance;
        LinearTransferMapReference exit;
        IndexMap::value_t active;
    };

    std::vector<BoundaryEvent> events;
    auto previousSet = itsOpalBeamline_m.getElements(referenceSamples_m.front().state.position);
    validateActiveSet(previousSet);

    for (std::size_t sample = 1; sample < referenceSamples_m.size(); ++sample) {
        const auto& before = referenceSamples_m[sample - 1].state;
        const auto& after  = referenceSamples_m[sample].state;
        auto currentSet    = itsOpalBeamline_m.getElements(after.position);
        validateActiveSet(currentSet);

        IndexMap::value_t exited;
        std::set_difference(
                previousSet.begin(), previousSet.end(), currentSet.begin(), currentSet.end(),
                std::inserter(exited, exited.begin()));
        for (const auto& element : exited) {
            events.push_back({refineBoundary(element, before, after, false), element, false});
        }

        IndexMap::value_t entered;
        std::set_difference(
                currentSet.begin(), currentSet.end(), previousSet.begin(), previousSet.end(),
                std::inserter(entered, entered.begin()));
        for (const auto& element : entered) {
            events.push_back({refineBoundary(element, before, after, true), element, true});
        }
        previousSet = std::move(currentSet);
    }

    std::sort(events.begin(), events.end(), [](const auto& left, const auto& right) {
        if (left.state.pathLength != right.state.pathLength) {
            return left.state.pathLength < right.state.pathLength;
        }
        // At an exactly shared boundary, close the old set before opening the new one.
        if (left.entering != right.entering) return left.entering < right.entering;
        return left.element->getName() < right.element->getName();
    });

    std::vector<Segment> segments;
    IndexMap::value_t active =
            itsOpalBeamline_m.getElements(referenceSamples_m.front().state.position);
    LinearTransferMapReference segmentEntrance = referenceSamples_m.front().state;
    constexpr double minimumSegmentLength      = 1.0e-12;
    for (const auto& event : events) {
        if (event.state.pathLength - segmentEntrance.pathLength > minimumSegmentLength) {
            segments.push_back({segmentEntrance, event.state, active});
        }
        if (event.entering) {
            active.insert(event.element);
        } else {
            active.erase(event.element);
        }
        segmentEntrance = event.state;
    }
    const auto& finalState = referenceSamples_m.back().state;
    if (finalState.pathLength - segmentEntrance.pathLength > minimumSegmentLength) {
        segments.push_back({segmentEntrance, finalState, active});
    }

    result.segments.reserve(segments.size());
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        const auto& interval = segments[segment];
        validateActiveSet(interval.active);
        LinearTransferMap map = makeLinearTransferMap(interval.entrance, interval.exit, segment);
        map.segment           = segment;
        map.includesOverlappingFields = interval.active.size() > 1;
        for (const auto& element : interval.active) {
            map.activeElements.push_back(element->getName());
        }
        std::sort(map.activeElements.begin(), map.activeElements.end());
        result.segments.push_back({std::move(map), interval.active});
    }

    if (result.segments.empty()) return result;
    matrix6x6_t combined;
    for (const auto& segment : result.segments) {
        combined = prod(segment.map.matrix, combined);
    }
    result.combined            = combined;
    result.determinantResidual = determinantResidual(combined);
    result.symplecticResidual  = symplecticResidual(combined);
    return result;
}
