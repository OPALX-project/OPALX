#include "Algorithms/LinearTransferMapBuilder.h"
#include "Algorithms/CompensatedSum.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

namespace {
    constexpr double boundaryTolerance = 1.0e-12;

    ExternalFieldRayTracker::State rayState(const LinearTransferMapReference& reference) {
        return {reference.position, reference.momentum, reference.time, reference.pathLength,
                reference.positionCorrection, reference.timeCorrection, reference.pathLengthCorrection};
    }

    void assignRayState(LinearTransferMapReference& reference, const ExternalFieldRayTracker::State& ray) {
        reference.position = ray.position;
        reference.momentum = ray.momentum;
        reference.time = ray.time;
        reference.pathLength = ray.pathLength;
        reference.positionCorrection = ray.positionCorrection;
        reference.timeCorrection = ray.timeCorrection;
        reference.pathLengthCorrection = ray.pathLengthCorrection;
    }

    double flightTime(const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit) {
        return compensated::difference(exit.time, exit.timeCorrection, entrance.time, entrance.timeCorrection);
    }

    Vector_t<double, 3> displacement(const ExternalFieldRayTracker::State& ray,
                                     const LinearTransferMapReference& reference) {
        Vector_t<double, 3> result;
        for (unsigned component = 0; component < 3; ++component)
            result(component) = compensated::difference(
                    ray.position(component), ray.positionCorrection(component),
                    reference.position(component), reference.positionCorrection(component));
        return result;
    }

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
    : LinearTransferMapBuilder(beamline, reference, dt, Settings{}) {}

LinearTransferMapBuilder::LinearTransferMapBuilder(
        OpalBeamline& beamline, const PartData& reference, const double dt,
        const Settings& settings)
    : itsOpalBeamline_m(beamline),
      tracker_m(beamline, reference, settings.integrationMethod),
      dt_m(dt),
      settings_m(settings) {
    settings_m.validate();
    if (!std::isfinite(dt) || dt == 0.0)
        throw OpalException(
                "LinearTransferMapBuilder", "The ray time step must be finite and nonzero.");
}

void LinearTransferMapBuilder::Settings::validate() const {
    if (richardsonLevels > maximumRichardsonLevels)
        throw OpalException(
                "LinearTransferMapBuilder::Settings",
                "LINEARTRANSFERMAPRICHARDSON must be an integer from 0 through 4.");
    for (const double step : finiteDifferenceSteps) {
        if (!std::isfinite(step) || step <= 0.0 || std::ldexp(step, -int(richardsonLevels)) == 0.0)
            throw OpalException(
                    "LinearTransferMapBuilder::Settings",
                    "LINEARTRANSFERMAPSTEPS must contain six finite positive amplitudes "
                    "that remain nonzero after refinement.");
    }
    if (finiteDifferenceSteps[5] >= 1.0)
        throw OpalException(
                "LinearTransferMapBuilder::Settings",
                "The delta amplitude in LINEARTRANSFERMAPSTEPS must be less than 1 "
                "to keep both perturbed momenta forward and nonzero.");
    ExternalFieldRayTracker::integrationMethodName(integrationMethod);
}

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
        return element->isInsideBody(itsOpalBeamline_m.transformToLocalCS(element, ray.position));
    };

    const RayState start = rayState(before);
    double lower   = 0.0;
    double upper   = flightTime(before, after);
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
    LinearTransferMapReference boundary = transportFrame(before, trial.momentum);
    assignRayState(boundary, trial);
    return boundary;
}

std::array<double, 6> LinearTransferMapBuilder::coordinates(
        const RayState& ray, const LinearTransferMapReference& reference) {
    const Vector_t<double, 3> offset = displacement(ray, reference);
    const double longitudinalMomentum      = dot(ray.momentum, reference.sAxis);
    if (std::abs(longitudinalMomentum) < 1.0e-14) {
        throw OpalException(
                "LinearTransferMapBuilder::coordinates",
                "A shadow ray has zero longitudinal momentum.");
    }
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    return {dot(offset, reference.xAxis),
            dot(ray.momentum, reference.xAxis) / longitudinalMomentum,
            dot(offset, reference.yAxis),
            dot(ray.momentum, reference.yAxis) / longitudinalMomentum,
            -beta * Physics::c * compensated::difference(
                    ray.time, ray.timeCorrection, reference.time, reference.timeCorrection),
            euclidean_norm(ray.momentum) / referenceMomentum - 1.0};
}

LinearTransferMapBuilder::RayState LinearTransferMapBuilder::rayFromCoordinates(
        const std::array<double, 6>& coordinate, const LinearTransferMapReference& reference) {
    RayState ray = rayState(reference);
    for (unsigned component = 0; component < 3; ++component) {
        compensated::add(coordinate[0] * reference.xAxis(component), ray.position(component),
                         ray.positionCorrection(component));
        compensated::add(coordinate[2] * reference.yAxis(component), ray.position(component),
                         ray.positionCorrection(component));
    }
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double momentum          = referenceMomentum * (1.0 + coordinate[5]);
    Vector_t<double, 3> direction =
            coordinate[1] * reference.xAxis + coordinate[3] * reference.yAxis + reference.sAxis;
    ray.momentum      = momentum * normalized(direction);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    compensated::add(-coordinate[4] / (beta * Physics::c), ray.time, ray.timeCorrection);
    return ray;
}

LinearTransferMapBuilder::RayState LinearTransferMapBuilder::trackRayToExit(
        const RayState& initial, const LinearTransferMapReference& exit,
        const double referenceFlightTime) {
    RayState previous                        = initial;
    Vector_t<double, 3> previousDisplacement = displacement(previous, exit);
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
        const Vector_t<double, 3> nextDisplacement = displacement(next, exit);
        const double nextDistance                  = dot(nextDisplacement, exit.sAxis);
        const double elapsed = std::abs(compensated::difference(
                next.time, next.timeCorrection, initial.time, initial.timeCorrection));
        // A one-turn entrance and exit can be the same plane. Do not accept the launch
        // crossing (or the opposite-side crossing) as the ray's return. This is a local
        // linear-map search about the supplied reference flight, not an arbitrary orbit search.
        if (elapsed >= 0.5 * std::abs(referenceFlightTime) && crossed(previousDistance, nextDistance)) {
            double lower   = 0.0;
            double upper   = dt_m;
            RayState trial = next;
            for (int iteration = 0; iteration < 60; ++iteration) {
                const double middle                    = 0.5 * (lower + upper);
                trial                                  = tracker_m.advance(previous, middle);
                const Vector_t<double, 3> offset = displacement(trial, exit);
                const double distance                  = dot(offset, exit.sAxis);
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

matrix6x6_t LinearTransferMapBuilder::makeCenteredMap(
        const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
        const std::array<double, 6>& steps, double& condition) {
    matrix6x6_t inputDifferences(0.0);
    matrix6x6_t outputDifferences(0.0);
    for (int column = 0; column < 6; ++column) {
        std::array<double, 6> plus{};
        std::array<double, 6> minus{};
        plus[column]                    = steps[column];
        minus[column]                   = -steps[column];
        const RayState plusEntrance  = rayFromCoordinates(plus, entrance);
        const RayState minusEntrance = rayFromCoordinates(minus, entrance);
        const RayState plusExit  = trackRayToExit(plusEntrance, exit, flightTime(entrance, exit));
        const RayState minusExit = trackRayToExit(minusEntrance, exit, flightTime(entrance, exit));
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
    condition = invert6x6(inputDifferences, inverseInput);
    if (!std::isfinite(condition) || condition > 1.0e12) {
        throw OpalException(
                "LinearTransferMapBuilder::makeLinearTransferMap",
                "The finite-difference input matrix is ill-conditioned (condition number "
                        + std::to_string(condition) + ").");
    }

    return prod(outputDifferences, inverseInput);
}

LinearTransferMap LinearTransferMapBuilder::makeLinearTransferMap(
        const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
        const std::size_t pass) {
    LinearTransferMap result;
    result.finiteDifferenceSteps = settings_m.finiteDifferenceSteps;
    result.richardsonLevels      = settings_m.richardsonLevels;
    result.integrationMethod =
            ExternalFieldRayTracker::integrationMethodName(settings_m.integrationMethod);
    std::vector<matrix6x6_t> previous;
    for (unsigned level = 0; level <= settings_m.richardsonLevels; ++level) {
        std::array<double, 6> steps;
        for (unsigned column = 0; column < 6; ++column)
            steps[column] = std::ldexp(settings_m.finiteDifferenceSteps[column], -int(level));
        double condition = 0.0;
        std::vector<matrix6x6_t> current;
        current.push_back(makeCenteredMap(entrance, exit, steps, condition));
        result.inputConditionNumber = std::max(result.inputConditionNumber, condition);
        for (unsigned order = 1; order <= level; ++order) {
            matrix6x6_t refined(0.0);
            const double denominator = std::ldexp(1.0, 2 * int(order)) - 1.0;
            for (unsigned row = 0; row < 6; ++row)
                for (unsigned column = 0; column < 6; ++column) {
                    const double fine = current.back()(row, column);
                    refined(row, column) =
                            fine + (fine - previous[order - 1](row, column)) / denominator;
                }
            current.push_back(refined);
        }
        if (level > 0) {
            std::array<double, 6> change{};
            for (unsigned row = 0; row < 6; ++row)
                for (unsigned column = 0; column < 6; ++column)
                    change[column] = std::max(
                            change[column],
                            std::abs(current.back()(row, column) - previous.back()(row, column)));
            result.richardsonCorrection = change;
        }
        result.finestFiniteDifferenceSteps = steps;
        previous                           = std::move(current);
    }
    result.matrix                    = previous.back();
    result.entrance                  = entrance;
    result.exit                      = exit;
    result.pass                      = pass;
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
                return compensated::difference(sample.state.pathLength,
                        sample.state.pathLengthCorrection, pathLength, 0.0) < 0.0;
            });
    if (firstAtOrAfterStart == referenceSamples_m.end()) return result;
    if (firstAtOrAfterStart != referenceSamples_m.begin()
        && compensated::difference(firstAtOrAfterStart->state.pathLength,
                firstAtOrAfterStart->state.pathLengthCorrection, transferMapStartPathLength_m, 0.0) > 0.0) {
        const auto& before    = std::prev(firstAtOrAfterStart)->state;
        const auto& after     = firstAtOrAfterStart->state;
        const RayState clipped = tracker_m.advanceToPathLength(
                rayState(before), flightTime(before, after), transferMapStartPathLength_m);
        LinearTransferMapReference start = transportFrame(before, clipped.momentum);
        assignRayState(start, clipped);
        start.pathLength                 = transferMapStartPathLength_m;
        start.pathLengthCorrection       = 0.0;

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
        // "active" here means nominal body owners, not active field supports.
        IndexMap::value_t active;
    };

    std::vector<BoundaryEvent> events;
    auto previousSet = itsOpalBeamline_m.getBodyElements(referenceSamples_m.front().state.position);
    validateActiveSet(itsOpalBeamline_m.getElements(referenceSamples_m.front().state.position));

    for (std::size_t sample = 1; sample < referenceSamples_m.size(); ++sample) {
        const auto& before = referenceSamples_m[sample - 1].state;
        const auto& after  = referenceSamples_m[sample].state;
        auto currentSet    = itsOpalBeamline_m.getBodyElements(after.position);
        validateActiveSet(itsOpalBeamline_m.getElements(after.position));

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
            itsOpalBeamline_m.getBodyElements(referenceSamples_m.front().state.position);
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
    std::size_t sampleIndex = 1;
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        const auto& interval = segments[segment];
        validateActiveSet(interval.active);
        LinearTransferMap map = makeLinearTransferMap(interval.entrance, interval.exit, segment);
        map.segment           = segment;
        // Ownership follows bodies, while contributors can change inside one body interval.
        // Contributor names are sampled at reference-bracket chord midpoints for diagnostics.
        // They are not a field mask: each shadow ray queries its own summed field at every
        // stage, including support-boundary subdivisions. A bracket may cross a body boundary.
        IndexMap::value_t contributors;
        map.includesOverlappingFields = interval.active.size() > 1;
        for (; sampleIndex < referenceSamples_m.size(); ++sampleIndex) {
            const auto& before = referenceSamples_m[sampleIndex - 1].state;
            const auto& after = referenceSamples_m[sampleIndex].state;
            if (before.pathLength >= interval.exit.pathLength) break;
            if (after.pathLength <= interval.entrance.pathLength) continue;
            auto fields = itsOpalBeamline_m.getElements(0.5 * (before.position + after.position));
            validateActiveSet(fields);
            map.includesOverlappingFields = map.includesOverlappingFields || fields.size() > 1;
            contributors.insert(fields.begin(), fields.end());
        }
        if (sampleIndex > 1) --sampleIndex; // a sample bracket may straddle a nominal boundary
        for (const auto& element : contributors) {
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
