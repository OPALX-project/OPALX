// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/TpsLinearTransferMap.h"
#include <cmath>
#include <stdexcept>
#include <utility>
#include "Utilities/OpalException.h"

namespace opalx::maps {
    namespace {
        tps::PhaseSpaceChart boundaryChart(const LinearTransferMapReference& reference) {
            for (unsigned i = 0; i < 3; ++i) {
                if (!std::isfinite(reference.position(i))
                    || !std::isfinite(reference.positionCorrection(i))
                    || !std::isfinite(reference.momentum(i)))
                    throw std::invalid_argument("Nonfinite boundary position/momentum.");
            }
            if (!std::isfinite(reference.time) || !std::isfinite(reference.timeCorrection)
                || !std::isfinite(reference.pathLength)
                || !std::isfinite(reference.pathLengthCorrection))
                throw std::invalid_argument("Nonfinite boundary clock/path.");
            const double momentum =
                    std::hypot(reference.momentum(0), reference.momentum(1), reference.momentum(2));
            const double longitudinal = dot(reference.momentum, reference.sAxis);
            if (!std::isfinite(momentum) || momentum <= 0 || !std::isfinite(longitudinal)
                || longitudinal == 0)
                throw std::invalid_argument("Invalid boundary reference momentum.");

            tps::PhaseSpaceChart chart;
            chart.name              = "mechanical slope perturbations";
            chart.coordinates       = {{"x", "m", 1},      {"xprime", "1", 1}, {"y", "m", 1},
                                       {"yprime", "1", 1}, {"zeta", "m", 1},   {"delta", "1", 1}};
            chart.reference         = {0, dot(reference.momentum, reference.xAxis) / longitudinal,
                                       0, dot(reference.momentum, reference.yAxis) / longitudinal,
                                       0, 0};
            chart.convention        = tps::CoordinateConvention::mechanical;
            chart.referenceTime     = reference.time - reference.timeCorrection;
            chart.referenceMomentum = momentum;
            chart.momentumUnit      = "p/(mc)";
            chart.frame.name        = "laboratory to reference section";
            for (unsigned i = 0; i < 3; ++i) {
                chart.frame.origin[i]   = reference.position(i) - reference.positionCorrection(i);
                chart.frame.axes[i]     = reference.xAxis(i);
                chart.frame.axes[3 + i] = reference.yAxis(i);
                chart.frame.axes[6 + i] = reference.sAxis(i);
            }
            return chart;
        }
    }  // namespace

    tps::Map<double> toTpsLinearTransferMap(
            const LinearTransferMapDraft& segment, const std::string& occurrenceId) {
        // Validate all matrix entries before allocating any coefficient storage.
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                if (!std::isfinite(segment.matrix(i, j)))
                    throw OpalException("toTpsLinearTransferMap", "Nonfinite segment matrix.");
        try {
            tps::MapContext context(
                    boundaryChart(segment.entrance), boundaryChart(segment.exit), 1,
                    {occurrenceId});
            tps::VpsMap<double> polynomial(6, 6);
            for (unsigned i = 0; i < 6; ++i) {
                polynomial[i] = tps::Tps<double>(1, 6);
                for (unsigned j = 0; j < 6; ++j)
                    polynomial[i][j + 1] = segment.matrix(i, j);
            }
            return tps::Map<double>(std::move(context), std::move(polynomial));
        } catch (const std::invalid_argument& error) {
            throw OpalException("toTpsLinearTransferMap", error.what());
        }
    }
}  // namespace opalx::maps
