// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/TpsClosedOrbitMap.h"
#include <cmath>
#include <stdexcept>
#include <utility>
#include "Utilities/OpalException.h"

namespace opalx::maps {
    tps::Map<double> toTpsClosedOrbitMap(
            const ClosedOrbitSolver::Result& orbit, const CoordinateSystemTrafo& section,
            double referenceMomentum) {
        const auto require = [](bool valid, const char* message) {
            if (!valid) throw OpalException("toTpsClosedOrbitMap", message);
        };
        require(orbit.status == ClosedOrbitSolver::Status::Converged,
                "A converged closed-orbit result is required.");
        require(std::isfinite(referenceMomentum) && referenceMomentum > 0
                        && std::isfinite(referenceMomentum * referenceMomentum),
                "Invalid fixed total momentum.");
        for (unsigned i = 0; i < 4; ++i) {
            require(std::isfinite(orbit.coordinates[i]) && std::isfinite(orbit.residual[i]),
                    "Nonfinite closed-orbit reference/residual.");
            for (unsigned j = 0; j < 4; ++j)
                require(std::isfinite(orbit.matrix[i][j]), "Nonfinite closed-orbit matrix.");
        }
        const double longitudinalSquared = referenceMomentum * referenceMomentum
                                           - orbit.coordinates[1] * orbit.coordinates[1]
                                           - orbit.coordinates[3] * orbit.coordinates[3];
        require(std::isfinite(longitudinalSquared) && longitudinalSquared > 0,
                "The reference has no real forward longitudinal momentum.");

        tps::PhaseSpaceChart chart;
        chart.name        = "COF mechanical transverse perturbations";
        chart.coordinates = {
                {"x", "m", 1}, {"px", "p/(mc)", 1}, {"y", "m", 1}, {"py", "p/(mc)", 1}};
        chart.reference.assign(orbit.coordinates.begin(), orbit.coordinates.end());
        chart.convention        = tps::CoordinateConvention::mechanical;
        chart.referenceMomentum = referenceMomentum;
        chart.momentumUnit      = "p/(mc)";
        chart.frame.name        = "laboratory to fixed COF section";
        const auto origin       = section.getOrigin();
        const auto rotation     = section.getRotationMatrix();
        for (unsigned i = 0; i < 3; ++i) {
            chart.frame.origin[i] = origin(i);
            for (unsigned j = 0; j < 3; ++j)
                chart.frame.axes[3 * i + j] = rotation(i, j);
        }
        try {
            tps::MapContext context(chart, chart, 1);
            tps::VpsMap<double> polynomial(4, 4);
            for (unsigned i = 0; i < 4; ++i) {
                polynomial[i]    = tps::Tps<double>(1, 4);
                polynomial[i][0] = orbit.residual[i];
                for (unsigned j = 0; j < 4; ++j)
                    polynomial[i][j + 1] = orbit.matrix[i][j];
            }
            return tps::Map<double>(std::move(context), std::move(polynomial));
        } catch (const std::invalid_argument& error) {
            throw OpalException("toTpsClosedOrbitMap", error.what());
        }
    }
}  // namespace opalx::maps
