// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_TPS_CLOSED_ORBIT_MAP_H
#define OPAL_TPS_CLOSED_ORBIT_MAP_H

#include "Algorithms/ClosedOrbitSolver.h"
#include "tps/Map.hpp"

namespace opalx::maps {
    /**
     * @brief Snapshot a converged COF linearization without any additional ray tracking.
     *
     * The local variables are deviations from the solved (x,px,y,py), with positions
     * in metres and mechanical momenta in p/(mc), not slopes. Input and output use
     * the same fixed section and expansion reference. Unit numerical chart scales
     * give F(xi)=residual+matrix*xi. The measured residual is retained exactly, not
     * replaced by zero or absorbed into a shifted output reference.
     *
     * referenceMomentum is the positive fixed total p/(mc). The solved transverse
     * momenta must admit positive longitudinal momentum on that energy shell.
     * Frame rows map lab vectors into section coordinates. Reference time is unset:
     * this autonomous 4D Poincare map has independent perturbed-ray return times.
     * T0, return duration, mass/charge and solver diagnostics remain with COF.
     *
     * Requires an active Kokkos runtime and serialized host TPS algebra. This is a
     * first-order affine approximation, not nonlinear map construction or a new
     * closure/convergence proof. No solver, field, ray or report state is changed.
     * @throws OpalException Nonconverged status, nonfinite data, invalid section or
     * invalid fixed-momentum launch.
     */
    tps::Map<double> toTpsClosedOrbitMap(
            const ClosedOrbitSolver::Result&, const CoordinateSystemTrafo& section,
            double referenceMomentum);
}  // namespace opalx::maps
#endif
