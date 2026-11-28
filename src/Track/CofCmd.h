// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_COF_CMD_H
#define OPAL_COF_CMD_H
#include <optional>
#include "AbstractObjects/Action.h"
#include "Algorithms/ClosedOrbitOptics.h"

/** Named, single-statement static-magnetic 4D closed-orbit solve.
 * A successful execution stores an independent launch snapshot on this object.
 * OUTPUT optionally names one JSON file; no orbit or text sidecars are produced.
 * The version-1 opalx.closed_orbits export contains species metadata and one
 * beam-identified orbit with section, launch, solver, optics and return diagnostics.
 * Numerical stability diagnostics do not certify map convergence or confinement.
 */
class CofCmd : public Action {
public:
    CofCmd();
    CofCmd* clone(const std::string& name) override;
    void execute() override;
    /// Resolve a completed named solve; templates and unrelated objects are rejected.
    static const ClosedOrbitInitialState& findResult(const std::string& name);

    /// Typed optics for a particular successful execution; no implicit ray solve.
    static const ClosedOrbitOptics& findOptics(const std::string& name);

private:
    CofCmd(const std::string& name, CofCmd* parent);
    std::optional<ClosedOrbitOptics> result_m;
};
#endif
