// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_COF_CMD_H
#define OPAL_COF_CMD_H
#include <optional>
#include "AbstractObjects/Action.h"
#include "Algorithms/ClosedOrbitInitialState.h"

/** Named, single-statement static-magnetic 4D closed-orbit solve.
 * A successful execution stores an independent launch snapshot on this object.
 * OUTPUT optionally names one JSON file; no orbit or text sidecars are produced.
 * Numerical stability diagnostics do not certify map convergence or confinement.
 */
class CofCmd : public Action {
public:
    CofCmd();
    CofCmd* clone(const std::string& name) override;
    void execute() override;
    /// Resolve a completed named solve; templates and unrelated objects are rejected.
    static const ClosedOrbitInitialState& findResult(const std::string& name);

private:
    CofCmd(const std::string& name, CofCmd* parent);
    std::optional<ClosedOrbitInitialState> result_m;
};
#endif
