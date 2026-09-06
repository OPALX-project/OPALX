// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_COF_CMD_H
#define OPAL_COF_CMD_H
#include "AbstractObjects/Action.h"

/** Scoped COF / RUN / ENDCOF input block for static magnetic 4D ring optics.
 * Uses the shared ray tracker, directed fixed section, Newton solver and LAPACK.
 * Does not modify BEAM energy or seed a subsequent TRACK. Numerical stability
 * diagnostics do not certify convergence of the differentiated tracking map.
 * One empty setup bunch preserves existing element reference-charge semantics;
 * no particles are emitted and no collective solve or production tracker runs.
 * Conventional aperture checks use the nearest finite design centreline,
 * independently of aperture size (all roundoff-level ties are checked). This
 * local-orbit assignment prevents opposite-arm slab losses; intersecting pipe
 * solids are not modelled. Cyclotron sectors retain native domain predicates.
 */
class CofCmd : public Action {
public:
    CofCmd();
    CofCmd* clone(const std::string& name) override;
    void execute() override;

private:
    CofCmd(const std::string& name, CofCmd* parent);
};
#endif
