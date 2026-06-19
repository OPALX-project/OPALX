#ifndef OPAL_SOLVER_CAPABILITIES_H
#define OPAL_SOLVER_CAPABILITIES_H

struct SolverCapabilities {
    bool isNoOp                 = false;
    bool supportsShiftedGreens = false;
    bool providesPotential     = false;
    bool providesGradient      = false;
    bool requiresPotentialBCs  = false;
    bool usesSeparatePotentialField = false;

    // These preserve the legacy OPALX_FIELD_DEBUG dump schedule independently
    // from the physical result fields exposed by the backend.
    bool debugDumpRhoBeforeSolve    = false;
    bool debugDumpScalarAfterSolve  = false;
    bool debugDumpVectorAfterSolve  = false;
};

#endif  // OPAL_SOLVER_CAPABILITIES_H
