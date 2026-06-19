#ifndef OPAL_SOLVER_CAPABILITIES_H
#define OPAL_SOLVER_CAPABILITIES_H

/**
 * @brief Backend feature and policy flags consumed by OPALX field-solver orchestration.
 *
 * New Poisson backends should describe their storage, diagnostic, and density-normalization
 * requirements here so `PartBunch`, binned algorithms, and corrections do not branch on solver
 * names.
 */
struct SolverCapabilities {
    bool isNoOp                     = false;  //!< Backend intentionally performs no solve.
    bool supportsShiftedGreens      = false;  //!< Backend accepts `SolveRequest::greensShift`.
    bool providesPotential          = false;  //!< Backend exposes an electrostatic potential.
    bool providesGradient           = false;  //!< Backend exposes electric-field gradients.
    bool requiresPotentialBCs       = false;  //!< Backend needs potential field BCs applied.
    bool usesSeparatePotentialField = false;

    bool normalizeRhoByCellVolume = true;  //!< Include mesh cell volume in rho normalization.
    bool subtractNeutralizingBackground =
            true;  //!< Subtract the uniform background charge before solving.

    // These preserve the legacy OPALX_FIELD_DEBUG dump schedule independently from the physical
    // result fields exposed by the backend.
    bool debugDumpRhoBeforeSolve   = false;
    bool debugDumpScalarAfterSolve = false;
    bool debugDumpVectorAfterSolve = false;
};

#endif  // OPAL_SOLVER_CAPABILITIES_H
