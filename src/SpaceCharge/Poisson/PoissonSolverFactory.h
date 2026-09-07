#ifndef OPALX_POISSON_SOLVER_FACTORY_H
#define OPALX_POISSON_SOLVER_FACTORY_H

#include "SpaceCharge/Poisson/PoissonSolver.h"

namespace opalx::spacecharge {

    /** @brief Validate configuration and construct the selected 3D adapter without warming it. */
    [[nodiscard]] std::unique_ptr<PoissonSolver> makePoissonSolver(
            PoissonSolverConfig config, PoissonFieldBinding fields);

    namespace detail {

        inline ippl::ParameterList commonFftParameters() {
            ippl::ParameterList parameters;
            parameters.add("use_heffte_defaults", false);
            parameters.add("use_pencils", true);
            parameters.add("use_reorder", false);
            parameters.add("use_gpu_aware", true);
            parameters.add("comm", ippl::p2p_pl);
            parameters.add("r2c_direction", 0);
            return parameters;
        }

        /** @brief Bind a fully constructed native solver, preserving IPPL's initialization order.
         */
        template <typename NativeBackend>
        void bindFields(NativeBackend& backend, PoissonFieldBinding fields) {
            // setRhs initializes native fields and CUDA/FFT resources; it must precede setLhs.
            backend.setRhs(*fields.chargeDensity);
            backend.setLhs(*fields.electricField);
        }

    }  // namespace detail

}  // namespace opalx::spacecharge

#endif  // OPALX_POISSON_SOLVER_FACTORY_H
