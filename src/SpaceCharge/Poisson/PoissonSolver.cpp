/**
 * @file PoissonSolver.cpp
 * @brief Implements variant-based IPPL Poisson backend dispatch.
 */

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        const PoissonSolverCapabilities nullCapabilities{
                .isNoOp                         = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true};

        const PoissonSolverCapabilities fftCapabilities{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        const PoissonSolverCapabilities openCapabilities{
                .supportsShiftedGreenFunction   = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        const PoissonSolverCapabilities p3mCapabilities{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        void requireCommonFields(const PoissonFieldBinding& fields, const char* where) {
            if (fields.chargeDensity == nullptr || fields.electricField == nullptr) {
                throw OpalException(where, "Charge-density and electric fields must be bound.");
            }
        }

        ippl::ParameterList commonFftParameters() {
            ippl::ParameterList parameters;
            parameters.add("use_heffte_defaults", false);
            parameters.add("use_pencils", true);
            parameters.add("use_reorder", false);
            parameters.add("use_gpu_aware", true);
            parameters.add("comm", ippl::p2p_pl);
            parameters.add("r2c_direction", 0);
            return parameters;
        }

        int openGreenFunctionValue(GreenFunctionType greenFunction) {
            switch (greenFunction) {
                case GreenFunctionType::Standard:
                    return OpenSolver_t<double, 3>::STANDARD;
                case GreenFunctionType::Integrated:
                    return OpenSolver_t<double, 3>::INTEGRATED;
            }
            throw OpalException(
                    "PoissonSolver::PoissonSolver", "Unknown open-solver Green function.");
        }

        template <typename NativeBackend>
        void bindFields(NativeBackend& backend, PoissonFieldBinding fields) {
            // IPPL requires the deposited charge field before the output field. Keep this order
            // during initial construction and every layout-driven reconstruction.
            backend.setRhs(*fields.chargeDensity);
            backend.setLhs(*fields.electricField);
        }

#ifdef OPALX_FIELD_DEBUG
        void dumpVectorField(VField_t<double, 3>& field, std::string what, std::size_t solveIndex) {
            Inform m("PoissonSolver::dumpVectorField");
            if (ippl::Comm->size() > 1 || solveIndex < 2) {
                m << level5 << "Skipping vector field dump for multiple ranks or first call."
                  << endl;
                return;
            }

            std::string type;
            std::string unit;
            if (Util::toUpper(what) == "EF") {
                type = "vector";
            }

            auto localIdx  = field.getOwned();
            auto* mesh     = &field.get_mesh();
            auto spacing   = mesh->getMeshSpacing();
            auto origin    = mesh->getOrigin();
            const int halo = field.getNghost();
            auto view      = field.getView();
            auto hostView  = field.getHostMirror();
            Kokkos::deep_copy(hostView, view);

            std::filesystem::path file("data/");
            std::ostringstream name;
            name << OpalData::getInstance()->getInputBasename() << "-" << what << "_" << type << "-"
                 << std::setfill('0') << std::setw(6) << solveIndex << ".dat";
            file /= name.str();
            std::ofstream output(file.string(), std::ios::out);
            output << std::setprecision(9);
            output << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
                   << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing
                   << std::endl
                   << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
                   << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17)
                   << "z [m]" << std::setw(10) << what << "x [" << unit << "]" << std::setw(10)
                   << what << "y [" << unit << "]" << std::setw(10) << what << "z [" << unit << "]"
                   << std::endl;

            for (int i = localIdx[0].first() + halo; i <= localIdx[0].last() + halo; ++i) {
                for (int j = localIdx[1].first() + halo; j <= localIdx[1].last() + halo; ++j) {
                    for (int k = localIdx[2].first() + halo; k <= localIdx[2].last() + halo; ++k) {
                        const double x = (i - halo) * spacing[0] + origin[0];
                        const double y = (j - halo) * spacing[1] + origin[1];
                        const double z = (k - halo) * spacing[2] + origin[2];
                        output << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                               << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                               << std::scientific << "\t" << hostView(i, j, k)[0] << "\t"
                               << hostView(i, j, k)[1] << "\t" << hostView(i, j, k)[2] << std::endl;
                    }
                }
            }
            output.close();
            m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to "
              << file.string() << endl;
        }

        void dumpScalarField(Field_t<3>& field, std::string what, std::size_t solveIndex) {
            Inform m("PoissonSolver::dumpScalarField");
            m << level5 << "Dumping scalar field: " << what << endl;
            if (ippl::Comm->size() > 1) {
                m << level5 << "Skipping scalar field dump for multiple ranks or first call."
                  << endl;
                return;
            }

            std::string type;
            std::string unit;
            if (Util::toUpper(what) == "RHO") {
                type = "scalar";
                unit = "Cb/m^3";
            } else if (Util::toUpper(what) == "PHI") {
                type = "scalar";
                unit = "V";
            }

            const ippl::NDIndex<3> localIdx = field.getLayout().getLocalNDIndex();
            const int halo                  = field.getNghost();
            auto* mesh                      = &field.get_mesh();
            auto spacing                    = mesh->getMeshSpacing();
            auto origin                     = mesh->getOrigin();
            auto view                       = field.getView();
            auto hostView                   = field.getHostMirror();
            Kokkos::deep_copy(hostView, view);

            std::filesystem::path file("data/");
            std::ostringstream name;
            name << OpalData::getInstance()->getInputBasename() << "-" << what << "_" << type << "-"
                 << std::setfill('0') << std::setw(6) << solveIndex << ".dat";
            file /= name.str();
            std::ofstream output(file.string(), std::ios::out);
            output << std::setprecision(9);
            output << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
                   << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing
                   << " nghosts=" << halo << std::endl
                   << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
                   << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17)
                   << "z [m]" << std::setw(13) << what << " [" << unit << "]" << std::endl;

            for (int i = localIdx[0].first() + halo; i <= localIdx[0].last() + halo; ++i) {
                for (int j = localIdx[1].first() + halo; j <= localIdx[1].last() + halo; ++j) {
                    for (int k = localIdx[2].first() + halo; k <= localIdx[2].last() + halo; ++k) {
                        const double x = (i - halo) * spacing[0] + origin[0];
                        const double y = (j - halo) * spacing[1] + origin[1];
                        const double z = (k - halo) * spacing[2] + origin[2];
                        output << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                               << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                               << std::scientific << "\t" << hostView(i, j, k) << std::endl;
                    }
                }
            }
            output.close();
            m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to "
              << file.string() << endl;
        }
#endif

    }  // namespace

    PoissonSolver::PoissonSolver(PoissonSolverConfig config, PoissonFieldBinding fields)
        : config_m(std::move(config)), fields_m(fields) {
        requireCommonFields(fields_m, "PoissonSolver::PoissonSolver");
        constructBackend();
        bindBackendFields(fields_m);
    }

    void PoissonSolver::constructBackend() {
        switch (config_m.type) {
            case PoissonSolverType::None: {
                auto& backend = backend_m.emplace<NullBackend>();
                backend.mergeParameters(ippl::ParameterList{});
                break;
            }
            case PoissonSolverType::PeriodicFFT: {
                auto& backend                  = backend_m.emplace<PeriodicBackend>();
                ippl::ParameterList parameters = commonFftParameters();
                parameters.add("output_type", PeriodicBackend::GRAD);
                backend.mergeParameters(parameters);
                break;
            }
            case PoissonSolverType::Open: {
                auto& backend                  = backend_m.emplace<OpenBackend>();
                ippl::ParameterList parameters = commonFftParameters();
                parameters.add("output_type", OpenBackend::SOL_AND_GRAD);
                parameters.add("algorithm", OpenBackend::HOCKNEY);
                parameters.add("greens_function", openGreenFunctionValue(config_m.greenFunction));
                backend.mergeParameters(parameters);
                break;
            }
            case PoissonSolverType::P3M: {
                if (!(config_m.p3mCutoff > 0.0)) {
                    throw OpalException(
                            "PoissonSolver::constructBackend",
                            "The P3M cutoff radius must be positive.");
                }
                const bool allPeriodic = std::all_of(
                        config_m.boundaryConditions.begin(), config_m.boundaryConditions.end(),
                        [](FieldBoundaryCondition kind) {
                            return kind == FieldBoundaryCondition::Periodic;
                        });
                const bool allOpen = std::all_of(
                        config_m.boundaryConditions.begin(), config_m.boundaryConditions.end(),
                        [](FieldBoundaryCondition kind) {
                            return kind == FieldBoundaryCondition::Open;
                        });
                if (!allPeriodic && !allOpen) {
                    throw OpalException(
                            "PoissonSolver::constructBackend",
                            "P3M requires uniform OPEN or PERIODIC boundary conditions.");
                }
                auto& backend                  = backend_m.emplace<P3MBackend>();
                ippl::ParameterList parameters = commonFftParameters();
                parameters.add("output_type", P3MBackend::GRAD);
                parameters.add("alpha", 2.0 / config_m.p3mCutoff);
                parameters.add("force_constant", -1.0 / (4.0 * Physics::pi));
                parameters.add("regularization_cutoff", 1.0e-9);
                parameters.add(
                        "boundary_type", allPeriodic ? P3MBackend::PERIODIC : P3MBackend::OPEN);
                backend.mergeParameters(parameters);
                break;
            }
            case PoissonSolverType::ConjugateGradient:
                throw OpalException(
                        "PoissonSolver::constructBackend",
                        "The CG Poisson backend is recognized but not implemented.");
        }
        throw OpalException(
                "PoissonSolver::constructBackend",
                "No IPPL backend matches the Poisson configuration.");
    }

    void PoissonSolver::bindBackendFields(PoissonFieldBinding fields) {
        std::visit(
                [&fields](auto& backend) {
                    using Backend = std::decay_t<decltype(backend)>;
                    if constexpr (std::is_same_v<Backend, std::monostate>) {
                        throw OpalException(
                                "PoissonSolver::bindBackendFields",
                                "The Poisson backend has not been constructed.");
                    } else {
                        bindFields(backend, fields);
                    }
                },
                backend_m);
    }

    void PoissonSolver::solve(
            const PoissonSolveRequest& request, const PoissonSolveOptions& options) {
        Inform m("PoissonSolver::solve");
        m << level3 << "Running solver with type: " << name()
          << ". Force skip field dump: " << options.suppressFieldDump << endl;

        if (request.hasShiftedGreenFunction() && !capabilities().supportsShiftedGreenFunction) {
            throw OpalException(
                    "PoissonSolver::solve",
                    "The selected backend does not support shifted Green functions.");
        }

        [[maybe_unused]] const std::size_t solveIndex                         = runtimeSolveCount_m;
        [[maybe_unused]] const PoissonSolverCapabilities& backendCapabilities = capabilities();
#ifdef OPALX_FIELD_DEBUG
        if (!options.suppressFieldDump && backendCapabilities.debugDumpChargeBeforeSolve) {
            dumpScalarField(*fields_m.chargeDensity, "rho", solveIndex);
        }
#endif

        std::visit(
                [&request](auto& backend) {
                    using Backend = std::decay_t<decltype(backend)>;
                    if constexpr (std::is_same_v<Backend, std::monostate>) {
                        throw OpalException(
                                "PoissonSolver::solve",
                                "The Poisson backend has not been constructed.");
                    } else if constexpr (std::is_same_v<Backend, OpenBackend>) {
                        if (request.hasShiftedGreenFunction()) {
                            backend.shiftedGreensFunction(*request.greenFunctionShift);
                            backend.solve();
                            backend.greensFunction();
                        } else {
                            backend.solve();
                        }
                    } else {
                        backend.solve();
                    }
                },
                backend_m);

#ifdef OPALX_FIELD_DEBUG
        if (!options.suppressFieldDump && backendCapabilities.debugDumpScalarAfterSolve) {
            dumpScalarField(*fields_m.chargeDensity, "phi", solveIndex);
        }
        if (!options.suppressFieldDump && backendCapabilities.debugDumpVectorAfterSolve) {
            dumpVectorField(*fields_m.electricField, "ef", solveIndex);
        }
#endif
        ++runtimeSolveCount_m;
    }

    void PoissonSolver::warmup() {
        // Exercise backend allocation and FFT planning without consuming runtime numbering.
        Kokkos::deep_copy(fields_m.chargeDensity->getView(), 0.0);
        solve({}, {.suppressFieldDump = true});
        runtimeSolveCount_m = 0;
    }

    void PoissonSolver::rebuildAfterLayoutChange(PoissonFieldBinding fields) {
        requireCommonFields(fields, "PoissonSolver::rebuildAfterLayoutChange");
        // Layout refresh can resize device fields while the previous FFT still owns buffers and
        // plans for the old extents. Complete that work before destroying the typed backend, then
        // reconstruct it so IPPL allocates matching internal fields as well as a matching plan.
        Kokkos::fence();
        constructBackend();
        bindBackendFields(fields);
        fields_m = fields;
    }

    std::string_view PoissonSolver::name() const {
        switch (config_m.type) {
            case PoissonSolverType::None:
                return "NONE";
            case PoissonSolverType::PeriodicFFT:
                return "FFT";
            case PoissonSolverType::Open:
                return "OPEN";
            case PoissonSolverType::P3M:
                return "P3M";
            case PoissonSolverType::ConjugateGradient:
                return "CG";
        }
        throw OpalException("PoissonSolver::name", "Unknown Poisson backend.");
    }

    const PoissonSolverCapabilities& PoissonSolver::capabilities() const {
        switch (config_m.type) {
            case PoissonSolverType::None:
                return nullCapabilities;
            case PoissonSolverType::PeriodicFFT:
                return fftCapabilities;
            case PoissonSolverType::Open:
                return openCapabilities;
            case PoissonSolverType::P3M:
                return p3mCapabilities;
            case PoissonSolverType::ConjugateGradient:
                break;
        }
        throw OpalException(
                "PoissonSolver::capabilities", "The CG Poisson backend is not implemented.");
    }

    double PoissonSolver::couplingConstant() const {
        if (config_m.type == PoissonSolverType::None) {
            return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
        }
        if (config_m.type == PoissonSolverType::PeriodicFFT
            || config_m.type == PoissonSolverType::Open
            || config_m.type == PoissonSolverType::P3M) {
            return 1.0 / Physics::epsilon_0;
        }
        throw OpalException(
                "PoissonSolver::couplingConstant", "The CG Poisson backend is not implemented.");
    }

}  // namespace opalx::spacecharge
