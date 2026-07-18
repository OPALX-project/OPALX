/**
 * @file IpplPoissonAdapter.cpp
 * @brief Implements typed IPPL Poisson backend construction and dispatch.
 */

#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        void requireCommonFields(const IpplPoissonFields& fields, const char* where) {
            if (fields.chargeDensity == nullptr || fields.electricField == nullptr) {
                throw OpalException(where, "Charge-density and electric fields must be bound.");
            }
        }

        const IpplPoissonCapabilities nullCapabilities{
                .isNoOp                         = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true};

        const IpplPoissonCapabilities fftCapabilities{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = true,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        const IpplPoissonCapabilities openCapabilities{
                .supportsShiftedGreenFunction   = true,
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        const IpplPoissonCapabilities cgCapabilities{
                .usesSeparatePotentialField          = true,
                .requiresPotentialBoundaryConditions = true,
                .normalizeChargeByCellVolume         = true,
                .subtractNeutralizingBackground      = true,
                .debugDumpChargeBeforeSolve          = true,
                .debugDumpScalarAfterSolve           = true};

        int openGreenFunctionValue(GreenFunctionKind greenFunction) {
            switch (greenFunction) {
                case GreenFunctionKind::Standard:
                    return OpenSolver_t<double, 3>::STANDARD;
                case GreenFunctionKind::Integrated:
                    return OpenSolver_t<double, 3>::INTEGRATED;
            }
            throw OpalException(
                    "IpplPoissonAdapter::create", "Unknown open-solver Green function.");
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

        const char* backendName(PoissonBackendKind kind) {
            switch (kind) {
                case PoissonBackendKind::None:
                    return "NONE";
                case PoissonBackendKind::FftPeriodic:
                    return "FFT";
                case PoissonBackendKind::Open:
                    return "OPEN";
                case PoissonBackendKind::ConjugateGradient:
                    return "CG";
            }
            throw OpalException(
                    "IpplPoissonAdapter::backendName",
                    "No IPPL backend matches the configuration.");
        }

        std::unique_ptr<detail::IpplPoissonBackend> makeBackend(
                PoissonBackendKind kind, GreenFunctionKind greenFunction,
                Solver_t<double, 3>& solverStorage, IpplPoissonFields fields);

    }  // namespace

    namespace detail {

        class IpplPoissonBackend {
        public:
            virtual ~IpplPoissonBackend() = default;

            virtual void solve(const IpplPoissonSolveRequest& request)                = 0;
            virtual void refresh(IpplPoissonFields fields)                            = 0;
            [[nodiscard]] virtual const IpplPoissonCapabilities& capabilities() const = 0;
            [[nodiscard]] virtual double couplingConstant() const                     = 0;
        };

        class NullIpplPoissonBackend final : public IpplPoissonBackend {
        public:
            NullIpplPoissonBackend(Solver_t<double, 3>& solverStorage, IpplPoissonFields fields) {
                requireCommonFields(fields, "NullIpplPoissonBackend::NullIpplPoissonBackend");
                solverStorage.emplace<NullSolver_t<double, 3>>();
                solver_m = &std::get<NullSolver_t<double, 3>>(solverStorage);
                solver_m->mergeParameters(ippl::ParameterList{});
                refresh(fields);
            }

            void solve(const IpplPoissonSolveRequest& request) override {
                if (request.hasShiftedGreenFunction()) {
                    throw OpalException(
                            "NullIpplPoissonBackend::solve",
                            "The null backend does not support shifted Green functions.");
                }
                solver_m->solve();
            }

            void refresh(IpplPoissonFields fields) override {
                requireCommonFields(fields, "NullIpplPoissonBackend::refresh");
                solver_m->setRhs(*fields.chargeDensity);
                solver_m->setLhs(*fields.electricField);
            }

            const IpplPoissonCapabilities& capabilities() const override {
                return nullCapabilities;
            }

            double couplingConstant() const override {
                return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
            }

        private:
            NullSolver_t<double, 3>* solver_m = nullptr;
        };

        class FftIpplPoissonBackend final : public IpplPoissonBackend {
        public:
            FftIpplPoissonBackend(Solver_t<double, 3>& solverStorage, IpplPoissonFields fields) {
                requireCommonFields(fields, "FftIpplPoissonBackend::FftIpplPoissonBackend");
                solverStorage.emplace<FFTSolver_t<double, 3>>();
                solver_m = &std::get<FFTSolver_t<double, 3>>(solverStorage);

                ippl::ParameterList parameters = commonFftParameters();
                parameters.add("output_type", FFTSolver_t<double, 3>::GRAD);
                solver_m->mergeParameters(parameters);
                refresh(fields);
            }

            void solve(const IpplPoissonSolveRequest& request) override {
                if (request.hasShiftedGreenFunction()) {
                    throw OpalException(
                            "FftIpplPoissonBackend::solve",
                            "The periodic FFT backend does not support shifted Green functions.");
                }
                solver_m->solve();
            }

            void refresh(IpplPoissonFields fields) override {
                requireCommonFields(fields, "FftIpplPoissonBackend::refresh");
                solver_m->setRhs(*fields.chargeDensity);
                solver_m->setLhs(*fields.electricField);
            }

            const IpplPoissonCapabilities& capabilities() const override { return fftCapabilities; }

            double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

        private:
            FFTSolver_t<double, 3>* solver_m = nullptr;
        };

        class OpenIpplPoissonBackend final : public IpplPoissonBackend {
        public:
            OpenIpplPoissonBackend(
                    GreenFunctionKind greenFunction, Solver_t<double, 3>& solverStorage,
                    IpplPoissonFields fields) {
                requireCommonFields(fields, "OpenIpplPoissonBackend::OpenIpplPoissonBackend");
                solverStorage.emplace<OpenSolver_t<double, 3>>();
                solver_m = &std::get<OpenSolver_t<double, 3>>(solverStorage);

                ippl::ParameterList parameters = commonFftParameters();
                parameters.add("output_type", OpenSolver_t<double, 3>::SOL_AND_GRAD);
                parameters.add("algorithm", OpenSolver_t<double, 3>::HOCKNEY);
                parameters.add("greens_function", openGreenFunctionValue(greenFunction));
                solver_m->mergeParameters(parameters);
                refresh(fields);
            }

            void solve(const IpplPoissonSolveRequest& request) override {
                if (!request.hasShiftedGreenFunction()) {
                    solver_m->solve();
                    return;
                }

                try {
                    solver_m->shiftedGreensFunction(*request.greenFunctionShift);
                    solver_m->solve();
                } catch (...) {
                    try {
                        solver_m->greensFunction();
                    } catch (...) {
                        // Preserve the original backend failure.
                    }
                    throw;
                }
                solver_m->greensFunction();
            }

            void refresh(IpplPoissonFields fields) override {
                requireCommonFields(fields, "OpenIpplPoissonBackend::refresh");
                solver_m->setRhs(*fields.chargeDensity);
                solver_m->setLhs(*fields.electricField);
            }

            const IpplPoissonCapabilities& capabilities() const override {
                return openCapabilities;
            }

            double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

        private:
            OpenSolver_t<double, 3>* solver_m = nullptr;
        };

        class CgIpplPoissonBackend final : public IpplPoissonBackend {
        public:
            CgIpplPoissonBackend(Solver_t<double, 3>& solverStorage, IpplPoissonFields fields) {
                requireCommonFields(fields, "CgIpplPoissonBackend::CgIpplPoissonBackend");
                if (fields.potential == nullptr) {
                    throw OpalException(
                            "CgIpplPoissonBackend::CgIpplPoissonBackend",
                            "The CG backend requires a potential field.");
                }

                solverStorage.emplace<CGSolver_t<double, 3>>();
                solver_m = &std::get<CGSolver_t<double, 3>>(solverStorage);
                ippl::ParameterList parameters;
                parameters.add("output_type", CGSolver_t<double, 3>::GRAD);
                parameters.add("tolerance", 1e-12);
                solver_m->mergeParameters(parameters);
                solver_m->setRhs(*fields.chargeDensity);

                throw OpalException(
                        "CgIpplPoissonBackend::CgIpplPoissonBackend",
                        "Cannot use CGSolver yet, not fully implemented.");
            }

            void solve(const IpplPoissonSolveRequest&) override { solver_m->solve(); }

            void refresh(IpplPoissonFields fields) override {
                requireCommonFields(fields, "CgIpplPoissonBackend::refresh");
                if (fields.potential == nullptr) {
                    throw OpalException(
                            "CgIpplPoissonBackend::refresh",
                            "The CG backend requires a potential field.");
                }
                solver_m->setRhs(*fields.chargeDensity);
                solver_m->setLhs(*fields.potential);
                solver_m->setGradient(*fields.electricField);
            }

            const IpplPoissonCapabilities& capabilities() const override { return cgCapabilities; }

            double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

        private:
            CGSolver_t<double, 3>* solver_m = nullptr;
        };

    }  // namespace detail

    namespace {

        std::unique_ptr<detail::IpplPoissonBackend> makeBackend(
                PoissonBackendKind kind, GreenFunctionKind greenFunction,
                Solver_t<double, 3>& solverStorage, IpplPoissonFields fields) {
            std::unique_ptr<detail::IpplPoissonBackend> backend;
            switch (kind) {
                case PoissonBackendKind::None:
                    backend =
                            std::make_unique<detail::NullIpplPoissonBackend>(solverStorage, fields);
                    break;
                case PoissonBackendKind::FftPeriodic:
                    backend =
                            std::make_unique<detail::FftIpplPoissonBackend>(solverStorage, fields);
                    break;
                case PoissonBackendKind::Open:
                    backend = std::make_unique<detail::OpenIpplPoissonBackend>(
                            greenFunction, solverStorage, fields);
                    break;
                case PoissonBackendKind::ConjugateGradient:
                    backend = std::make_unique<detail::CgIpplPoissonBackend>(solverStorage, fields);
                    break;
            }
            if (backend == nullptr) {
                throw OpalException(
                        "IpplPoissonAdapter::create", "No IPPL backend matches the configuration.");
            }
            return backend;
        }

#ifdef OPALX_FIELD_DEBUG
        void dumpVectorField(VField_t<double, 3>& field, std::string what, std::size_t solveIndex) {
            Inform m("IpplPoissonAdapter::dumpVectorField");
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
            Inform m("IpplPoissonAdapter::dumpScalarField");
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

    std::unique_ptr<IpplPoissonAdapter> IpplPoissonAdapter::create(
            PoissonBackendKind kind, GreenFunctionKind greenFunction, IpplPoissonFields fields) {
        return std::unique_ptr<IpplPoissonAdapter>(
                new IpplPoissonAdapter(kind, greenFunction, fields));
    }

    const IpplPoissonCapabilities& IpplPoissonAdapter::capabilitiesFor(PoissonBackendKind kind) {
        switch (kind) {
            case PoissonBackendKind::None:
                return nullCapabilities;
            case PoissonBackendKind::FftPeriodic:
                return fftCapabilities;
            case PoissonBackendKind::Open:
                return openCapabilities;
            case PoissonBackendKind::ConjugateGradient:
                return cgCapabilities;
        }
        throw OpalException(
                "IpplPoissonAdapter::capabilitiesFor",
                "No IPPL backend matches the configuration.");
    }

    double IpplPoissonAdapter::couplingConstantFor(PoissonBackendKind kind) {
        switch (kind) {
            case PoissonBackendKind::None:
                return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
            case PoissonBackendKind::FftPeriodic:
            case PoissonBackendKind::Open:
            case PoissonBackendKind::ConjugateGradient:
                return 1.0 / Physics::epsilon_0;
        }
        throw OpalException(
                "IpplPoissonAdapter::couplingConstantFor",
                "No IPPL backend matches the configuration.");
    }

    IpplPoissonAdapter::IpplPoissonAdapter(
            PoissonBackendKind kind, GreenFunctionKind greenFunction, IpplPoissonFields fields)
        : backend_m(makeBackend(kind, greenFunction, solverStorage_m, fields)),
          fields_m(fields),
          kind_m(kind) {}

    IpplPoissonAdapter::~IpplPoissonAdapter() = default;

    void IpplPoissonAdapter::solve(
            const IpplPoissonSolveRequest& request, const IpplPoissonSolveOptions& options) {
        Inform m("IpplPoissonAdapter::solve");
        m << level3 << "Running solver with type: " << backendName(kind_m)
          << ". Force skip field dump: " << options.suppressFieldDump << endl;

        if (request.hasShiftedGreenFunction()
            && !backend_m->capabilities().supportsShiftedGreenFunction) {
            throw OpalException(
                    "IpplPoissonAdapter::solve",
                    "The selected backend does not support shifted Green functions.");
        }

        std::optional<SelfFieldDiagnostics::ScopedEvent> backendEvent;
        if (options.diagnostics != nullptr) {
            backendEvent.emplace(
                    *options.diagnostics, SelfFieldEventKind::BackendSolve,
                    request.hasShiftedGreenFunction() ? "OPEN-shifted" : "backend");
        }

        [[maybe_unused]] const std::size_t solveIndex =
                options.diagnostics != nullptr ? options.diagnostics->backendSolveCount() : 0;
        [[maybe_unused]] const IpplPoissonCapabilities& backendCapabilities = capabilities();
#ifdef OPALX_FIELD_DEBUG
        if (!options.suppressFieldDump && backendCapabilities.debugDumpChargeBeforeSolve) {
            dumpScalarField(*fields_m.chargeDensity, "rho", solveIndex);
        }
#endif

        backend_m->solve(request);

#ifdef OPALX_FIELD_DEBUG
        if (!options.suppressFieldDump && backendCapabilities.debugDumpScalarAfterSolve) {
            Field_t<3>& scalarOutput = backendCapabilities.usesSeparatePotentialField
                                               ? *fields_m.potential
                                               : *fields_m.chargeDensity;
            dumpScalarField(scalarOutput, "phi", solveIndex);
        }
        if (!options.suppressFieldDump && backendCapabilities.debugDumpVectorAfterSolve) {
            dumpVectorField(*fields_m.electricField, "ef", solveIndex);
        }
#endif

        if (options.diagnostics != nullptr) {
            options.diagnostics->completeBackendSolve();
        }
    }

    void IpplPoissonAdapter::warmup() {
        Kokkos::deep_copy(fields_m.chargeDensity->getView(), 0.0);
        solve({}, {.suppressFieldDump = true});
    }

    void IpplPoissonAdapter::setPotentialBoundaryConditions(
            const std::array<BoundaryConditionKind, 3>& boundaryConditions) {
        if (!capabilities().requiresPotentialBoundaryConditions) {
            return;
        }
        ippl::BConds<Field_t<3>, 3> translated;
        for (unsigned face = 0; face < 6; ++face) {
            switch (boundaryConditions[face / 2]) {
                case BoundaryConditionKind::Open:
                    translated[face] = std::make_shared<ippl::NoBcFace<Field_t<3>>>(face);
                    break;
                case BoundaryConditionKind::Dirichlet:
                    translated[face] = std::make_shared<ippl::ZeroFace<Field_t<3>>>(face);
                    break;
                case BoundaryConditionKind::Periodic:
                    translated[face] = std::make_shared<ippl::PeriodicFace<Field_t<3>>>(face);
                    break;
            }
        }
        setPotentialBoundaryConditions(std::move(translated));
    }

    void IpplPoissonAdapter::setPotentialBoundaryConditions(
            ippl::BConds<Field_t<3>, 3> boundaryConditions) {
        if (!capabilities().requiresPotentialBoundaryConditions) {
            return;
        }
        if (fields_m.potential == nullptr) {
            throw OpalException(
                    "IpplPoissonAdapter::setPotentialBoundaryConditions",
                    "The selected backend requires a potential field.");
        }
        fields_m.potential->setFieldBC(boundaryConditions);
    }

    void IpplPoissonAdapter::refresh(IpplPoissonFields fields) {
        backend_m->refresh(fields);
        fields_m = fields;
    }

    const IpplPoissonCapabilities& IpplPoissonAdapter::capabilities() const {
        return backend_m->capabilities();
    }

    double IpplPoissonAdapter::couplingConstant() const { return backend_m->couplingConstant(); }

}  // namespace opalx::spacecharge
