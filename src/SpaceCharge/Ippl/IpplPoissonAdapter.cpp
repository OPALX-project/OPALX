/**
 * @file IpplPoissonAdapter.cpp
 * @brief Implements typed IPPL Poisson backend construction and dispatch.
 */

#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

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

    std::unique_ptr<IpplPoissonAdapter> IpplPoissonAdapter::create(
            PoissonBackendKind kind, GreenFunctionKind greenFunction,
            Solver_t<double, 3>& solverStorage, IpplPoissonFields fields) {
        std::unique_ptr<detail::IpplPoissonBackend> backend;
        switch (kind) {
            case PoissonBackendKind::None:
                backend = std::make_unique<detail::NullIpplPoissonBackend>(solverStorage, fields);
                break;
            case PoissonBackendKind::FftPeriodic:
                backend = std::make_unique<detail::FftIpplPoissonBackend>(solverStorage, fields);
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
        return std::unique_ptr<IpplPoissonAdapter>(new IpplPoissonAdapter(std::move(backend)));
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

    IpplPoissonAdapter::IpplPoissonAdapter(std::unique_ptr<detail::IpplPoissonBackend> backend)
        : backend_m(std::move(backend)) {}

    IpplPoissonAdapter::~IpplPoissonAdapter() = default;

    void IpplPoissonAdapter::solve(const IpplPoissonSolveRequest& request) {
        if (request.hasShiftedGreenFunction()
            && !backend_m->capabilities().supportsShiftedGreenFunction) {
            throw OpalException(
                    "IpplPoissonAdapter::solve",
                    "The selected backend does not support shifted Green functions.");
        }
        backend_m->solve(request);
    }

    void IpplPoissonAdapter::refresh(IpplPoissonFields fields) { backend_m->refresh(fields); }

    const IpplPoissonCapabilities& IpplPoissonAdapter::capabilities() const {
        return backend_m->capabilities();
    }

    double IpplPoissonAdapter::couplingConstant() const { return backend_m->couplingConstant(); }

}  // namespace opalx::spacecharge
