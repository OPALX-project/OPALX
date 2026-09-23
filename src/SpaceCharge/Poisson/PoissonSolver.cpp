/**
 * @file PoissonSolver.cpp
 * @brief Implements shared 3D Poisson request handling, diagnostics and lifecycle.
 */

#include "SpaceCharge/Poisson/PoissonSolver.h"

#include "SpaceCharge/Poisson/NullPoissonAdapter.h"
#include "SpaceCharge/Poisson/OpenPoissonAdapter.h"
#include "SpaceCharge/Poisson/P3MAdapters.h"
#include "SpaceCharge/Poisson/PeriodicPoissonAdapter.h"

#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <Kokkos_Core.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

namespace opalx::spacecharge {
    namespace {

        void requireCommonFields(const PoissonFieldBinding& fields, const char* where) {
            if (fields.chargeDensity == nullptr || fields.electricField == nullptr) {
                throw OpalException(where, "Charge-density and electric fields must be bound.");
            }
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

    PoissonSolver::PoissonSolver(
            PoissonSolverConfig config, PoissonFieldBinding fields, PoissonSolverType expectedType)
        : config_m(std::move(config)), fields_m(fields) {
        validatePoissonSolverConfig(config_m);
        requireCommonFields(fields_m, "PoissonSolver::PoissonSolver");
        if (config_m.type != expectedType) {
            throw OpalException(
                    "PoissonSolver::PoissonSolver",
                    "The Poisson configuration does not match the concrete adapter type.");
        }
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

        solveImpl(request);

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
        // Finish work using the old extents before rebuilding fields and FFT plans.
        Kokkos::fence();
        rebuildImpl(fields);
        fields_m = fields;
    }

    std::unique_ptr<PoissonSolver> makePoissonSolver(
            PoissonSolverConfig config, PoissonFieldBinding fields) {
        validatePoissonSolverConfig(config);
        switch (config.type) {
            case PoissonSolverType::None:
                return std::make_unique<NullPoissonAdapter>(std::move(config), fields);
            case PoissonSolverType::PeriodicFFT:
                return std::make_unique<PeriodicPoissonAdapter>(std::move(config), fields);
            case PoissonSolverType::Open:
                return std::make_unique<OpenPoissonAdapter>(std::move(config), fields);
            case PoissonSolverType::P3M:
                return std::make_unique<P3MMeshPoissonAdapter>(std::move(config), fields);
            case PoissonSolverType::ConjugateGradient:
                break;
        }
        throw OpalException(
                "makePoissonSolver", "No implemented Poisson adapter matches the configuration.");
    }

}  // namespace opalx::spacecharge
