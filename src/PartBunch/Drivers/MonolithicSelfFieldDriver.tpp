#ifndef OPALX_MONOLITHIC_SELF_FIELD_DRIVER_TPP
#define OPALX_MONOLITHIC_SELF_FIELD_DRIVER_TPP

#include "PartBunch/FieldOps.hpp"

#include <functional>
#include <numeric>

template <typename T, unsigned Dim>
void MonolithicSelfFieldDriver<T, Dim>::compute(
        typename MonolithicSelfFieldDriver<T, Dim>::Solver_t& solver,
        typename MonolithicSelfFieldDriver<T, Dim>::PartBunch_t& bunch,
        const BoundaryCorrection<T, Dim>* activeCorrection) {
    using ParticleCtr_t = typename Solver_t::ParticleCtr_t;

    std::shared_ptr<ParticleCtr_t> pc = bunch.getParticleContainer();
    if (!pc) {
        throw OpalException(
                "MonolithicSelfFieldDriver::compute",
                "Bunch particle container is not available.");
    }

    Inform m("MonolithicSelfFieldDriver::compute");
    m << level4 << "Legacy mode entry: localP=" << pc->getLocalNum()
      << ", totalP=" << pc->getTotalNum() << ", stype=" << solver.getStype() << endl;

    if (solver.scatterAttribute_m != Solver_t::ScatterAttribute::ChargeQ) {
        throw OpalException(
                "MonolithicSelfFieldDriver::compute",
                "Unsupported scatter attribute in legacy solver.");
    }

    typename Solver_t::particle_position_type* R = &pc->R;

    // Legacy step 1: scatter the whole bunch to rho in one mesh solve.
    Field_t<Dim>& rho = *(solver.getRho());
    opalx::fieldops::setScalarField(rho, 0.0);

    if (activeCorrection && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge) {
        solver.imageScatterController_m.scatterPrimaryAndImage(pc, *R, rho);
    } else {
        solver.imageScatterController_m.scatterPrimaryOnly(pc, *R, rho);
    }

    // Legacy step 2: apply the same rho normalization policy as the binned path.
    double normalizer = bunch.getdT();
    if (solver.config_m.normalizeRhoByCellVolume) {
        const double cellVolume =
                std::reduce(bunch.hr_m.begin(), bunch.hr_m.end(), 1.0, std::multiplies<double>());
        normalizer *= cellVolume;
    }

    double shift = 0.0;
    if (solver.config_m.subtractNeutralizingBackground) {
        double size = 1.0;
        for (size_t d = 0; d < Dim; ++d) {
            size *= bunch.rmax_m[d] - bunch.rmin_m[d];
        }

        const double totalQ = bunch.getParticleContainer()->getTotalCharge();
        shift               = -(totalQ / size) * solver.getCouplingConstant();
    }

    opalx::fieldops::scaleAndShiftScalarField(
            rho, solver.getCouplingConstant() / normalizer, shift);

    // Legacy step 3: run the Poisson backend once and gather E directly to particles.
    opalx::fieldops::setVectorField(*(solver.getE()), Vector_t<T, Dim>(0.0));

    m << level4 << "Legacy mode: runSolver() start" << endl;
    solver.runSolver();
    if (activeCorrection && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge) {
        solver.dumpDirichletPlaneDiagnosticsIfRequested(bunch, "legacy");
    }
    m << level4 << "Legacy mode: gather E->particles" << endl;

    if (solver.gatherAttribute_m == Solver_t::GatherAttribute::ElectricFieldE) {
        gather(pc->E, *solver.getE(), *R);
    } else {
        throw OpalException(
                "MonolithicSelfFieldDriver::compute",
                "Unsupported gather attribute in legacy solver.");
    }
}

#endif  // OPALX_MONOLITHIC_SELF_FIELD_DRIVER_TPP
