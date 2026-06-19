#ifndef OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_TPP
#define OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_TPP

#include "PartBunch/FieldOps.hpp"

#include <vector>

template <typename T, unsigned Dim>
void BinnedLorentzSelfFieldDriver<T, Dim>::compute(
        typename BinnedLorentzSelfFieldDriver<T, Dim>::Solver_t& solver,
        typename BinnedLorentzSelfFieldDriver<T, Dim>::PartBunch_t& bunch,
        const BoundaryCorrection<T, Dim>* activeCorrection) {
    using AdaptBins_t    = typename Solver_t::AdaptBins_t;
    using bin_index_type = typename Solver_t::bin_index_type;
    using size_type      = typename Solver_t::size_type;
    using BinStatsRow    = typename Solver_t::BinStatsRow;

    std::shared_ptr<AdaptBins_t> bins = bunch.getBins();
    if (!bins) {
        solver.monolithicDriver_m->compute(solver, bunch, activeCorrection);
        return;
    }

    solver.rebinAndPrepare(bunch, bins);

    std::shared_ptr<VField_t<T, Dim>> EtmpSP = bunch.getTempEField();
    if (!EtmpSP) {
        throw OpalException(
                "BinnedLorentzSelfFieldDriver::compute",
                "Temporary E field (Etmp) is not initialized.");
    }
    std::shared_ptr<VField_t<T, Dim>> BtmpSP = bunch.getTempBField();
    if (!BtmpSP) {
        throw OpalException(
                "BinnedLorentzSelfFieldDriver::compute",
                "Temporary B field (Btmp) is not initialized.");
    }

    solver.fieldAccumulator_m.bind(EtmpSP, BtmpSP);
    solver.fieldAccumulator_m.clear();

    const bin_index_type nBins = bins->getCurrentBinCount();

    Inform m("BinnedLorentzSelfFieldDriver::compute");
    m << level4 << "Binned mode: nBins=" << static_cast<int>(nBins)
      << ", stype=" << solver.getStype() << endl;

    std::vector<BinStatsRow> binStats;
    binStats.reserve(static_cast<size_t>(nBins));

    bool dumpedDirichletPlaneThisStep = false;

    for (bin_index_type binIndex = 0; binIndex < nBins; ++binIndex) {
        const size_type nPartGlobal = bins->getNPartInBin(binIndex, true);
        if (nPartGlobal == 0) {
            continue;
        }

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

        const typename Solver_t::BinKinematics kinematics =
                solver.computeGammaBinGlobal(bunch, bins, binIndex, nPartGlobal);
        if (kinematics.gammaBin <= 0.0) {
            throw OpalException(
                    "BinnedLorentzSelfFieldDriver::compute",
                    "Computed non-positive gamma for bin.");
        }

        m << level4 << "binIndex=" << static_cast<int>(binIndex)
          << " gammaBin=" << std::setprecision(10) << kinematics.gammaBin << endl;

        binStats.push_back(
                BinStatsRow{
                        static_cast<long long>(binIndex),
                        static_cast<unsigned long long>(nPartGlobal), kinematics.gammaBin});

        auto& mesh        = solver.getRho()->get_mesh();
        const auto hrOrig = mesh.getMeshSpacing();
        auto hrStretched  = hrOrig;
        hrStretched[Dim - 1] *= kinematics.gammaBin;

        BinContext<T, Dim> context;
        context.binIndex            = static_cast<long long>(binIndex);
        context.particleCountGlobal = static_cast<std::uint64_t>(nPartGlobal);
        context.gamma               = kinematics.gammaBin;
        context.pmean               = kinematics.pmean;
        context.meshSpacing         = hrOrig;
        context.solveSpacing        = hrStretched;
        context.imageCorrectionActive =
                activeCorrection
                && activeCorrection->kind() == BoundaryCorrectionKind::ImageCharge;
        context.shiftedGreensCorrectionActive =
                activeCorrection
                && activeCorrection->kind() == BoundaryCorrectionKind::ShiftedGreens;

        {
            const ImageScatterMode scatterMode = solver.primaryScatterModeForStep(activeCorrection);
            solver.prepareRhoForBin(bunch, bins, context, scatterMode);

            opalx::fieldops::setVectorField(*(solver.getE()), Vector_t<T, Dim>(0.0));
            mesh.setMeshSpacing(context.solveSpacing);

            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) start"
              << " (hr_z stretched by gamma=" << context.gamma << ")" << endl;
            solver.runSolver(true);
            m << level4 << "binIndex=" << static_cast<int>(binIndex)
              << " primary runSolver(true) done; accumulate->Etmp" << endl;

            solver.fieldAccumulator_m.accumulate(
                    *(solver.getE()), context.gamma, context.pmean, +1.0);

            mesh.setMeshSpacing(context.meshSpacing);
        }

        if (activeCorrection) {
            const auto correctionPass = activeCorrection->makePass(
                    context, mesh.getOrigin(),
                    solver.getRho()->getLayout().getDomain()[Dim - 1].length());
            if (correctionPass) {
                solver.applyBoundaryCorrectionPass(
                        bunch, bins, context, *correctionPass, dumpedDirichletPlaneThisStep);
            }
        }
    }

    solver.gatherFromTempToParticles(bunch);

    if (solver.tablePrintFrequency_m > 0) {
        const long long step = bunch.getGlobalTrackStep();
        if (step >= 0 && (step % solver.tablePrintFrequency_m) == 0) {
            solver.printBinStatsTable(bins->getBinningCmdName(), binStats);
        }
    }
}

#endif  // OPALX_BINNED_LORENTZ_SELF_FIELD_DRIVER_TPP
