/**
 * @file CorrectionPlan.tpp
 * @brief Implements ordered Cartesian PIC correction-pass planning.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_CORRECTION_PLAN_TPP
#define OPALX_SPACE_CHARGE_PIC_CORRECTION_PLAN_TPP

#include "Utilities/OpalException.h"

#include <utility>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    CorrectionPlan<T, Dim>::CorrectionPlan(
            CorrectionConfig config, PoissonBackendKind backend, IterationKind iteration)
        : config_m(std::move(config)), backend_m(backend), iteration_m(iteration) {
        switch (config_m.kind()) {
            case CorrectionKind::None:
            case CorrectionKind::ImageCharge:
                break;
            case CorrectionKind::ShiftedGreen:
                if (backend_m != PoissonBackendKind::Open) {
                    throw OpalException(
                            "CorrectionPlan::CorrectionPlan",
                            "The shifted-Green correction requires the Open Poisson backend.");
                }
                break;
            default:
                throw OpalException(
                        "CorrectionPlan::CorrectionPlan",
                        "The configured correction kind is invalid.");
        }

        if (iteration_m != IterationKind::WholeBunch && iteration_m != IterationKind::Binning) {
            throw OpalException("CorrectionPlan::CorrectionPlan", "The iteration kind is invalid.");
        }
    }

    template <typename T, unsigned Dim>
    typename CorrectionPlan<T, Dim>::Prepared CorrectionPlan<T, Dim>::prepare(
            const RequestedPhysics& requested, std::size_t step) const {
        const bool withinStepBudget =
                config_m.maximumSteps() == 0 || step < config_m.maximumSteps();
        const bool correctionActive = config_m.enabled() && withinStepBudget;
        validateRequest(requested, correctionActive);

        Prepared prepared;
        prepared.configuredCorrection = {config_m.kind(), config_m.planeZ()};
        prepared.maximumSteps         = config_m.maximumSteps();
        prepared.correctionExpired    = config_m.enabled() && !withinStepBudget;
        if (correctionActive) {
            prepared.activeCorrection = requested.correction;
        }

        const ImagePolicy disabledImage{};
        const ImagePolicy enabledImage{true, config_m.planeZ()};

        if (iteration_m == IterationKind::WholeBunch) {
            const bool imageActive  = requested.correction.kind == CorrectionKind::ImageCharge;
            prepared.shiftedIgnored = requested.correction.kind == CorrectionKind::ShiftedGreen;

            append(prepared,
                   makePass(
                           SolvePassKind::PrimaryAndImage,
                           imageActive ? DepositKind::PrimaryAndImage : DepositKind::Primary,
                           imageActive ? enabledImage : disabledImage, BackendSolveRule::Standard,
                           false, FieldOutputRule::DirectGather, 1.0, FieldSourceRule::Direct,
                           imageActive));
            return prepared;
        }

        const bool correctionPassActive = requested.correction.kind != CorrectionKind::None;
        append(prepared,
               makePass(
                       SolvePassKind::Primary,
                       correctionPassActive ? DepositKind::Primary : DepositKind::PrimaryAndImage,
                       disabledImage, BackendSolveRule::Standard, true,
                       FieldOutputRule::LorentzAccumulation, 1.0, FieldSourceRule::Direct, false));

        if (requested.correction.kind == CorrectionKind::ImageCharge) {
            append(prepared,
                   makePass(
                           SolvePassKind::ImageCharge, DepositKind::Image, enabledImage,
                           BackendSolveRule::Standard, true, FieldOutputRule::LorentzAccumulation,
                           -1.0, FieldSourceRule::Direct, true));
        } else if (requested.correction.kind == CorrectionKind::ShiftedGreen) {
            append(prepared, makePass(
                                     SolvePassKind::ShiftedGreen, DepositKind::Primary,
                                     disabledImage, BackendSolveRule::ShiftedGreen, false,
                                     FieldOutputRule::LorentzAccumulation, -1.0,
                                     FieldSourceRule::ShiftedGreenImageZ, false));
        }

        return prepared;
    }

    template <typename T, unsigned Dim>
    typename CorrectionPlan<T, Dim>::Pass CorrectionPlan<T, Dim>::makePass(
            SolvePassKind kind, DepositKind depositKind, ImagePolicy imagePolicy,
            BackendSolveRule backendRule, bool suppressFieldDump, FieldOutputRule outputRule,
            double magneticSign, FieldSourceRule sourceRule, bool dumpDirichletPlaneAfter) const {
        Pass pass;
        pass.kind                    = kind;
        pass.depositKind             = depositKind;
        pass.imagePolicy             = imagePolicy;
        pass.backendRule             = backendRule;
        pass.suppressFieldDump       = suppressFieldDump;
        pass.outputRule              = outputRule;
        pass.magneticSign            = magneticSign;
        pass.sourceRule              = sourceRule;
        pass.planeZ                  = config_m.planeZ();
        pass.dumpDirichletPlaneAfter = dumpDirichletPlaneAfter;
        return pass;
    }

    template <typename T, unsigned Dim>
    void CorrectionPlan<T, Dim>::append(Prepared& prepared, const Pass& pass) const {
        if (prepared.passCount >= prepared.maximumPassCount) {
            throw OpalException(
                    "CorrectionPlan::prepare", "The correction pass capacity was exceeded.");
        }
        prepared.passes[prepared.passCount] = pass;
        ++prepared.passCount;
    }

    template <typename T, unsigned Dim>
    void CorrectionPlan<T, Dim>::validateRequest(
            const RequestedPhysics& requested, bool correctionActive) const {
        const bool expectsBinning = iteration_m == IterationKind::Binning;
        if (requested.useBinning != expectsBinning) {
            throw OpalException(
                    "CorrectionPlan::prepare",
                    "The requested binning mode differs from the immutable iteration plan.");
        }

        const CorrectionKind expectedKind =
                correctionActive ? config_m.kind() : CorrectionKind::None;
        if (requested.correction.kind != expectedKind) {
            throw OpalException(
                    "CorrectionPlan::prepare",
                    "The requested correction differs from the step-resolved configuration.");
        }
        if (correctionActive && requested.correction.planeZ != config_m.planeZ()) {
            throw OpalException(
                    "CorrectionPlan::prepare",
                    "The requested correction plane differs from the immutable configuration.");
        }

        const bool expectsPotential = correctionActive
                                      && config_m.kind() == CorrectionKind::ImageCharge
                                      && config_m.planeDumpFrequency() != 0;
        if (requested.writePotential != expectsPotential) {
            throw OpalException(
                    "CorrectionPlan::prepare",
                    "The requested potential output differs from the correction schedule.");
        }
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC_CORRECTION_PLAN_TPP
