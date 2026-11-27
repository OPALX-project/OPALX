#include "Algorithms/SpaceChargeTransport.h"

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include "Algorithms/CompensatedSum.h"
#include "Algorithms/PartData.h"
#include "Utilities/OpalException.h"

namespace opalx::maps {
    namespace {
        using Chart = SpaceChargeChart;
        void require(bool condition, const char* message) {
            if (!condition) throw OpalException("SpaceChargeTransport", message);
        }
        double path(const LinearTransferMapReference& r) {
            return r.pathLength - r.pathLengthCorrection;
        }
        SpaceChargeCoordinates coordinates(
                const LinearTransferMapReference& station, const PartData& species,
                const ExternalFieldRayTracker& tracker) {
            SpaceChargeReference reference;
            reference.orbit   = station;
            reference.massEV  = species.getM();
            reference.chargeE = species.getQ();
            ExternalFieldRayTracker::State ray{
                    station.position,
                    station.momentum,
                    station.time,
                    station.pathLength,
                    station.positionCorrection,
                    station.timeCorrection,
                    station.pathLengthCorrection};
            Vector_t<double, 3> e(0), b(0);
            tracker.sampleExternalFields(ray, e, b);
            for (unsigned i = 0; i < 3; ++i) {
                reference.electricFieldLab[i] = e(i);
                reference.magneticFieldLab[i] = b(i);
            }
            return SpaceChargeCoordinates(reference);
        }
        tps::Map<double> canonicalHalf(
                const LinearTransferMapDraft& draft, const SpaceChargeCoordinates& in,
                const SpaceChargeCoordinates& out, const std::string& id) {
            return tps::compose(
                    out.map(Chart::SectionSlopes, Chart::CommonTimeCanonical, id + "/output"),
                    tps::compose(
                            toTpsLinearTransferMap(draft, id + "/external"),
                            in.map(Chart::CommonTimeCanonical, Chart::SectionSlopes,
                                   id + "/input")),
                    1);
        }
    }  // namespace

    SpaceChargeTransportPlan::SpaceChargeTransportPlan(
            SpaceChargeCoordinates entrance, SpaceChargeCoordinates exit,
            std::vector<SpaceChargeTransportSlice> slices, std::string id, std::size_t ns,
            double pathLength)
        : entrance_m(std::move(entrance)),
          exit_m(std::move(exit)),
          slices_m(std::move(slices)),
          id_m(std::move(id)),
          ns_m(ns),
          pathLength_m(pathLength) {}

    SpaceChargeTransportPlan buildSpaceChargeTransportPlan(
            OpalBeamline& beamline, const PartData& species, const TransferMapPass& pass, double dt,
            const SpaceChargeTransportSettings& settings, const std::string& occurrenceId) {
        require(std::isfinite(dt) && dt > 0, "External DT must be positive finite seconds.");
        require(std::isfinite(species.getM()) && species.getM() > 0
                        && std::isfinite(species.getQ()),
                "Invalid reference species.");
        require(pass.available() && !pass.segments().empty(),
                "A nonempty available reference pass is required.");
        require(!occurrenceId.empty(), "Empty transport provenance namespace.");
        require(settings.ns > 0 && settings.maximumSlices > 0
                        && settings.ns <= settings.maximumSlices / pass.segments().size()
                        && settings.ns <= std::numeric_limits<std::size_t>::max() / 2,
                "NS must be positive and the complete plan must fit maximumSlices.");
        auto external                  = settings.external;
        external.requireIsolatedFields = true;
        external.validate();
        const auto method =
                ExternalFieldRayTracker::integrationMethodName(external.integrationMethod);
        // Preflight all retained records before constructing any half-maps.
        const LinearTransferMapMetadata* previous = nullptr;
        for (const auto& handle : pass.segments()) {
            require(handle && handle->available(), "Unavailable source segment.");
            const auto& record = handle->metadata();
            require(!record.includesOverlappingFields,
                    "Overlapping elements/fields remain deferred.");
            require(record.integrationMethod == method,
                    "Reference-pass and half-map integrators differ.");
            require(std::isfinite(path(record.entrance)) && std::isfinite(path(record.exit))
                            && path(record.exit) > path(record.entrance),
                    "Nonpositive source interval.");
            if (previous) {
                require(linearTransferMapChart(previous->exit)
                                        == linearTransferMapChart(record.entrance)
                                && compensated::difference(
                                           record.entrance.pathLength,
                                           record.entrance.pathLengthCorrection,
                                           previous->exit.pathLength,
                                           previous->exit.pathLengthCorrection)
                                           == 0,
                        "Source pass has mismatched or noncontiguous boundaries.");
            }
            previous = &record;
        }
        LinearTransferMapBuilder builder(beamline, species, dt, external);
        ExternalFieldRayTracker tracker(beamline, species, external.integrationMethod, true);
        const auto& first     = pass.segments().front()->metadata().entrance;
        auto currentReference = first;
        auto current          = coordinates(first, species, tracker);
        const auto entrance   = current;
        const auto id = pass.id() + "/" + occurrenceId + "/NS:" + std::to_string(settings.ns);
        std::vector<SpaceChargeTransportSlice> slices;
        slices.reserve(settings.ns * pass.segments().size());
        double total = 0, correction = 0;
        try {
            for (std::size_t source = 0; source < pass.segments().size(); ++source) {
                const auto& record  = pass.segments()[source]->metadata();
                const double begin  = path(record.entrance);
                const double length = compensated::difference(
                        record.exit.pathLength, record.exit.pathLengthCorrection,
                        record.entrance.pathLength, record.entrance.pathLengthCorrection);
                for (std::size_t slice = 0; slice < settings.ns; ++slice) {
                    const double middle =
                            std::fma((double(slice) + 0.5) / double(settings.ns), length, begin);
                    const double end = slice + 1 == settings.ns
                                               ? path(record.exit)
                                               : std::fma(
                                                         double(slice + 1) / double(settings.ns),
                                                         length, begin);
                    require(middle > path(currentReference) && end > middle,
                            "NS subdivision is below path-coordinate resolution.");
                    const auto midpointReference =
                            builder.referenceAtPath(currentReference, record.exit, middle);
                    const auto exitReference =
                            slice + 1 == settings.ns
                                    ? record.exit
                                    : builder.referenceAtPath(midpointReference, record.exit, end);
                    auto midpoint  = coordinates(midpointReference, species, tracker);
                    auto exit      = coordinates(exitReference, species, tracker);
                    auto leftDraft = builder.externalInterval(
                            currentReference, midpointReference, 2 * slices.size());
                    auto rightDraft = builder.externalInterval(
                            midpointReference, exitReference, 2 * slices.size() + 1);
                    const auto prefix = id + "/slice:" + std::to_string(slices.size());
                    auto left      = canonicalHalf(leftDraft, current, midpoint, prefix + "/left");
                    auto right     = canonicalHalf(rightDraft, midpoint, exit, prefix + "/right");
                    const double h = compensated::difference(
                            exitReference.pathLength, exitReference.pathLengthCorrection,
                            currentReference.pathLength, currentReference.pathLengthCorrection);
                    require(std::isfinite(h) && h > 0, "Invalid SC integration length.");
                    compensated::add(h, total, correction);
                    slices.push_back(
                            {source, slice, path(currentReference), path(midpointReference),
                             path(exitReference), h, std::move(midpoint), std::move(leftDraft),
                             std::move(rightDraft), std::move(left), std::move(right)});
                    currentReference = exitReference;
                    current          = std::move(exit);
                }
            }
        } catch (const std::invalid_argument& error) {
            throw OpalException("buildSpaceChargeTransportPlan", error.what());
        }
        return SpaceChargeTransportPlan(
                entrance, std::move(current), std::move(slices), id, settings.ns,
                total - correction);
    }

    SpaceChargeTransportResult transportSpaceCharge(
            const SpaceChargeTransportPlan& plan, const SpaceChargeMatrix& input,
            opalx::spacecharge::analytic3d::Profile profile, double charge,
            const opalx::spacecharge::analytic3d::IntegrationSettings& quadrature) {
        require(!plan.slices().empty(), "Empty transport plan.");
        auto covariance = plan.entrance().convertCovariance(
                input, Chart::CommonTimeCanonical, Chart::CommonTimeCanonical);
        std::vector<SpaceChargeTransportStep> steps;
        steps.reserve(plan.slices().size());
        std::optional<tps::Map<double>> product;
        try {
            for (std::size_t i = 0; i < plan.slices().size(); ++i) {
                const auto& slice   = plan.slices()[i];
                const auto midpoint = transportSpaceChargeCovariance(slice.left, covariance);
                auto kick           = makeSpaceChargeKick(
                        slice.midpoint, midpoint, Chart::CommonTimeCanonical, profile, charge,
                        slice.integrationLength,
                        plan.id() + "/slice:" + std::to_string(i) + "/kick", quadrature);
                const auto afterKick = transportSpaceChargeCovariance(kick.polynomial, midpoint);
                covariance           = transportSpaceChargeCovariance(slice.right, afterKick);
                auto interval =
                        tps::compose(slice.right, tps::compose(kick.polynomial, slice.left), 1);
                product = product ? tps::compose(interval, *product, 1) : std::move(interval);
                steps.push_back({midpoint, covariance, std::move(kick)});
            }
        } catch (const std::invalid_argument& error) {
            throw OpalException("transportSpaceCharge", error.what());
        }
        return {covariance, std::move(*product), std::move(steps)};
    }
}  // namespace opalx::maps
