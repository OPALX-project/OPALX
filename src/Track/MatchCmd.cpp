#include "Track/MatchCmd.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Track/CofCmd.h"
#include "Utilities/ResultJson.h"
#include "Utility/Inform.h"

extern Inform* gmsg;
namespace {
    enum { ORBIT, MODEEMITTANCES_N, OUTPUT, SIZE };
    void require(bool valid, const std::string& message) {
        if (!valid) throw OpalException("MATCH", message);
    }
    std::string serialize(const std::string& name, const MatchCmd::Result& result) {
        const auto& source = result.source;
        const auto& match  = result.transverse;
        opalx::report::ResultJsonWriter json;
        json.object(nullptr, [&] {
            json.text("schema", "opalx.transverse_match");
            json.number("schema_version", 1);
            json.text("match_id", name);
            json.text("scope", "static_magnetic_ring_zero_current_transverse_4d");
            json.object("source", [&] {
                json.text("cof_id", source.command);
                json.text("cof_generation", std::to_string(source.generation));
                json.text("ring", source.launch.line);
                json.text("species", source.launch.species);
                json.number("mass_eV", source.launch.massEV);
                json.number("charge_e", source.launch.chargeE);
                json.number("time_s", source.launch.time);
                const auto& position = source.launch.position;
                const auto& momentum = source.launch.momentum;
                json.numbers(
                        "launch_position_m", std::array{position[0], position[1], position[2]});
                json.numbers(
                        "launch_momentum_p_mc", std::array{momentum[0], momentum[1], momentum[2]});
                json.number(
                        "reference_momentum_p_mc",
                        std::hypot(momentum[0], momentum[1], momentum[2]));
                json.numbers("section_coordinates", source.coordinates);
                json.numbers("closure_residual", source.residual);
                const auto origin = source.launch.section.getOrigin();
                json.numbers("section_origin_m", std::array{origin[0], origin[1], origin[2]});
                const auto axes = source.launch.section.getRotationMatrix();
                json.array("laboratory_to_section_rotation", [&] {
                    for (unsigned i = 0; i < 3; ++i)
                        json.numbers(nullptr, std::array{axes(i, 0), axes(i, 1), axes(i, 2)});
                });
                json.array("one_turn_matrix", [&] {
                    for (const auto& row : source.matrix)
                        json.numbers(nullptr, row);
                });
                json.numbers("coordinate_scales", source.scales);
                json.number("normal_magnetic_field_T", source.normalMagneticField);
                json.numbers(
                        "section_magnetic_field_T",
                        std::array{
                                source.transverseMagneticField[0],
                                source.transverseMagneticField[1], source.normalMagneticField});
                json.number("relative_energy_drift", source.relativeEnergyDrift);
            });
            json.text("coordinate_order", "x,px,y,py");
            json.text("coordinate_units", "m,p/(mc),m,p/(mc)");
            json.text("coordinate_convention", "mechanical_fixed_energy_section");
            json.text("covariance_units", "products_of_coordinate_units");
            json.text("emittance_convention", "normalized_rms_mode_eigenemittances_m_no_pi");
            json.text("mode_order", "ascending_positive_imaginary_fractional_phase");
            json.numbers("mode_emittances_n_m", match.emittances);
            json.numbers("mode_fractional_phases", match.tunes);
            json.array("covariance", [&] {
                for (const auto& row : match.covariance)
                    json.numbers(nullptr, row);
            });
            json.object("diagnostics", [&] {
                json.number("symplectic_residual", match.symplecticResidual);
                json.number("periodic_residual", match.periodicResidual);
                json.number("normalization_residual", match.normalizationResidual);
                const MatchedCovariance::Settings defaults;
                json.number("symplectic_tolerance", defaults.symplecticTolerance);
                json.number("periodic_tolerance", defaults.periodicTolerance);
                json.number("minimum_mode_separation", defaults.minimumModeSeparation);
            });
        });
        return json.str();
    }
    void report(const std::string& name, const MatchCmd::Result& result) {
        const auto& match = result.transverse;
        std::ostringstream out;
        out << std::setprecision(17) << "MATCH " << name << ": COF=" << result.source.command
            << " generation=" << result.source.generation << " ring=" << result.source.launch.line
            << "\nScope: static magnetic RING, zero-current transverse 4D; no longitudinal "
               "equilibrium."
            << "\nCoordinates: (x,px,y,py), metres and mechanical p/(mc), fixed-energy section."
            << "\nNormalized RMS mode emittances [m], no pi; ascending positive-imaginary phase:\n";
        for (unsigned k = 0; k < 2; ++k)
            out << "  mode " << k + 1 << " phase=" << match.tunes[k]
                << " emittance=" << match.emittances[k] << '\n';
        out << "Covariance (entries have products of coordinate units):\n";
        MatchedCovariance::writeCovariance(out, match.covariance);
        out << "Symplectic residual=" << match.symplecticResidual
            << " periodic residual=" << match.periodicResidual
            << " normalization residual=" << match.normalizationResidual
            << "\nResiduals are numerical diagnostics; refine COF DT and FDSTEP to establish "
               "accuracy.\n";
        *gmsg << level4 << out.str() << endl;
    }
}  // namespace

MatchCmd::MatchCmd()
    : Action(SIZE, "MATCH", "Match zero-current transverse covariance to a named RING COF.") {
    itsAttr[ORBIT] = Attributes::makeString("ORBIT", "Required successful named COF snapshot.");
    itsAttr[MODEEMITTANCES_N] = Attributes::makeRealArray(
            "MODEEMITTANCES_N", "Two positive normalized RMS mode eigenemittances [m], no pi.");
    itsAttr[OUTPUT] = Attributes::makeString(
            "OUTPUT", "Optional exact JSON filename; atomic replacement on success.");
    registerOwnership(AttributeHandler::COMMAND);
}
MatchCmd::MatchCmd(const std::string& name, MatchCmd* parent) : Action(name, parent) {}
MatchCmd* MatchCmd::clone(const std::string& name) { return new MatchCmd(name, this); }

void MatchCmd::execute() {
    result_m.reset();
    require(!OpalData::getInstance()->inRestartRun(), "MATCH does not support restart mode.");
    Result candidate;
    candidate.source  = CofCmd::findOptics(Attributes::getString(itsAttr[ORBIT]));
    const auto values = Attributes::getRealArray(itsAttr[MODEEMITTANCES_N]);
    require(values.size() == 2, "MODEEMITTANCES_N requires exactly two positive values [m].");
    const std::array<double, 2> eps{values[0], values[1]};
    const auto output = Attributes::getString(itsAttr[OUTPUT]);
    require(itsAttr[OUTPUT].defaultUsed() || !output.empty(), "Explicit OUTPUT must not be empty.");
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string error;
    if (rank == 0) try {
            candidate.transverse = MatchedCovariance::compute(candidate.source, eps);
            candidate.json       = serialize(getOpalName(), candidate);
            if (!output.empty()) opalx::report::writeResult("MATCH", output, candidate.json);
        } catch (const std::exception& ex) {
            error = ex.what();
        } catch (const OpalException& ex) {
            error = ex.what();
        }
    const auto broadcast = [](std::string& text) {
        int n = static_cast<int>(text.size());
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        text.resize(n);
        MPI_Bcast(text.data(), n, MPI_CHAR, 0, MPI_COMM_WORLD);
    };
    broadcast(error);
    require(error.empty(), error);
    broadcast(candidate.json);
    // Numeric results are broadcast as values, never recovered by parsing JSON.
    std::array<double, 23> numbers{};
    if (rank == 0) {
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                numbers[4 * i + j] = candidate.transverse.covariance[i][j];
        for (unsigned k = 0; k < 2; ++k) {
            numbers[16 + k] = eps[k];
            numbers[18 + k] = candidate.transverse.tunes[k];
        }
        numbers[20] = candidate.transverse.symplecticResidual;
        numbers[21] = candidate.transverse.periodicResidual;
        numbers[22] = candidate.transverse.normalizationResidual;
    }
    MPI_Bcast(numbers.data(), numbers.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            candidate.transverse.covariance[i][j] = numbers[4 * i + j];
    for (unsigned k = 0; k < 2; ++k) {
        candidate.transverse.emittances[k] = numbers[16 + k];
        candidate.transverse.tunes[k]      = numbers[18 + k];
    }
    candidate.transverse.symplecticResidual    = numbers[20];
    candidate.transverse.periodicResidual      = numbers[21];
    candidate.transverse.normalizationResidual = numbers[22];
    if (rank == 0) report(getOpalName(), candidate);
    result_m = std::move(candidate);
}
const MatchCmd::Result& MatchCmd::findResult(const std::string& name) {
    auto* command = dynamic_cast<MatchCmd*>(OpalData::getInstance()->find(name));
    require(command && command->result_m.has_value(),
            "Result must name a successfully executed MATCH: " + name);
    return *command->result_m;
}
