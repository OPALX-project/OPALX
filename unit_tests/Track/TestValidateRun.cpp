//
// Tests for validate(): run descriptions built by hand, no input file.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include "Physics/Physics.h"
#include "Track/ValidateRun.h"

using namespace opalx::run;

namespace {

    AttributeEntry entry(const std::string& name, AttributeValue value, bool set = true) {
        AttributeEntry e;
        e.name        = name;
        e.hasValue    = true;
        e.defaultUsed = !set;
        e.value       = std::move(value);
        return e;
    }

    ElementDescription element(
            const std::string& name, const std::string& type, const std::string& keyword) {
        ElementDescription e;
        e.object.name    = name;
        e.object.keyword = keyword;
        e.elementType    = type;
        return e;
    }

    /// A valid run: one proton beam with one source, a field solver and a LINE of a drift.
    RunDescription validRun() {
        RunDescription run;
        run.track.beamNames = {"B1"};
        run.track.dT        = {1e-11};
        run.track.maxSteps  = {100};
        run.track.zStop     = {1.0};
        run.track.run.attributes.push_back(entry("FIELDSOLVER", std::string("FS1")));
        run.track.run.attributes.push_back(entry("TURNS", 1.0, false));
        run.track.run.attributes.push_back(entry("SPECTRALTUNES", false, false));
        run.track.fieldSolverType = "NONE";
        run.fieldSolver           = ObjectDescription{"FIELDSOLVER", "FS1", {}};

        BeamDescription beam;
        beam.object.name             = "B1";
        beam.particleName            = "PROTON";
        beam.reference.mass          = Physics::m_p * 1e9;
        beam.reference.charge        = 1;
        beam.reference.beta          = 0.4;
        beam.emissionSourceList      = "SRC";
        beam.emissionSourceListFound = true;
        beam.hasEmissionSources      = true;
        SourceDescription source;
        source.name              = "ES1";
        source.distribution      = "DI";
        source.distributionFound = true;
        source.distributionType  = "GAUSS";
        beam.sources.push_back(source);
        run.beams.push_back(beam);

        run.line.name     = "L1";
        run.line.topology = "LINE";
        run.line.elements.push_back(element("D1", "Drift", "DRIFT"));
        return run;
    }

    bool contains(const std::vector<std::string>& errors, const std::string& text) {
        return std::any_of(errors.begin(), errors.end(), [&](const std::string& e) {
            return e.find(text) != std::string::npos;
        });
    }

    RunDescription ringRun() {
        RunDescription run = validRun();
        run.line.topology  = "RING";
        run.line.length    = 10.0;
        run.track.zStop    = {1e6};
        return run;
    }

}  // namespace

TEST(ValidateRun, ValidRunHasNoErrors) {
    EXPECT_TRUE(validate(validRun()).empty());
    EXPECT_TRUE(validate(ringRun()).empty());
}

TEST(ValidateRun, MissingFieldSolverAndBeam) {
    RunDescription run                   = validRun();
    run.track.run.attributes[0].hasValue = false;  // FIELDSOLVER, as freeze() records it unset
    run.track.beamNames.push_back("B2");
    const auto errors = validate(run);
    EXPECT_TRUE(contains(errors, "\"FIELDSOLVER\" must be set in \"RUN\" command."));
    EXPECT_TRUE(contains(errors, "Beam \"B2\" not found."));
}

TEST(ValidateRun, ReportsAllErrorsTogether) {
    RunDescription run    = validRun();
    run.beams[0].isPhoton = true;
    run.beams[0].sources.clear();
    run.beams[0].untilMaxSteps = true;
    run.track.zStopGiven       = true;
    const auto errors          = validate(run);
    ASSERT_EQ(errors.size(), 3u) << joinErrors(errors);
    EXPECT_TRUE(contains(errors, "PARTICLE=PHOTON"));
    EXPECT_TRUE(contains(errors, "must contain at least one EMISSIONSOURCE"));
    EXPECT_TRUE(contains(errors, "has both a ZSTOP (from TRACK) and UNTILMAXSTEPS"));
}

TEST(ValidateRun, RestartRefusesEmittingSourcesAndDecay) {
    RunDescription run              = validRun();
    run.context.restart             = true;
    run.beams[0].sources[0].emitted = true;
    run.beams[0].globalProcesses    = {"DECAY"};
    run.beams[0].daughterBeam       = "";
    const auto errors               = validate(run);
    EXPECT_TRUE(contains(errors, "time-dependent or delayed emission sources (beam 'B1'"));
    EXPECT_TRUE(contains(errors, "stochastic global processes such as DECAY (beam 'B1')"));
}

TEST(ValidateRun, RingSettings) {
    RunDescription run = ringRun();
    run.line.length    = 0;
    EXPECT_TRUE(contains(validate(run), "A RING requires a positive circumference."));

    run                            = validRun();
    run.track.kineticEnergyStopGeV = 0.1;
    const auto errors              = validate(run);
    EXPECT_TRUE(contains(errors, "EKINSTOP requires a non-restarted RING"));
    EXPECT_TRUE(contains(errors, "EKINSTOP requires the single-container cyclotron gap path."));

    run                         = ringRun();
    run.track.run.attributes[1] = entry("TURNS", 2.5);
    EXPECT_TRUE(contains(validate(run), "TURNS must be a positive integer for a RING."));
}

TEST(ValidateRun, TurnsElementTypes) {
    RunDescription run          = ringRun();
    run.track.run.attributes[1] = entry("TURNS", 3.0);
    run.line.elements.push_back(element("C1", "RFCavity", "RFCAVITY"));
    run.line.elements.push_back(element("V1", "Any", "VARIABLE_RF_CAVITY"));
    run.line.elements.push_back(element("E1", "Any", "ENGE"));
    run.line.elements.push_back(element("P1", "Probe", "PROBE"));
    const auto errors = validate(run);
    ASSERT_EQ(errors.size(), 1u) << joinErrors(errors);
    EXPECT_TRUE(contains(errors, "unsupported element C1"));
}

TEST(ValidateRun, ProbeNeedsRingOrCyclotron) {
    RunDescription run = validRun();
    run.line.elements.push_back(element("P1", "Probe", "PROBE"));
    EXPECT_TRUE(contains(validate(run), "PROBE requires cyclotron-sector or RING"));
    run.line.elements.push_back(element("S1", "CyclotronSector", "CYCLOTRONSECTOR"));
    EXPECT_TRUE(validate(run).empty());
    run.context.restart = true;
    EXPECT_TRUE(contains(validate(run), "PROBE does not yet support restart"));
}

TEST(ValidateRun, SingleGapNeedsOneUnpolarizedProton) {
    RunDescription run     = ringRun();
    ElementDescription gap = element("G1", "RFCavity", "RFCAVITY");
    gap.object.attributes.push_back(entry("TYPE", std::string("SINGLEGAP")));
    gap.object.attributes.push_back(entry("FMAPFN", std::string("gap.txt")));
    run.line.elements.push_back(gap);
    EXPECT_TRUE(validate(run).empty());
    run.beams[0].hasPolarization = true;
    EXPECT_TRUE(contains(validate(run), "SINGLEGAP requires an unpolarized positive proton"));
}

TEST(ValidateRun, DecayDaughterBeams) {
    RunDescription run           = validRun();
    BeamDescription mu           = run.beams[0];
    mu.object.name               = "MU";
    mu.hasExplicitEnergy         = true;
    run.beams[0].globalProcesses = {"DECAY"};
    run.beams[0].daughterBeam    = "MU";
    run.beams.push_back(mu);
    run.track.beamNames.push_back("MU");
    const auto errors = validate(run);
    EXPECT_TRUE(
            contains(errors, "\"GAMMA\", \"ENERGY\" and \"PC\" are not allowed on beam \"MU\""));
    EXPECT_TRUE(contains(errors, "\"SOURCES\" is not allowed on beam \"MU\""));
}
