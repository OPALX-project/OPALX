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
        run.track.run.attributes.push_back(entry("METHOD", std::string("PARALLEL")));
        run.track.fieldSolverType = "NONE";
        run.fieldSolver           = ObjectDescription{"FIELDSOLVER", "FS1", {}};

        BeamDescription beam;
        beam.object.name  = "B1";
        beam.particleName = "PROTON";
        beam.reference.setM(Physics::m_p * 1e9);
        beam.reference.setQ(1);
        beam.reference.setBeta(0.4);
        beam.emissionSourceList      = "SRC";
        beam.emissionSourceListFound = true;
        beam.hasEmissionSources      = true;
        beam.hasExplicitEnergy       = true;
        SourceDescription source;
        source.name              = "ES1";
        source.distribution      = "DI";
        source.distributionFound = true;
        source.distributionType  = "GAUSS";
        beam.sources.push_back(source);
        run.beams.push_back(beam);

        run.line.name = "L1";
        run.line.elements.push_back(element("D1", "Drift", "DRIFT"));
        return run;
    }

    bool contains(const std::vector<std::string>& errors, const std::string& text) {
        return std::any_of(errors.begin(), errors.end(), [&](const std::string& e) {
            return e.find(text) != std::string::npos;
        });
    }

}  // namespace

TEST(ValidateRun, ValidRunHasNoErrors) { EXPECT_TRUE(validate(validRun()).empty()); }

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

TEST(ValidateRun, MissingDistributionOnFreshRunAndRestart) {
    RunDescription run                        = validRun();
    run.beams[0].sources[0].distributionFound = false;
    EXPECT_TRUE(contains(validate(run), "Distribution \"DI\" not found"));
    run.context.restart = true;
    EXPECT_TRUE(contains(validate(run), "Distribution \"DI\" not found"));
}

TEST(ValidateRun, GeneratedDistributionNeedsBeamEnergy) {
    RunDescription run             = validRun();
    run.beams[0].hasExplicitEnergy = false;
    EXPECT_TRUE(contains(validate(run), "The energy of beam \"B1\" hasn't been set."));
    run.beams[0].sources[0].distributionType = "FROMFILE";
    EXPECT_TRUE(validate(run).empty());
}

TEST(ValidateRun, ReferenceDirectionNeedsFileDistribution) {
    RunDescription run                 = validRun();
    run.beams[0].hasReferenceDirection = true;
    EXPECT_TRUE(contains(validate(run), "\"REFDIRECTION\" on beam \"B1\" needs a FROMFILE"));
    run.beams[0].sources[0].distributionType = "EMITTEDFROMFILE";
    EXPECT_TRUE(validate(run).empty());
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
