//
// Unit tests for BeamBeam HDF5 diagnostics metadata normalization.
//

#include <gtest/gtest.h>

#include "Structure/H5BeamBeamDiagnosticsWriter.h"

#include <cmath>
#include <limits>

TEST(BeamBeamDiagnosticsWriterTest, PreservesPhysicalParticleMeanS) {
    H5BeamBeamDiagnosticsWriter::StepMetadata meta;
    meta.pathLengthS      = 0.0128910489;
    meta.particleMeanR[2] = 0.0133407366;
    meta.particleMeanS    = 0.0133407366;

    EXPECT_DOUBLE_EQ(H5BeamBeamDiagnosticsWriter::normalizeParticleMeanS(meta), meta.particleMeanS);
}

TEST(BeamBeamDiagnosticsWriterTest, CorrectsLegacyDoubleCountedParticleMeanS) {
    H5BeamBeamDiagnosticsWriter::StepMetadata meta;
    meta.pathLengthS      = 0.0128910489;
    meta.particleMeanR[2] = 0.0133407366;
    meta.particleMeanS    = meta.pathLengthS + meta.particleMeanR[2];

    EXPECT_DOUBLE_EQ(
            H5BeamBeamDiagnosticsWriter::normalizeParticleMeanS(meta), meta.particleMeanR[2]);
}

TEST(BeamBeamDiagnosticsWriterTest, PreservesNonFiniteDiagnosticsForCallerVisibility) {
    H5BeamBeamDiagnosticsWriter::StepMetadata meta;
    meta.pathLengthS      = 0.0128910489;
    meta.particleMeanR[2] = 0.0133407366;
    meta.particleMeanS    = std::numeric_limits<double>::quiet_NaN();

    EXPECT_TRUE(std::isnan(H5BeamBeamDiagnosticsWriter::normalizeParticleMeanS(meta)));
}

TEST(BeamBeamDiagnosticsWriterTest, BuildsOneFilePerDiagnosticQuantityName) {
    EXPECT_EQ(
            H5BeamBeamDiagnosticsWriter::makeOutputFileName(
                    "gamma_gamma_pairs", "RHO", "scalar", "beambeam_rho_pre"),
            "data/gamma_gamma_pairs-RHO_scalar-beambeam_rho_pre.h5");
    EXPECT_EQ(
            H5BeamBeamDiagnosticsWriter::makeOutputFileName(
                    "gamma_gamma_pairs", "EF", "vector", "beambeam_e"),
            "data/gamma_gamma_pairs-EF_vector-beambeam_e.h5");
}
