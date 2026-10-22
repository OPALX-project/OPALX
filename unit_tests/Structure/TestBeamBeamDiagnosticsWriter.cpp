//
// Unit tests for BeamBeam HDF5 diagnostics metadata and collective I/O.
//

#include <gtest/gtest.h>

#include "PartBunch/CartesianDomain.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "Structure/H5BeamBeamDiagnosticsWriter.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <type_traits>

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

namespace {
    class BeamBeamDiagnosticsIOTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            char pattern[64] = "beambeam-diagnostics-XXXXXX";
            if (ippl::Comm->rank() == 0 && !mkdtemp(pattern)) {
                pattern[0] = '\0';
            }
            MPI_Bcast(pattern, sizeof(pattern), MPI_CHAR, 0, ippl::Comm->getCommunicator());
            ASSERT_NE(pattern[0], '\0');
            directory = pattern;
        }
        void TearDown() override {
            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0 && !directory.empty()) {
                std::filesystem::remove_all(directory);
            }
        }
        opalx::spacecharge::CartesianDomainConfig3D config() const {
            opalx::spacecharge::CartesianDomainConfig3D result;
            result.meshSize      = {4, 4, 4};
            result.decomposition = {true, true, true};
            return result;
        }
        std::filesystem::path directory;
    };

    TEST_F(BeamBeamDiagnosticsIOTest, ScalarOutputAppendsAndPreservesOwnedValues) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        auto& rho = fields.chargeDensity();
        rho       = 3.5;
        H5BeamBeamDiagnosticsWriter::StepMetadata meta;
        meta.shape      = {4, 4, 4};
        meta.spacing    = {1.0, 1.0, 1.0};
        const auto name = (directory / "rho.h5").string();
        H5BeamBeamDiagnosticsWriter::writeScalarQuantity(name, "rho", meta, rho);
        rho = -2.25;
        H5BeamBeamDiagnosticsWriter::writeScalarQuantity(name, "rho", meta, rho);

        auto props    = H5CreateFileProp();
        MPI_Comm comm = ippl::Comm->getCommunicator();
        EXPECT_EQ(H5SetPropFileMPIOCollective(props, &comm), H5_SUCCESS);
        auto file = H5OpenFile(name.c_str(), H5_O_RDONLY, props);
        EXPECT_EQ(H5CloseProp(props), H5_SUCCESS);
        EXPECT_EQ(H5GetNumSteps(file), 2);
        const auto owned = rho.getLayout().getLocalNDIndex();
        std::vector<h5_float64_t> values(owned[0].length() * owned[1].length() * owned[2].length());
        for (int step = 0; step < 2; ++step) {
            EXPECT_EQ(H5SetStep(file, step), H5_SUCCESS);
            EXPECT_EQ(
                    H5Block3dSetView(
                            file, owned[0].first(), owned[0].last(), owned[1].first(),
                            owned[1].last(), owned[2].first(), owned[2].last()),
                    H5_SUCCESS);
            EXPECT_EQ(H5Block3dReadScalarFieldFloat64(file, "rho", values.data()), H5_SUCCESS);
            for (const auto value : values) {
                EXPECT_DOUBLE_EQ(value, step == 0 ? 3.5 : -2.25);
            }
        }
        EXPECT_EQ(H5CloseFile(file), H5_SUCCESS);
    }

    TEST_F(BeamBeamDiagnosticsIOTest, PreparationFailureIsReportedOnEveryRankAndRetryWorks) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        fields.chargeDensity() = 1.0;
        const auto parent      = directory / "not-a-directory";
        if (ippl::Comm->rank() == 0) {
            std::ofstream file(parent);
            file << "sentinel";
        }
        ippl::Comm->barrier();
        const auto name = (parent / "rho.h5").string();
        EXPECT_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "rho", {}, fields.chargeDensity()),
                OpalException);
        EXPECT_THROW(H5BeamBeamDiagnosticsWriter writer(name), OpalException);
        if (ippl::Comm->rank() == 0) {
            std::filesystem::remove(parent);
        }
        ippl::Comm->barrier();
        EXPECT_NO_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "rho", {}, fields.chargeDensity()));
    }

    TEST_F(BeamBeamDiagnosticsIOTest, WriteFailureClosesHandleAndAllowsRetry) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        fields.chargeDensity() = 1.0;
        const auto name        = (directory / "retry.h5").string();
        // The missing intermediate group makes HDF5 reject the field after open.
        EXPECT_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "missing/rho", {}, fields.chargeDensity()),
                OpalException);
        EXPECT_NO_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "rho", {}, fields.chargeDensity()));
    }

    TEST_F(BeamBeamDiagnosticsIOTest, FailedHdf5OpenIsCheckedWithoutAssertions) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        fields.chargeDensity() = 1.0;
        const auto name        = (directory / "corrupted.h5").string();
        H5BeamBeamDiagnosticsWriter::writeScalarQuantity(name, "rho", {}, fields.chargeDensity());
        const auto backup = directory / "valid.h5";
        // A valid signature alone is not sufficient: also reject a truncated header.
        for (const std::string& damaged :
             {std::string("not an HDF5 file"), std::string("\211HDF\r\n\032\n", 8)}) {
            if (ippl::Comm->rank() == 0) {
                std::filesystem::copy_file(
                        name, backup, std::filesystem::copy_options::overwrite_existing);
                std::ofstream corrupt(name, std::ios::binary);
                corrupt << damaged;
            }
            ippl::Comm->barrier();
            H5E_auto2_t handler = nullptr;
            void* handlerData   = nullptr;
            EXPECT_GE(H5Eget_auto2(H5E_DEFAULT, &handler, &handlerData), 0);
#if H5_VERSION_GE(1, 14, 5)
            bool pausedBefore = false;
            EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &pausedBefore), 0);
#endif
            const auto openFiles = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_FILE);
            EXPECT_THROW(
                    H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                            name, "rho", {}, fields.chargeDensity()),
                    OpalException);
#if H5_VERSION_GE(1, 14, 5)
            // The GH200 workaround must suppress recording, not merely printing.
            // Check before another HDF5 API call can clear the error stack.
            EXPECT_EQ(H5Eget_num(H5E_DEFAULT), 0);
#endif
            EXPECT_THROW(
                    H5BeamBeamDiagnosticsWriter::writeVectorQuantity(
                            name, "E", {}, fields.electricField()),
                    OpalException);
#if H5_VERSION_GE(1, 14, 5)
            EXPECT_EQ(H5Eget_num(H5E_DEFAULT), 0);
#endif
            H5E_auto2_t restoredHandler = nullptr;
            void* restoredHandlerData   = nullptr;
            EXPECT_GE(H5Eget_auto2(H5E_DEFAULT, &restoredHandler, &restoredHandlerData), 0);
            EXPECT_EQ(restoredHandler, handler);
            EXPECT_EQ(restoredHandlerData, handlerData);
            EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_FILE), openFiles);
#if H5_VERSION_GE(1, 14, 5)
            bool pausedAfter = false;
            EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &pausedAfter), 0);
            EXPECT_EQ(pausedAfter, pausedBefore);
#endif
            if (ippl::Comm->rank() == 0) {
                std::ifstream rejected(name, std::ios::binary);
                const std::string contents{
                        std::istreambuf_iterator<char>(rejected), std::istreambuf_iterator<char>()};
                EXPECT_EQ(contents, damaged);  // No silent truncation/replacement on failure.
                rejected.close();
                std::filesystem::copy_file(
                        backup, name, std::filesystem::copy_options::overwrite_existing);
            }
            ippl::Comm->barrier();
            EXPECT_NO_THROW(
                    H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                            name, "rho", {}, fields.chargeDensity()));
        }
    }

// HDF5 error-stack pause/resume was introduced in 1.14.5. Older installations
// still exercise the corrupt-file rejection and callback restoration above.
#if H5_VERSION_GE(1, 14, 5)
    TEST_F(BeamBeamDiagnosticsIOTest, AppendProbePreservesNestedErrorStackPause) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        fields.chargeDensity() = 1.0;
        const auto name        = (directory / "paused.h5").string();
        H5BeamBeamDiagnosticsWriter::writeScalarQuantity(name, "rho", {}, fields.chargeDensity());

        bool paused = true;
        // Do not return early on just one rank before the collective writes.
        EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &paused), 0);
        EXPECT_FALSE(paused);
        EXPECT_GE(H5Epause_stack(H5E_DEFAULT), 0);
        // A successful append must not resume the caller's surrounding pause.
        EXPECT_NO_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "rho", {}, fields.chargeDensity()));
        EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &paused), 0);
        EXPECT_TRUE(paused);

        if (ippl::Comm->rank() == 0) {
            std::ofstream corrupt(name, std::ios::binary);
            corrupt << "not an HDF5 file";
        }
        ippl::Comm->barrier();
        EXPECT_THROW(
                H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
                        name, "rho", {}, fields.chargeDensity()),
                OpalException);
        EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &paused), 0);
        EXPECT_TRUE(paused);
        EXPECT_GE(H5Eresume_stack(H5E_DEFAULT), 0);
        EXPECT_GE(H5Eis_paused(H5E_DEFAULT, &paused), 0);
        EXPECT_FALSE(paused);  // No leaked extra pause after success or rejection.
    }
#endif

    TEST_F(BeamBeamDiagnosticsIOTest, VectorAndTwoStageOutputsPreserveValues) {
        using namespace opalx::spacecharge;
        CartesianDomain<double, 3> domain(config());
        CartesianPIC3DFieldStorage<double, 3> fields(domain);
        fields.initializeFields(PoissonSolverType::None);
        fields.chargeDensity() = 3.5;
        fields.electricField() = Vector_t<double, 3>(1.0, -2.0, 4.0);
        H5BeamBeamDiagnosticsWriter::StepMetadata meta;
        meta.shape            = {4, 4, 4};
        meta.spacing          = {1.0, 1.0, 1.0};
        const auto vectorName = (directory / "vector.h5").string();
        const auto stepName   = (directory / "step.h5").string();
        H5BeamBeamDiagnosticsWriter::writeVectorQuantity(
                vectorName, "E", meta, fields.electricField());
        {
            H5BeamBeamDiagnosticsWriter writer(stepName);
            EXPECT_THROW(
                    writer.endStep(fields.chargeDensity(), fields.electricField()), OpalException);
            writer.beginStep(meta, fields.chargeDensity());
            fields.chargeDensity() = -0.5;  // Cached rho must stay 3.5; use -0.5 as phi.
            writer.endStep(fields.chargeDensity(), fields.electricField());
            writer.close();
        }
        const auto owned = fields.chargeDensity().getLayout().getLocalNDIndex();
        const auto size  = owned[0].length() * owned[1].length() * owned[2].length();
        for (const bool vector : {true, false}) {
            const auto& name = vector ? vectorName : stepName;
            auto props       = H5CreateFileProp();
            MPI_Comm comm    = ippl::Comm->getCommunicator();
            EXPECT_EQ(H5SetPropFileMPIOCollective(props, &comm), H5_SUCCESS);
            auto file = H5OpenFile(name.c_str(), H5_O_RDONLY, props);
            EXPECT_EQ(H5CloseProp(props), H5_SUCCESS);
            EXPECT_EQ(H5GetNumSteps(file), 1);
            EXPECT_EQ(H5SetStep(file, 0), H5_SUCCESS);
            EXPECT_EQ(
                    H5Block3dSetView(
                            file, owned[0].first(), owned[0].last(), owned[1].first(),
                            owned[1].last(), owned[2].first(), owned[2].last()),
                    H5_SUCCESS);
            std::vector<h5_float64_t> x(size), y(size), z(size);
            if (vector) {
                EXPECT_EQ(
                        H5Block3dReadVector3dFieldFloat64(file, "E", x.data(), y.data(), z.data()),
                        H5_SUCCESS);
                for (std::size_t i = 0; i < x.size(); ++i) {
                    EXPECT_DOUBLE_EQ(x[i], 1.0);
                    EXPECT_DOUBLE_EQ(y[i], -2.0);
                    EXPECT_DOUBLE_EQ(z[i], 4.0);
                }
            } else {
                for (const auto& [field, value] : std::vector<std::pair<const char*, double>>{
                             {"rho", 3.5}, {"phi", -0.5}, {"Ex", 1.0}, {"Ey", -2.0}, {"Ez", 4.0}}) {
                    EXPECT_EQ(H5Block3dReadScalarFieldFloat64(file, field, x.data()), H5_SUCCESS);
                    for (const auto actual : x)
                        EXPECT_DOUBLE_EQ(actual, value);
                }
            }
            EXPECT_EQ(H5CloseFile(file), H5_SUCCESS);
        }
    }

    TEST_F(BeamBeamDiagnosticsIOTest, ExplicitCloseIsIdempotentAndDestructorIsNonthrowing) {
        static_assert(std::is_nothrow_destructible_v<H5BeamBeamDiagnosticsWriter>);
        const auto name = (directory / "instance.h5").string();
        {
            H5BeamBeamDiagnosticsWriter writer(name);
            EXPECT_NO_THROW(writer.close());
            EXPECT_NO_THROW(writer.close());
        }
        // Reopen after destructor-owned cleanup too; handles must not remain live.
        {
            H5BeamBeamDiagnosticsWriter writer(name);
        }
        {
            H5BeamBeamDiagnosticsWriter writer(name);
        }
    }
}  // namespace
