#include <gtest/gtest.h>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "SpaceCharge/SpaceChargeDiagnostics.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

namespace opalx::spacecharge {
    namespace {

        class FFT2D5AlgorithmTest : public ::testing::Test {
        protected:
            using ParticleContainer = FFT2D5Algorithm::ParticleContainer;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                pathFile_m = std::filesystem::current_path() / "fft2d5-design-path.dat";
                std::ofstream path(pathFile_m);
                path << "s rx ry rz px py pz ex ey ez bx by bz ke t\n";
                path << "0 0 0 0 0 0 1 0 0 0 0 0 0 0 0\n";
                path << "1 0 0 1 0 0 1 0 0 0 0 0 0 0 0\n";
                path.close();

                ippl::NDIndex<3> domain;
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    domain[dimension] = ippl::Index(4);
                }
                mesh_m = std::make_unique<Mesh_t<3>>(
                        domain, Vector_t<double, 3>(0.25), Vector_t<double, 3>(-0.5));
                layout_m = std::make_unique<FieldLayout_t<3>>(
                        MPI_COMM_WORLD, domain, std::array<bool, 3>{false, false, false}, false);
                particles_m = std::make_shared<ParticleContainer>(*mesh_m, *layout_m);
                particles_m->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            }

            void TearDown() override {
                particles_m.reset();
                layout_m.reset();
                mesh_m.reset();
                std::filesystem::remove(pathFile_m);
            }

            ParticleFieldBinding3D binding() { return makeParticleFieldBinding(*particles_m); }

            FFT2D5Config config(FFT2D5LongitudinalFieldMode mode) const {
                FFT2D5Config::Parameters values;
                values.meshSize              = {4, 4, 4};
                values.longitudinalFieldMode = mode;
                values.pipeSizeX             = 1.0;
                values.pipeSizeY             = 1.0;
                values.beamRadius            = 0.1;
                values.referencePathFile     = pathFile_m.string();
                return FFT2D5Config(std::move(values));
            }

            std::filesystem::path pathFile_m;
            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::shared_ptr<ParticleContainer> particles_m;
            CoordinateSystemTrafo trackerToSolve_m;
            CoordinateSystemTrafo solveToTracker_m;
        };

        TEST_F(FFT2D5AlgorithmTest, LazyInitializationBuildsPersistentSliceSolvers) {
            std::array bindings{binding()};
            FFT2D5Algorithm solver(config(FFT2D5LongitudinalFieldMode::Open), bindings);
            EXPECT_FALSE(solver.initialized());

            ParticleFieldSet particleView(bindings, 0);
            CoordinateFrameTransforms frames{trackerToSolve_m, solveToTracker_m};
            SpaceChargeStepState step{
                    0, 0.0, 1.0e-12, false, 1.0, ippl::Comm->size(), std::move(frames)};
            SpaceChargeSolveContext solveContext(std::move(particleView), std::move(step));
            SpaceChargeDiagnostics diagnostics;

            EXPECT_NO_THROW(solver.solve(solveContext, diagnostics));
            EXPECT_TRUE(solver.initialized());
            EXPECT_EQ(solver.sliceCount(), 4u);
            EXPECT_EQ(solver.referencePathView().extent(0), 2u);
            EXPECT_EQ(diagnostics.backendSolveCount(), 4u);
        }

        TEST_F(FFT2D5AlgorithmTest, EveryMasterLongitudinalModeIsAccepted) {
            for (const FFT2D5LongitudinalFieldMode mode :
                 {FFT2D5LongitudinalFieldMode::Open, FFT2D5LongitudinalFieldMode::Cylindrical,
                  FFT2D5LongitudinalFieldMode::Plates, FFT2D5LongitudinalFieldMode::None}) {
                EXPECT_NO_THROW(config(mode));
            }
        }

        TEST_F(FFT2D5AlgorithmTest, SelectionModesUsePrimaryOrEveryTrackingActiveBinding) {
            std::array bindings{binding(), binding(), binding()};
            bindings[0].trackingActive = true;
            bindings[1].trackingActive = true;
            bindings[2].trackingActive = false;
            ParticleFieldSet particles(bindings, 0);

            particles.selectForSolve(ParticleSelectionMode::AllTrackingActive);
            EXPECT_TRUE(bindings[0].selectedForSolve);
            EXPECT_TRUE(bindings[1].selectedForSolve);
            EXPECT_FALSE(bindings[2].selectedForSolve);

            particles.selectForSolve(ParticleSelectionMode::PrimaryOnly);
            EXPECT_TRUE(bindings[0].selectedForSolve);
            EXPECT_FALSE(bindings[1].selectedForSolve);
            EXPECT_FALSE(bindings[2].selectedForSolve);
        }

    }  // namespace
}  // namespace opalx::spacecharge
