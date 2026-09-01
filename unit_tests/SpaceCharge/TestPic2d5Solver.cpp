#include <gtest/gtest.h>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/Pic2d5/Pic2d5Solver.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

namespace opalx::spacecharge {
    namespace {

        class Pic2d5SolverTest : public ::testing::Test {
        protected:
            using ParticleContainer = Pic2d5Solver::ParticleContainer;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void SetUp() override {
                pathFile_m = std::filesystem::current_path() / "pic2d5-design-path.dat";
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

            ParticleContainerAttributes attributes() {
                ParticleContainerAttributes result;
                result.position = ParticleAttributeHandle::writable(
                        ParticleAttribute::Position, particles_m->R);
                result.momentum = ParticleAttributeHandle::readable(
                        ParticleAttribute::Momentum, particles_m->P);
                result.charge =
                        ParticleAttributeHandle::readable(ParticleAttribute::Charge, *particles_m);
                result.mass =
                        ParticleAttributeHandle::readable(ParticleAttribute::Mass, *particles_m);
                result.timeStep = ParticleAttributeHandle::writable(
                        ParticleAttribute::TimeStep, particles_m->dt);
                result.electricField = ParticleAttributeHandle::writable(
                        ParticleAttribute::ElectricField, particles_m->E);
                result.magneticField = ParticleAttributeHandle::writable(
                        ParticleAttribute::MagneticField, particles_m->B);
                result.invalidMask = ParticleAttributeHandle::readable(
                        ParticleAttribute::InvalidMask, particles_m->InvalidMask);
                result.bin =
                        ParticleAttributeHandle::writable(ParticleAttribute::Bin, particles_m->Bin);
                return result;
            }

            Pic2d5Config config(Pic2d5LongitudinalFieldMode mode) const {
                Pic2d5ConfigValues values;
                values.meshSize              = {4, 4, 4};
                values.longitudinalFieldMode = mode;
                values.pipeSizeX             = 1.0;
                values.pipeSizeY             = 1.0;
                values.beamRadius            = 0.1;
                values.referencePathFile     = pathFile_m.string();
                return Pic2d5Config(std::move(values));
            }

            std::filesystem::path pathFile_m;
            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::shared_ptr<ParticleContainer> particles_m;
            CoordinateSystemTrafo trackerToSolve_m;
            CoordinateSystemTrafo solveToTracker_m;
        };

        TEST_F(Pic2d5SolverTest, LazyInitializationBuildsPersistentSliceBackends) {
            std::array containers{particles_m};
            Pic2d5Solver solver(config(Pic2d5LongitudinalFieldMode::Open), containers);
            EXPECT_FALSE(solver.initialized());

            std::array views{ParticleContainerView("primary", attributes(), true, true)};
            ParticleSetView particleView(views, 0, 0);
            FrameState frames;
            frames.trackerToSolve = BorrowedHostObject::reference(trackerToSolve_m);
            frames.solveToTracker = BorrowedHostObject::reference(solveToTracker_m);
            CommunicatorView communicator{
                    BorrowedHostObject::reference(*ippl::Comm), ippl::Comm->rank(),
                    ippl::Comm->size()};
            StepState step{
                    0, 0.0, 1.0e-12, 0, false, std::move(communicator), {}, std::move(frames), 1.0};
            SolveContext solveContext(std::move(particleView), std::move(step));
            SelfFieldDiagnostics diagnostics;

            EXPECT_NO_THROW(solver.execute(solveContext, diagnostics));
            EXPECT_TRUE(solver.initialized());
            EXPECT_EQ(solver.sliceCount(), 4u);
            EXPECT_EQ(solver.referencePathView().extent(0), 2u);
            EXPECT_EQ(diagnostics.backendSolveCount(), 4u);
        }

        TEST_F(Pic2d5SolverTest, EveryMasterLongitudinalModeIsAccepted) {
            for (const Pic2d5LongitudinalFieldMode mode :
                 {Pic2d5LongitudinalFieldMode::Open, Pic2d5LongitudinalFieldMode::Cylindrical,
                  Pic2d5LongitudinalFieldMode::Plates, Pic2d5LongitudinalFieldMode::None}) {
                EXPECT_NO_THROW(config(mode));
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
