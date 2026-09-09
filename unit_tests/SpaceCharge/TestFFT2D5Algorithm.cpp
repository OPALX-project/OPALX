#include <gtest/gtest.h>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"

#include <array>
#include <cmath>
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
                particles_m  = std::make_shared<ParticleContainer>(*mesh_m, *layout_m);
                bunchState_m = std::make_shared<BunchStateHandler>();
                particles_m->setBunchStateHandler(bunchState_m);
            }

            void TearDown() override {
                particles_m.reset();
                layout_m.reset();
                mesh_m.reset();
                std::filesystem::remove(pathFile_m);
            }

            FFT2D5Config config(FFT2D5LongitudinalFieldMode mode) const {
                FFT2D5Config values;
                values.grid.meshSize         = {4, 4, 4};
                values.longitudinalFieldMode = mode;
                values.pipeSizeX             = 1.0;
                values.pipeSizeY             = 1.0;
                values.beamRadius            = 0.1;
                values.referencePathFile     = pathFile_m.string();
                return values;
            }

            std::filesystem::path pathFile_m;
            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::shared_ptr<ParticleContainer> particles_m;
            std::shared_ptr<BunchStateHandler> bunchState_m;
            CoordinateSystemTrafo trackerToSolve_m;
            CoordinateSystemTrafo solveToTracker_m;
        };

        TEST_F(FFT2D5AlgorithmTest, LazyInitializationBuildsPersistentSliceSolvers) {
            std::array particles{particles_m.get()};
            FFT2D5Algorithm solver(
                    config(FFT2D5LongitudinalFieldMode::Open), particles, bunchState_m);
            EXPECT_FALSE(solver.initialized());

            CoordinateFrameTransforms frames{trackerToSolve_m, solveToTracker_m};
            SpaceChargeStepState step{
                    0, 0.0, 1.0e-12, false, 1.0, ippl::Comm->size(), std::move(frames)};
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeSolveContext solveContext(activity, std::move(step));

            const SpaceChargeSolveResult result = solver.solve(solveContext);
            EXPECT_TRUE(solver.initialized());
            EXPECT_EQ(solver.sliceCount(), 4u);
            EXPECT_EQ(solver.referencePathView().extent(0), 2u);
            EXPECT_EQ(result.backendSolves, 4u);
        }

        TEST_F(FFT2D5AlgorithmTest, RejectsFixedCartesianDomain) {
            std::array particles{particles_m.get()};
            FFT2D5Algorithm solver(
                    config(FFT2D5LongitudinalFieldMode::Open), particles, bunchState_m);
            bunchState_m->setFixedCartesianDomain({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeSolveContext context(activity, {});
            EXPECT_THROW(static_cast<void>(solver.solve(context)), OpalException);
        }

        TEST_F(FFT2D5AlgorithmTest, RestoresTrackerFrameMoments) {
            particles_m->createParticles(2);
            particles_m->setQ(1.0e-15);
            particles_m->setM(Physics::m_e);
            auto positions = particles_m->R.getHostMirror();
            auto momenta   = particles_m->P.getHostMirror();
            auto timeSteps = particles_m->dt.getHostMirror();
            positions(0)   = Vector_t<double, 3>(0.0, 0.0, 0.25);
            positions(1)   = Vector_t<double, 3>(0.0, 0.0, 0.75);
            momenta(0)     = Vector_t<double, 3>(0.0, 0.0, 0.1);
            momenta(1)     = Vector_t<double, 3>(0.0, 0.0, 0.1);
            timeSteps(0)   = 1.0e-12;
            timeSteps(1)   = 1.0e-12;
            Kokkos::deep_copy(particles_m->R.getView(), positions);
            Kokkos::deep_copy(particles_m->P.getView(), momenta);
            Kokkos::deep_copy(particles_m->dt.getView(), timeSteps);

            std::array particles{particles_m.get()};
            FFT2D5Algorithm solver(
                    config(FFT2D5LongitudinalFieldMode::Open), particles, bunchState_m);
            const CoordinateSystemTrafo trackerToSolve(
                    Vector_t<double, 3>(0.0, 0.0, 0.1), Quaternion(1.0, 0.0, 0.0, 0.0));
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            step.frames   = {trackerToSolve, trackerToSolve.inverted()};
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeSolveContext context(activity, step);

            const SpaceChargeSolveResult result = solver.solve(context);
            EXPECT_EQ(result.backendSolves, 4u);
            EXPECT_FALSE(particles_m->isMomentsDirty());
            const auto restored = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->R.getView());
            for (std::size_t index = 0; index < 2; ++index) {
                for (unsigned dimension = 0; dimension < 3; ++dimension) {
                    EXPECT_NEAR(restored(index)[dimension], positions(index)[dimension], 1.0e-14);
                }
            }
            EXPECT_NEAR(particles_m->getMeanR()[2], 0.5, 1.0e-14);
        }

        TEST_F(FFT2D5AlgorithmTest, RepeatedSolveReplacesFieldsAndPreservesInactiveContainer) {
            particles_m->createParticles(2);
            particles_m->setQ(1.0e-15);
            particles_m->setM(Physics::m_e);
            auto positions = particles_m->R.getHostMirror();
            positions(0)   = Vector_t<double, 3>(-0.1, 0.0, 0.25);
            positions(1)   = Vector_t<double, 3>(0.1, 0.0, 0.75);
            Kokkos::deep_copy(particles_m->R.getView(), positions);
            Kokkos::deep_copy(particles_m->P.getView(), Vector_t<double, 3>(0.0, 0.0, 0.1));
            Kokkos::deep_copy(particles_m->dt.getView(), 1.0e-12);
            ParticleContainer secondary(*mesh_m, *layout_m);
            secondary.setBunchStateHandler(bunchState_m);
            secondary.createParticles(1);
            Kokkos::deep_copy(secondary.R.getView(), Vector_t<double, 3>(100.0));
            Kokkos::deep_copy(secondary.E.getView(), Vector_t<double, 3>(7.0));
            Kokkos::deep_copy(secondary.B.getView(), Vector_t<double, 3>(9.0));
            std::array particles{particles_m.get(), &secondary};
            FFT2D5Algorithm solver(
                    config(FFT2D5LongitudinalFieldMode::Open), particles, bunchState_m);
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            const std::array<std::uint8_t, 2> activity{1, 0};
            SpaceChargeSolveContext context(activity, step);
            EXPECT_EQ(solver.solve(context).backendSolves, 4u);
            auto firstE = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->E.getView());
            auto firstB = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->B.getView());
            auto* fieldAllocation = particles_m->E.getView().data();
            Kokkos::deep_copy(particles_m->E.getView(), Vector_t<double, 3>(123.0));
            Kokkos::deep_copy(particles_m->B.getView(), Vector_t<double, 3>(456.0));
            EXPECT_EQ(solver.solve(context).backendSolves, 4u);
            EXPECT_EQ(particles_m->E.getView().data(), fieldAllocation);
            auto repeatedE = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->E.getView());
            auto repeatedB = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), particles_m->B.getView());
            auto otherE =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), secondary.E.getView());
            auto otherB =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), secondary.B.getView());
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_DOUBLE_EQ(otherE(0)[d], 7.0);
                EXPECT_DOUBLE_EQ(otherB(0)[d], 9.0);
                for (std::size_t i = 0; i < 2; ++i) {
                    EXPECT_TRUE(std::isfinite(repeatedE(i)[d]));
                    EXPECT_NEAR(repeatedE(i)[d], firstE(i)[d], 1.0e-12);
                    EXPECT_NEAR(repeatedB(i)[d], firstB(i)[d], 1.0e-12);
                }
            }
        }

    }  // namespace
}  // namespace opalx::spacecharge
