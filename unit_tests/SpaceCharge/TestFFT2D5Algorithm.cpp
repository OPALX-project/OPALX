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
#include <iomanip>
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

            using Vector                            = Vector_t<double, 3>;
            static constexpr std::size_t sourceSize = 8;
            using Vectors                           = std::array<Vector, sourceSize>;

            static Vectors sourcePositions() {
                return {
                        {{-0.12, 0.03, 0.25},
                         {0.08, -0.04, 0.25},
                         {-0.06, -0.07, 0.40},
                         {0.10, 0.05, 0.40},
                         {-0.09, 0.06, 0.55},
                         {0.07, -0.05, 0.55},
                         {-0.04, -0.02, 0.70},
                         {0.11, 0.04, 0.70}}};
            }

            void writePath(const CoordinateSystemTrafo& labRotation) const {
                std::ofstream path(pathFile_m);
                path << std::setprecision(17);
                path << "s rx ry rz px py pz ex ey ez bx by bz ke t\n";
                for (const double s : {0.0, 1.0}) {
                    const Vector point = labRotation.transformTo(Vector(0.0, 0.0, s));
                    path << s << " " << point[0] << " " << point[1] << " " << point[2]
                         << " 0 0 1 0 0 0 0 0 0 0 0\n";
                }
            }

            void setPhysicalSource(
                    const Vectors& labPositions, const Vector& labMomentum,
                    const CoordinateSystemTrafo& referenceToLab,
                    ParticleContainer* container = nullptr) {
                auto& particles = container == nullptr ? *particles_m : *container;
                if (particles.getLocalNum() == 0) particles.createParticles(sourceSize);
                particles.setQ(1.0e-12);
                particles.setM(Physics::m_e);
                particles.setToLabTrafo(referenceToLab);
                auto positions = particles.R.getHostMirror();
                auto momenta   = particles.P.getHostMirror();
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    positions(i) = referenceToLab.transformFrom(labPositions[i]);
                    momenta(i)   = referenceToLab.rotateFrom(labMomentum);
                }
                Kokkos::deep_copy(particles.R.getView(), positions);
                Kokkos::deep_copy(particles.P.getView(), momenta);
                Kokkos::deep_copy(particles.dt.getView(), 1.0e-12);
                Kokkos::deep_copy(particles.InvalidMask.getView(), false);
                particles.markMomentsDirty();
                particles.updateMoments();
            }

            SpaceChargeStepState alignedStep(const double originZ) const {
                SpaceChargeStepState step;
                step.timeStep = 1.0e-12;
                const CoordinateSystemTrafo trackerToSolve(
                        Vector(0.0, 0.0, originZ),
                        getQuaternion(particles_m->getMeanP(), Vector(0.0, 0.0, 1.0)));
                step.frames = {trackerToSolve, trackerToSolve.inverted()};
                return step;
            }

            Vectors physicalFields(
                    const bool electric, ParticleContainer* container = nullptr) const {
                const auto& particles = container == nullptr ? *particles_m : *container;
                const auto fields     = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(),
                        electric ? particles.E.getView() : particles.B.getView());
                Vectors result;
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    result[i] = particles.getToLabTrafo().rotateTo(fields(i));
                }
                return result;
            }

            void expectDepositedCharge(
                    const FFT2D5Algorithm& solver, const double sourceCopies = 1.0) const {
                const auto* rho = solver.getRho();
                const auto host =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rho->getView());
                const auto spacing = rho->get_mesh().getMeshSpacing();
                const auto ng      = static_cast<std::size_t>(rho->getNghost());
                double densitySum  = 0.0;
                for (std::size_t i = ng; i < host.extent(0) - ng; ++i) {
                    for (std::size_t j = ng; j < host.extent(1) - ng; ++j) {
                        for (std::size_t k = ng; k < host.extent(2) - ng; ++k) {
                            densitySum += host(i, j, k);
                        }
                    }
                }
                const double charge = densitySum * spacing[0] * spacing[1] * spacing[2];
                const double expected =
                        sourceCopies * sourceSize * particles_m->getChargePerParticle();
                EXPECT_NEAR(charge, expected, std::abs(expected) * 1.0e-13);
            }

            void expectSourceRestored(
                    const Vectors& labPositions, const Vector& labMomentum,
                    ParticleContainer* container = nullptr) const {
                const auto& particles = container == nullptr ? *particles_m : *container;
                const auto positions  = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.R.getView());
                const auto momenta = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.P.getView());
                const auto timeSteps = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace(), particles.dt.getView());
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    const Vector r = particles.getToLabTrafo().transformTo(positions(i));
                    const Vector p = particles.getToLabTrafo().rotateTo(momenta(i));
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(r[d], labPositions[i][d], 1.0e-13);
                        EXPECT_NEAR(p[d], labMomentum[d], 1.0e-13);
                    }
                    EXPECT_NEAR(timeSteps(i), 1.0e-12, 1.0e-26);
                }
            }

            std::filesystem::path pathFile_m;
            std::unique_ptr<Mesh_t<3>> mesh_m;
            std::unique_ptr<FieldLayout_t<3>> layout_m;
            std::shared_ptr<ParticleContainer> particles_m;
            std::shared_ptr<BunchStateHandler> bunchState_m;
            CoordinateSystemTrafo trackerToSolve_m;
            CoordinateSystemTrafo solveToTracker_m;
        };

        TEST_F(FFT2D5AlgorithmTest, LazyInitializationSkipsEmptySlices) {
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
            EXPECT_EQ(result.backendSolves, 0u);
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

        TEST_F(FFT2D5AlgorithmTest, LabPathFieldsSurviveReferenceFrameChanges) {
            const Vectors positions = sourcePositions();
            const Vector momentum(0.0, 0.0, 0.4);
            setPhysicalSource(positions, momentum, CoordinateSystemTrafo());
            auto values          = config(FFT2D5LongitudinalFieldMode::Open);
            values.grid.meshSize = {16, 16, 16};
            std::array particles{particles_m.get()};
            FFT2D5Algorithm solver(values, particles, bunchState_m);
            const std::array<std::uint8_t, 1> activity{1};
            ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, alignedStep(0.0)})));
            expectDepositedCharge(solver);
            const Vectors expectedE = physicalFields(true);
            const Vectors expectedB = physicalFields(false);
            double electricNorm = 0.0, magneticNorm = 0.0, longitudinalNorm = 0.0;
            for (std::size_t i = 0; i < sourceSize; ++i) {
                electricNorm += expectedE[i].dot(expectedE[i]);
                magneticNorm += expectedB[i].dot(expectedB[i]);
                longitudinalNorm += expectedE[i][2] * expectedE[i][2];
            }
            ASSERT_GT(electricNorm, 0.0);
            ASSERT_GT(magneticNorm, 0.0);
            ASSERT_GT(longitudinalNorm, 0.0);

            // The physical particles and lab path stay fixed; only their stored reference
            // coordinates and the tracker's aligned solve frame change between calls.
            for (const double angle : {0.0, 0.45, 1.1}) {
                const CoordinateSystemTrafo referenceToLab(
                        Vector(0.2, 0.0, -0.3),
                        Quaternion(std::cos(angle / 2), 0.0, std::sin(angle / 2), 0.0));
                setPhysicalSource(positions, momentum, referenceToLab);
                ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, alignedStep(-0.37)})));
                expectSourceRestored(positions, momentum);
                expectDepositedCharge(solver);
                const Vectors actualE = physicalFields(true);
                const Vectors actualB = physicalFields(false);
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(actualE[i][d], expectedE[i][d], 1.0e-10);
                        EXPECT_NEAR(actualB[i][d], expectedB[i][d], 1.0e-18);
                    }
                }
            }
        }

        TEST_F(FFT2D5AlgorithmTest, RigidPlanarRotationPreservesLocalFields) {
            const Vectors source = sourcePositions();
            Vectors expectedE, expectedB;
            const std::array<std::uint8_t, 1> activity{1};
            for (const double degrees : {0.0, 45.0, 75.0}) {
                SCOPED_TRACE(degrees);
                const double angle = degrees * Physics::pi / 180;
                const CoordinateSystemTrafo labRotation(
                        Vector(0.0),
                        Quaternion(std::cos(angle / 2), 0.0, std::sin(angle / 2), 0.0));
                writePath(labRotation);
                Vectors positions;
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    positions[i] = labRotation.transformTo(source[i]);
                }
                const Vector momentum = labRotation.rotateTo(Vector(0.0, 0.0, 0.4));
                const CoordinateSystemTrafo referenceToLab =
                        labRotation * CoordinateSystemTrafo(Vector(0.0, 0.0, -0.5), Quaternion());
                setPhysicalSource(positions, momentum, referenceToLab);
                auto values          = config(FFT2D5LongitudinalFieldMode::Open);
                values.grid.meshSize = {16, 16, 16};
                std::array particles{particles_m.get()};
                FFT2D5Algorithm solver(values, particles, bunchState_m);
                ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, alignedStep(-0.5)})));
                expectSourceRestored(positions, momentum);
                expectDepositedCharge(solver);
                const Vectors actualE = physicalFields(true);
                const Vectors actualB = physicalFields(false);
                double fieldNorm      = 0.0;
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    const Vector e = labRotation.rotateFrom(actualE[i]);
                    const Vector b = labRotation.rotateFrom(actualB[i]);
                    fieldNorm += e.dot(e);
                    if (degrees == 0.0) {
                        expectedE[i] = e;
                        expectedB[i] = b;
                    } else {
                        for (unsigned d = 0; d < 3; ++d) {
                            EXPECT_NEAR(e[d], expectedE[i][d], 1.0e-10);
                            EXPECT_NEAR(b[d], expectedB[i][d], 1.0e-18);
                        }
                    }
                }
                EXPECT_GT(fieldNorm, 0.0);
            }
        }

        TEST_F(FFT2D5AlgorithmTest, ActiveContainersShareLabPathAcrossReferencePoses) {
            const Vectors positions = sourcePositions();
            const Vector momentum(0.0, 0.0, 0.4);
            auto values          = config(FFT2D5LongitudinalFieldMode::Open);
            values.grid.meshSize = {16, 16, 16};
            setPhysicalSource(positions, momentum, CoordinateSystemTrafo());
            particles_m->setQ(2.0e-12);
            Vectors expectedE, expectedB;
            {
                std::array particles{particles_m.get()};
                FFT2D5Algorithm solver(values, particles, bunchState_m);
                const std::array<std::uint8_t, 1> activity{1};
                ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, alignedStep(0.0)})));
                expectDepositedCharge(solver);
                expectedE = physicalFields(true);
                expectedB = physicalFields(false);
            }

            // Split the same physical source charge across two containers whose particles
            // are stored in different reference frames. The combined source is unchanged.
            ParticleContainer secondary(*mesh_m, *layout_m);
            secondary.setBunchStateHandler(bunchState_m);
            setPhysicalSource(
                    positions, momentum,
                    CoordinateSystemTrafo(Vector(0.0, 0.0, -0.1), Quaternion()));
            setPhysicalSource(
                    positions, momentum,
                    CoordinateSystemTrafo(
                            Vector(0.2, 0.0, -0.4),
                            Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0)),
                    &secondary);
            std::array particles{particles_m.get(), &secondary};
            FFT2D5Algorithm solver(values, particles, bunchState_m);
            const std::array<std::uint8_t, 2> activity{1, 1};
            ASSERT_NO_THROW(static_cast<void>(solver.solve({activity, alignedStep(-0.4)})));
            expectDepositedCharge(solver, 2.0);
            for (auto* container : particles) {
                expectSourceRestored(positions, momentum, container);
                const Vectors actualE = physicalFields(true, container);
                const Vectors actualB = physicalFields(false, container);
                for (std::size_t i = 0; i < sourceSize; ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        EXPECT_NEAR(actualE[i][d], expectedE[i][d], 1.0e-10);
                        EXPECT_NEAR(actualB[i][d], expectedB[i][d], 1.0e-18);
                    }
                }
            }
        }

        TEST_F(FFT2D5AlgorithmTest, CoverageFailureRestoresReferenceCoordinatesAndTimeSteps) {
            Vectors positions = sourcePositions();
            positions[0][2]   = 1.2;
            const Vector momentum(0.0, 0.0, 0.4);
            const CoordinateSystemTrafo referenceToLab(
                    Vector(0.1, 0.0, -0.2), Quaternion(std::cos(0.3), 0.0, std::sin(0.3), 0.0));
            setPhysicalSource(positions, momentum, referenceToLab);
            std::array particles{particles_m.get()};
            FFT2D5Algorithm solver(
                    config(FFT2D5LongitudinalFieldMode::Open), particles, bunchState_m);
            const std::array<std::uint8_t, 1> activity{1};
            EXPECT_THROW(
                    static_cast<void>(solver.solve({activity, alignedStep(-0.5)})), OpalException);
            expectSourceRestored(positions, momentum);
        }

        TEST_F(FFT2D5AlgorithmTest, PreparedProjectionsMatchUncachedStagesAfterSourceChanges) {
            // Nonuniform segments, a repeated vertex and a closed seam exercise geometry
            // and first-nearest selection. Both boundary modes use this same spatial path.
            {
                std::ofstream path(pathFile_m);
                path << std::setprecision(17);
                path << "s rx ry rz px py pz ex ey ez bx by bz ke t\n";
                for (const auto& point : std::array<Vector, 8>{
                             Vector(0, 0, 0), Vector(0, 0, 0.37), Vector(0, 0, 1), Vector(-1, 0, 1),
                             Vector(-1, 0, 1), Vector(-1, 0, 0), Vector(-0.4, 0, 0),
                             Vector(0, 0, 0)}) {
                    path << "0 " << point[0] << " " << point[1] << " " << point[2]
                         << " 0 0 1 0 0 0 0 0 0 0 0\n";
                }
            }
            setPhysicalSource(sourcePositions(), Vector(0, 0, 0.4), CoordinateSystemTrafo());
            const std::array<std::uint8_t, 1> activity{1};
            SpaceChargeStepState step;
            step.timeStep = 1.0e-12;
            SpaceChargeSolveContext context(activity, step);
            size_t iteration = 0;
            for (const bool closed : {false, true}) {
                auto values          = config(FFT2D5LongitudinalFieldMode::Open);
                values.grid.meshSize = {16, 16, 32};
                values.closedRing    = closed;
                std::array containers{particles_m.get()};
                FFT2D5Algorithm solver(values, containers, bunchState_m);
                for (unsigned repeat = 0; repeat < 2; ++repeat, ++iteration) {
                    const size_t count = sourceSize + 4 * iteration;
                    if (particles_m->getLocalNum() < count) {
                        particles_m->createParticles(count - particles_m->getLocalNum());
                    }
                    particles_m->setQ(1.0e-12);
                    particles_m->setM(Physics::m_e);
                    auto positions = particles_m->R.getHostMirror();
                    auto momenta   = particles_m->P.getHostMirror();
                    auto mask      = particles_m->InvalidMask.getHostMirror();
                    for (size_t n = 0; n < count; ++n) {
                        // Reordering, movement and momentum changes must refresh every solve.
                        positions(n) = sourcePositions()[(n + iteration) % sourceSize];
                        positions(n)[0] += 0.002 * iteration;
                        momenta(n) = Vector(0.002 * iteration, 0.001 * n, 0.4);
                        mask(n)    = n + 1 == count;
                    }
                    positions(0) = Vector(-0.1, 0.01, 0);    // Last segment near periodic seam.
                    positions(1) = Vector(0, -0.01, 0.001);  // First segment near periodic seam.
                    positions(2) = Vector(0, 0, 1);          // Equal-distance corner tie.
                    Kokkos::deep_copy(particles_m->R.getView(), positions);
                    Kokkos::deep_copy(particles_m->P.getView(), momenta);
                    Kokkos::deep_copy(particles_m->InvalidMask.getView(), mask);
                    Kokkos::deep_copy(particles_m->dt.getView(), step.timeStep);
                    particles_m->markMomentsDirty();
                    particles_m->updateMoments();
                    ASSERT_NO_THROW(static_cast<void>(solver.solve(context)));
                    // create_mirror always allocates; create_mirror_view could alias on CPU.
                    auto expectedE   = Kokkos::create_mirror(particles_m->E.getView());
                    auto expectedB   = Kokkos::create_mirror(particles_m->B.getView());
                    auto expectedRho = Kokkos::create_mirror(solver.getRho()->getView());
                    Kokkos::deep_copy(expectedE, particles_m->E.getView());
                    Kokkos::deep_copy(expectedB, particles_m->B.getView());
                    Kokkos::deep_copy(expectedRho, solver.getRho()->getView());
                    solver.scatterToGrid(context);  // Uncached original conversion is the oracle.
                    const auto actualRho = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), solver.getRho()->getView());
                    double rhoScale = 0.0;
                    for (size_t i = 0; i < actualRho.extent(0); ++i)
                        for (size_t j = 0; j < actualRho.extent(1); ++j)
                            for (size_t k = 0; k < actualRho.extent(2); ++k)
                                rhoScale = std::max(rhoScale, std::abs(expectedRho(i, j, k)));
                    for (size_t i = 0; i < actualRho.extent(0); ++i)
                        for (size_t j = 0; j < actualRho.extent(1); ++j)
                            for (size_t k = 0; k < actualRho.extent(2); ++k)
                                EXPECT_NEAR(
                                        actualRho(i, j, k), expectedRho(i, j, k),
                                        rhoScale * 1e-12 + 1e-30);
                    solver.solvePoissons();
                    solver.calculateLineDensity();
                    solver.gatherFromGrid(context);
                    const auto actualE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), particles_m->E.getView());
                    const auto actualB = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), particles_m->B.getView());
                    double eScale = 0.0, bScale = 0.0;
                    for (size_t n = 0; n < count; ++n)
                        for (unsigned d = 0; d < 3; ++d) {
                            eScale = std::max(eScale, std::abs(expectedE(n)[d]));
                            bScale = std::max(bScale, std::abs(expectedB(n)[d]));
                        }
                    for (size_t n = 0; n < count; ++n)
                        for (unsigned d = 0; d < 3; ++d) {
                            EXPECT_NEAR(actualE(n)[d], expectedE(n)[d], eScale * 1e-12 + 1e-25);
                            EXPECT_NEAR(actualB(n)[d], expectedB(n)[d], bScale * 1e-12 + 1e-25);
                        }
                }
            }
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
