/**
 * @file TestBendRep.cpp
 * @brief Regression tests for the OPALX-native SBEND/RBEND representation layer.
 *
 * The tests cover the first port of the analytic bend representation:
 * geometry ownership, dipole-field channels, and the curved placement hooks
 * used by the element-placement bridge.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include "AbsBeamline/BeamlineVisitor.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Channels/Channel.h"
#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Steppers/BorisPusher.h"

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace {
    struct ReferenceState {
        Vector_t<double, 3> r;
        Vector_t<double, 3> p;
    };

    using ParticleContainer_t = ParticleContainer<double, 3>;

    class BendRepParticleContainerTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        /**
         * @brief Create a deterministic particle container with fixed local coordinates.
         *
         * The helper constructs a small particle container in the bend-local chart
         * and populates it with the supplied positions. Momenta and fields are
         * initialized to zero so the bend field evaluation can be compared
         * particle-by-particle against the scalar host path.
         */
        std::shared_ptr<ParticleContainer_t> makeContainer(
                const std::vector<std::array<double, 3>>& positions) {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> rmin   = -4.0;
            ippl::Vector<double, 3> rmax   = 4.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; ++i) {
                domain[i] = ippl::Index(nr[i]);
            }

            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);

            auto pc = std::make_shared<ParticleContainer_t>(mesh, fl);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            pc->createParticles(positions.size());

            auto Rhost = pc->R.getHostMirror();
            auto Phost = pc->P.getHostMirror();
            auto Ehost = pc->E.getHostMirror();
            auto Bhost = pc->B.getHostMirror();

            for (std::size_t i = 0; i < positions.size(); ++i) {
                Rhost(i)[0] = positions[i][0];
                Rhost(i)[1] = positions[i][1];
                Rhost(i)[2] = positions[i][2];
                Phost(i)    = Vector_t<double, 3>(0.0);
                Ehost(i)    = Vector_t<double, 3>(0.0);
                Bhost(i)    = Vector_t<double, 3>(0.0);
            }

            Kokkos::deep_copy(pc->R.getView(), Rhost);
            Kokkos::deep_copy(pc->P.getView(), Phost);
            Kokkos::deep_copy(pc->E.getView(), Ehost);
            Kokkos::deep_copy(pc->B.getView(), Bhost);
            Kokkos::fence();

            return pc;
        }

        std::shared_ptr<ParticleContainer_t> makeContainer(const std::array<double, 3>& position) {
            return makeContainer(std::vector<std::array<double, 3>>{position});
        }
    };

    class BendDispatchVisitor : public BeamlineVisitor {
    public:
        void execute() override {}
        void visitBeamline(const Beamline&) override {}
        void visitComponent(const Component&) override { visitedComponent = true; }
        void visitConstantEFieldCavity(const ConstantEFieldCavity&) override {}
        void visitDrift(const Drift&) override {}
        void visitFlaggedElmPtr(const FlaggedElmPtr&) override {}
        void visitLaser(const Laser&) override {}
        void visitMarker(const Marker&) override {}
        void visitMonitor(const Monitor&) override {}
        void visitMultipole(const Multipole&) override {}
        void visitMultipoleT(const MultipoleT&) override {}
        void visitRFCavity(const RFCavity&) override {}
        void visitRBend(const RBend&) override { visitedRBend = true; }
        void visitVariableRFCavity(const VariableRFCavity&) override {}
        void visitScalingFFAMagnet(const ScalingFFAMagnet&) override {}
        void visitRing(const Ring&) override {}
        void visitSolenoid(const Solenoid&) override {}
        void visitSBend(const SBend&) override { visitedSBend = true; }
        void visitTravelingWave(const TravelingWave&) override {}
        void visitVerticalFFAMagnet(const VerticalFFAMagnet&) override {}
        void visitProbe(const Probe&) override {}

        bool visitedComponent = false;
        bool visitedSBend     = false;
        bool visitedRBend     = false;
    };

    /**
     * @brief Compute the dimensionless momentum \f$\beta\gamma\f$ from kinetic energy.
     *
     * The relativistic momentum variable used by OPALX tracking is
     * \f[
     * \mathbf{p} = \beta\gamma \hat{\mathbf{u}},
     * \f]
     * where \f$\gamma = 1 + E_\mathrm{kin}/m c^2\f$.
     */
    double betaGammaFromKineticEnergy(const double kineticEnergyEV, const double massEV) {
        const double gamma = 1.0 + kineticEnergyEV / massEV;
        return std::sqrt(gamma * gamma - 1.0);
    }

    /**
     * @brief Relativistic magnetic rigidity relation used for the bend field.
     *
     * For motion in a pure dipole field the design radius satisfies
     * \f[
     * \rho = \frac{\beta\gamma m}{|q| c |B|},
     * \f]
     * so the field magnitude required for a target radius is
     * \f[
     * |B| = \frac{\beta\gamma m}{|q| c \rho}.
     * \f]
     */
    double dipoleFieldMagnitude(
            const double betaGamma, const double massEV, const double charge, const double radius) {
        return betaGamma * massEV / (Physics::c * std::abs(charge) * radius);
    }

    Vector_t<double, 3> fieldLocalToEntryCartesian(
            const Vector_t<double, 3>& fieldLocal, const double curvature,
            const double bodyLength) {
        if (std::abs(curvature) <= 1.0e-15) {
            return fieldLocal;
        }

        if (fieldLocal(2) <= 0.0) {
            return fieldLocal;
        }

        if (fieldLocal(2) >= bodyLength) {
            const double exitPhi = curvature * bodyLength;
            const double c       = std::cos(exitPhi);
            const double sinExit = std::sin(exitPhi);
            const double ds      = fieldLocal(2) - bodyLength;
            return Vector_t<double, 3>(
                    (c - 1.0) / curvature + c * fieldLocal(0) - sinExit * ds,
                    fieldLocal(1), sinExit / curvature + sinExit * fieldLocal(0) + c * ds);
        }

        const double radius = 1.0 / curvature;
        const double phi    = curvature * fieldLocal(2);
        return Vector_t<double, 3>(
                (radius + fieldLocal(0)) * std::cos(phi) - radius, fieldLocal(1),
                (radius + fieldLocal(0)) * std::sin(phi));
    }

    Vector_t<double, 3> fieldLocalVectorToEntryCartesian(
            const Vector_t<double, 3>& fieldVector, const double curvature, const double s,
            const double bodyLength) {
        if (std::abs(curvature) <= 1.0e-15) {
            return fieldVector;
        }

        const double phi = curvature * ((s >= bodyLength) ? bodyLength : s);
        if (s <= 0.0) {
            return fieldVector;
        }
        return Vector_t<double, 3>(
                std::cos(phi) * fieldVector(0) - std::sin(phi) * fieldVector(2), fieldVector(1),
                std::sin(phi) * fieldVector(0) + std::cos(phi) * fieldVector(2));
    }

    Vector_t<double, 3> entryCartesianToFieldLocal(
            const Vector_t<double, 3>& entryCartesian, const double curvature,
            const double bodyLength) {
        if (std::abs(curvature) <= 1.0e-15) {
            return entryCartesian;
        }

        if (entryCartesian(2) <= 0.0) {
            return entryCartesian;
        }

        const double exitPhi = curvature * bodyLength;
        const double c       = std::cos(exitPhi);
        const double sinExit = std::sin(exitPhi);
        const double exitX   = (c - 1.0) / curvature;
        const double exitZ   = sinExit / curvature;
        const double dx      = entryCartesian(0) - exitX;
        const double dz      = entryCartesian(2) - exitZ;
        const double exitLocalX = c * dx + sinExit * dz;
        const double exitLocalZ = -sinExit * dx + c * dz;
        if (exitLocalZ >= 0.0) {
            return Vector_t<double, 3>(exitLocalX, entryCartesian(1), bodyLength + exitLocalZ);
        }

        const double radius = 1.0 / curvature;
        const double phi    = std::atan2(entryCartesian(2), entryCartesian(0) + radius);
        const double s      = phi / curvature;
        const double radialDistance =
                std::hypot(entryCartesian(0) + radius, entryCartesian(2)) - std::abs(radius);
        return Vector_t<double, 3>(radialDistance, entryCartesian(1), s);
    }

    ReferenceState evaluateDerivative(
            BendBase& bend, const ReferenceState& state, const double massEV, const double charge) {
        Vector_t<double, 3> E(0.0);
        Vector_t<double, 3> B(0.0);
        if (bend.applyToReferenceParticle(state.r, state.p, 0.0, E, B)) {
            throw std::runtime_error(
                    "Reference particle left bend field support during regression.");
        }

        const double gamma = std::sqrt(1.0 + dot(state.p, state.p));
        ReferenceState derivative{};
        derivative.r = Physics::c * state.p / gamma;
        derivative.p = charge * Physics::c / massEV * E
                       + charge * Physics::c * Physics::c / (gamma * massEV) * cross(state.p, B);
        return derivative;
    }

    ReferenceState integrateReferenceParticle(
            BendBase& bend, const double massEV, const double charge, const double kineticEnergyEV,
            const double integrationTime, const std::size_t steps, const double initialZ = 0.0) {
        ReferenceState state{};
        state.r = Vector_t<double, 3>(0.0, 0.0, initialZ);
        state.p =
                Vector_t<double, 3>(0.0, 0.0, betaGammaFromKineticEnergy(kineticEnergyEV, massEV));

        const double dt = integrationTime / static_cast<double>(steps);
        for (std::size_t i = 0; i < steps; ++i) {
            const ReferenceState k1 = evaluateDerivative(bend, state, massEV, charge);

            ReferenceState tmp = state;
            tmp.r += 0.5 * dt * k1.r;
            tmp.p += 0.5 * dt * k1.p;
            const ReferenceState k2 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += 0.5 * dt * k2.r;
            tmp.p += 0.5 * dt * k2.p;
            const ReferenceState k3 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += dt * k3.r;
            tmp.p += dt * k3.p;
            const ReferenceState k4 = evaluateDerivative(bend, tmp, massEV, charge);

            state.r += dt / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
            state.p += dt / 6.0 * (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p);
        }

        return state;
    }

    ReferenceState integrateReferenceParticle(
            BendBase& bend, const double massEV, const double charge,
            const ReferenceState& initialState, const double integrationTime,
            const std::size_t steps) {
        ReferenceState state = initialState;

        const double dt = integrationTime / static_cast<double>(steps);
        for (std::size_t i = 0; i < steps; ++i) {
            const ReferenceState k1 = evaluateDerivative(bend, state, massEV, charge);

            ReferenceState tmp = state;
            tmp.r += 0.5 * dt * k1.r;
            tmp.p += 0.5 * dt * k1.p;
            const ReferenceState k2 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += 0.5 * dt * k2.r;
            tmp.p += 0.5 * dt * k2.p;
            const ReferenceState k3 = evaluateDerivative(bend, tmp, massEV, charge);

            tmp = state;
            tmp.r += dt * k3.r;
            tmp.p += dt * k3.p;
            const ReferenceState k4 = evaluateDerivative(bend, tmp, massEV, charge);

            state.r += dt / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
            state.p += dt / 6.0 * (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p);
        }

        return state;
    }

    void kickParticleContainer(
            ParticleContainer_t& pc, const char* label, const double massEV, const double charge) {
        BorisPusher pusher;
        auto Pview  = pc.P.getView();
        auto dtview = pc.dt.getView();
        auto Eview  = pc.E.getView();
        auto Bview  = pc.B.getView();
        Kokkos::parallel_for(
                label, pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                    Vector_t<double, 3> p = Pview(i);
                    pusher.kick(
                            Vector_t<double, 3>(0.0), p, Eview(i), Bview(i), dtview(i), massEV,
                            charge);
                    Pview(i) = p;
                });
        Kokkos::fence();
    }

    void pushParticleContainer(ParticleContainer_t& pc, const char* label, const double dt) {
        BorisPusher pusher;
        auto Rview = pc.R.getView();
        auto Pview = pc.P.getView();
        Kokkos::parallel_for(
                label, pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                    Vector_t<double, 3> r = Rview(i);
                    pusher.push(r, Pview(i), dt);
                    Rview(i) = r;
                });
        Kokkos::fence();
    }

    TEST(SBendRep, GeometryAndDesignPathFollowPlanarArc) {
        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(2.0, M_PI / 12.0);
        bend.setElementLength(2.0);
        bend.setBendAngle(M_PI / 6.0);
        bend.setExitAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::SBEND);
        EXPECT_NEAR(bend.getGeometry().getElementLength(), 2.0, 1e-12);
        EXPECT_NEAR(bend.getBendAngle(), M_PI / 6.0, 1e-12);

        const double radius = 2.0 / (M_PI / 6.0);
        const double chord  = 2.0 * radius * std::sin(M_PI / 12.0);
        EXPECT_NEAR(bend.getChordLength(), chord, 1e-12);

        const auto path = bend.getDesignPath(7);
        ASSERT_GE(path.size(), 7u);

        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), entry(0), 1e-12);
        EXPECT_NEAR(path.front()(1), entry(1), 1e-12);
        EXPECT_NEAR(path.front()(2), entry(2), 1e-12);
        EXPECT_NEAR(path.back()(0), exit(0), 1e-12);
        EXPECT_NEAR(path.back()(1), exit(1), 1e-12);
        EXPECT_NEAR(path.back()(2), exit(2), 1e-12);
    }

    TEST(SBendRep, DipoleChannelUsesCurrentZeroBasedMultipoleConvention) {
        SBendRep bend("SBEND");

        Channel* by = bend.getChannel("BY", false);
        ASSERT_NE(by, nullptr);
        EXPECT_TRUE(by->set(1.75));
        EXPECT_DOUBLE_EQ(bend.getB(), 1.75);
        EXPECT_DOUBLE_EQ(bend.getField().getNormalComponent(0), 1.75);
        delete by;
    }

    TEST(SBendRep, AcceptDispatchesToVisitSBend) {
        SBendRep bend("SBEND");
        BendDispatchVisitor visitor;
        bend.accept(visitor);
        EXPECT_TRUE(visitor.visitedSBend);
        EXPECT_FALSE(visitor.visitedRBend);
        EXPECT_FALSE(visitor.visitedComponent);
    }

    TEST(RBendRep, ExitAngleTracksRectangularBendConvention) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(1.2);
        bend.getGeometry().setBendAngle(M_PI / 9.0);
        bend.setElementLength(1.2);
        bend.setBendAngle(M_PI / 9.0);
        bend.setEntranceAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::RBEND);
        EXPECT_NEAR(bend.getExitAngle(), M_PI / 18.0, 1e-12);
        EXPECT_NEAR(bend.getChordLength(), 1.2, 1e-12);

        const auto path = bend.getDesignPath(5);
        ASSERT_GE(path.size(), 5u);
        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), 0.0, 1e-12);
        EXPECT_NEAR(path.back()(0), 0.0, 1e-12);
        EXPECT_NEAR(entry(0), exit(0), 1e-12);
        EXPECT_LT(entry(0), 0.0);
    }

    TEST(RBendRep, ClonePreservesGeometryAndField) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(0.8);
        bend.getGeometry().setBendAngle(M_PI / 8.0);
        bend.setElementLength(0.8);
        bend.setBendAngle(M_PI / 8.0);
        bend.setB(0.9);

        std::unique_ptr<ElementBase> clone(bend.clone());
        auto* copy = dynamic_cast<RBendRep*>(clone.get());
        ASSERT_NE(copy, nullptr);
        EXPECT_NEAR(copy->getElementLength(), 0.8, 1e-12);
        EXPECT_NEAR(copy->getBendAngle(), M_PI / 8.0, 1e-12);
        EXPECT_NEAR(copy->getB(), 0.9, 1e-12);
    }

    TEST(RBendRep, AcceptDispatchesToVisitRBend) {
        RBendRep bend("RBEND");
        BendDispatchVisitor visitor;
        bend.accept(visitor);
        EXPECT_TRUE(visitor.visitedRBend);
        EXPECT_FALSE(visitor.visitedSBend);
        EXPECT_FALSE(visitor.visitedComponent);
    }

    TEST(RBendRep, EffectiveFieldLengthUsesReferencePathInsteadOfStraightBody) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double referencePathLength =
                bodyLength * (bendAngle / 2.0) / std::sin(bendAngle / 2.0);
        const double profileGap           = 2.0 * fringeHalfGap;
        const double fringeHalfWidth      = 5.0 * profileGap;
        const double expectedEntrySupport = fringeHalfWidth;
        const double expectedExitSupport =
                fringeHalfWidth / std::abs(std::cos(bend.getExitAngle()));

        EXPECT_NEAR(bend.getChordLength(), bodyLength, 1.0e-12);
        EXPECT_NEAR(bend.getGeometry().getArcLength(), referencePathLength, 1.0e-12);
        EXPECT_NEAR(bend.getEntryFringeSupportLength(), expectedEntrySupport, 1.0e-12);
        EXPECT_NEAR(bend.getExitFringeSupportLength(), expectedExitSupport, 1.0e-12);
        EXPECT_NEAR(bend.getEffectiveFieldLength(), referencePathLength, 2.0e-6);
        EXPECT_GT(bend.getEffectiveFieldLength(), bend.getElementLength());
    }

    TEST(SBendRep, TrackingSlicesCoverBodyAndFringeSupport) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr double e1             = Physics::pi / 10.0;
        constexpr double e2             = Physics::pi / 12.0;
        constexpr std::size_t nSlices   = 64;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(e1);
        bend.setExitAngle(e2);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);

        const auto slices = bend.buildTrackingSlices();
        ASSERT_EQ(slices.size(), nSlices);

        const double supportBegin = -bend.getEntryFringeSupportLength();
        const double supportEnd   = bodyLength + bend.getExitFringeSupportLength();
        const double ds           = (supportEnd - supportBegin) / static_cast<double>(nSlices);

        EXPECT_NEAR(slices.front().sBegin, supportBegin, 1.0e-12);
        EXPECT_NEAR(slices.back().sEnd, supportEnd, 1.0e-12);
        EXPECT_NEAR(slices.front().sCenter, supportBegin + 0.5 * ds, 1.0e-12);
        EXPECT_NEAR(slices.back().sCenter, supportEnd - 0.5 * ds, 1.0e-12);

        const Vector_t<double, 3> frontTangent =
                slices.front().entryToSliceLocal.rotateFrom(Vector_t<double, 3>(0.0, 0.0, 1.0));
        EXPECT_NEAR(frontTangent(0), 0.0, 1.0e-12);
        EXPECT_NEAR(frontTangent(1), 0.0, 1.0e-12);
        EXPECT_NEAR(frontTangent(2), 1.0, 1.0e-12);
        EXPECT_NEAR(slices.front().entryOrigin(0), 0.0, 1.0e-12);
        EXPECT_NEAR(slices.front().entryOrigin(2), slices.front().sCenter, 1.0e-12);

        const Vector_t<double, 3> backTangent =
                slices.back().entryToSliceLocal.rotateFrom(Vector_t<double, 3>(0.0, 0.0, 1.0));
        EXPECT_NEAR(backTangent(0), -std::sin(bendAngle), 1.0e-12);
        EXPECT_NEAR(backTangent(1), 0.0, 1.0e-12);
        EXPECT_NEAR(backTangent(2), std::cos(bendAngle), 1.0e-12);
    }

    TEST(RBendRep, TrackingSlicesUseReferenceArcLengthInsteadOfBodyChord) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 5;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);

        const auto slices = bend.buildTrackingSlices();
        ASSERT_EQ(slices.size(), nSlices);

        const double supportBegin = -bend.getEntryFringeSupportLength();
        const double supportEnd =
                bend.getGeometry().getArcLength() + bend.getExitFringeSupportLength();

        EXPECT_NEAR(slices.front().sBegin, supportBegin, 1.0e-12);
        EXPECT_NEAR(slices.back().sEnd, supportEnd, 1.0e-12);
        EXPECT_GT(slices.back().sEnd, bodyLength);
    }

    TEST_F(BendRepParticleContainerTest,
           SliceKernelMatchesHostFieldAtEquivalentReferenceCoordinate) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 8;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);
        bend.getField().setNormalComponent(1, 0.4);
        bend.getField().setSkewComponent(1, -0.2);
        const double curvature = bendAngle / bodyLength;

        const auto slices = bend.buildTrackingSlices();
        ASSERT_FALSE(slices.empty());
        const BendBase::TrackingSlice& slice = slices.front();
        const Vector_t<double, 3> fieldLocalPosition(1.5e-3, -0.8e-3, slice.sCenter);
        const Vector_t<double, 3> entryCartesian =
                fieldLocalToEntryCartesian(fieldLocalPosition, curvature, bodyLength);
        const Vector_t<double, 3> sliceLocalPosition =
                slice.entryToSliceLocal.transformTo(entryCartesian);

        auto pc = makeContainer(
                {sliceLocalPosition(0), sliceLocalPosition(1), sliceLocalPosition(2)});
        bend.applySlice(pc, slice);

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        Vector_t<double, 3> Eexpected(0.0);
        Vector_t<double, 3> Bexpected(0.0);
        bend.apply(fieldLocalPosition, Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
        const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                Bexpected, curvature, fieldLocalPosition(2), bodyLength);
        const Vector_t<double, 3> BexpectedSlice = slice.entryToSliceLocal.rotateTo(BexpectedEntry);

        EXPECT_NEAR(BhostView(0)(0), BexpectedSlice(0), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(1), BexpectedSlice(1), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(2), BexpectedSlice(2), 1.0e-12);
    }

    TEST_F(BendRepParticleContainerTest,
           SliceKernelMatchesHostFieldForDeterministicManyParticleOffsetSet) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 10;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);
        bend.getField().setNormalComponent(1, 0.4);
        bend.getField().setSkewComponent(1, -0.2);
        const double curvature = bendAngle / bodyLength;

        const auto slices = bend.buildTrackingSlices();
        ASSERT_GE(slices.size(), 3u);
        const BendBase::TrackingSlice& slice                       = slices.at(2);
        const std::vector<Vector_t<double, 3>> fieldLocalPositions = {
                {0.0, 0.0, slice.sCenter},
                {1.0e-3, 0.0, slice.sCenter},
                {-1.0e-3, 0.0, slice.sCenter},
                {0.0, 1.0e-3, slice.sCenter},
                {0.0, -1.0e-3, slice.sCenter}};
        std::vector<std::array<double, 3>> sliceLocalPositions;
        for (const auto& fieldLocalPosition : fieldLocalPositions) {
            const Vector_t<double, 3> entryCartesian = fieldLocalToEntryCartesian(
                    fieldLocalPosition, curvature, bodyLength);
            const Vector_t<double, 3> sliceLocal =
                    slice.entryToSliceLocal.transformTo(entryCartesian);
            sliceLocalPositions.push_back({sliceLocal(0), sliceLocal(1), sliceLocal(2)});
        }

        auto pc = makeContainer(sliceLocalPositions);
        bend.applySlice(pc, slice);

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        for (std::size_t i = 0; i < fieldLocalPositions.size(); ++i) {
            Vector_t<double, 3> Eexpected(0.0);
            Vector_t<double, 3> Bexpected(0.0);
            bend.apply(fieldLocalPositions[i], Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
            const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                    Bexpected, curvature, fieldLocalPositions[i](2), bodyLength);
            const Vector_t<double, 3> BexpectedSlice =
                    slice.entryToSliceLocal.rotateTo(BexpectedEntry);

            EXPECT_NEAR(BhostView(i)(0), BexpectedSlice(0), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(BhostView(i)(1), BexpectedSlice(1), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(BhostView(i)(2), BexpectedSlice(2), 1.0e-12) << "particle " << i;
        }
    }

    TEST_F(BendRepParticleContainerTest,
           EntryFringeVerticalOffsetProducesExpectedHorizontalFieldComponent) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 32;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);

        const double curvature = bendAngle / bodyLength;
        const double entryFringe = bend.getEntryFringeSupportLength();
        const auto slices        = bend.buildTrackingSlices();
        ASSERT_FALSE(slices.empty());

        const Vector_t<double, 3> fieldLocalPosition(0.0, 1.0e-3, -0.5 * entryFringe);
        const Vector_t<double, 3> entryCartesian =
                fieldLocalToEntryCartesian(fieldLocalPosition, curvature, bodyLength);

        const BendBase::TrackingSlice* containingSlice = nullptr;
        for (const auto& slice : slices) {
            if (fieldLocalPosition(2) >= slice.sBegin && fieldLocalPosition(2) <= slice.sEnd) {
                containingSlice = &slice;
                break;
            }
        }
        ASSERT_NE(containingSlice, nullptr);

        const Vector_t<double, 3> sliceLocalPosition =
                containingSlice->entryToSliceLocal.transformTo(entryCartesian);
        auto pc = makeContainer(
                {sliceLocalPosition(0), sliceLocalPosition(1), sliceLocalPosition(2)});

        bend.applySlice(pc, *containingSlice);

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        Vector_t<double, 3> Eexpected(0.0);
        Vector_t<double, 3> Bexpected(0.0);
        bend.apply(fieldLocalPosition, Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
        const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                Bexpected, curvature, fieldLocalPosition(2), bodyLength);
        const Vector_t<double, 3> BexpectedSlice =
                containingSlice->entryToSliceLocal.rotateTo(BexpectedEntry);

        EXPECT_GT(std::abs(BexpectedSlice(0)), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(0), BexpectedSlice(0), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(1), BexpectedSlice(1), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(2), BexpectedSlice(2), 1.0e-12);
    }

    TEST_F(BendRepParticleContainerTest,
           ProtonLikeEntryFringeVerticalOffsetProducesNegativeHorizontalFieldComponent) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        const double curvature = bendAngle / bodyLength;
        const double entryFringe = bend.getEntryFringeSupportLength();
        const auto slices        = bend.buildTrackingSlices();
        ASSERT_FALSE(slices.empty());

        const Vector_t<double, 3> fieldLocalPosition(0.0, 1.0e-3, -0.5 * entryFringe);
        const Vector_t<double, 3> entryCartesian =
                fieldLocalToEntryCartesian(fieldLocalPosition, curvature, bodyLength);

        const BendBase::TrackingSlice* containingSlice = nullptr;
        for (const auto& slice : slices) {
            if (fieldLocalPosition(2) >= slice.sBegin && fieldLocalPosition(2) <= slice.sEnd) {
                containingSlice = &slice;
                break;
            }
        }
        ASSERT_NE(containingSlice, nullptr);

        const Vector_t<double, 3> sliceLocalPosition =
                containingSlice->entryToSliceLocal.transformTo(entryCartesian);
        auto pc = makeContainer(
                {sliceLocalPosition(0), sliceLocalPosition(1), sliceLocalPosition(2)});

        bend.applySlice(pc, *containingSlice);

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        Vector_t<double, 3> Eexpected(0.0);
        Vector_t<double, 3> Bexpected(0.0);
        bend.apply(fieldLocalPosition, Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
        const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                Bexpected, curvature, fieldLocalPosition(2), bodyLength);
        const Vector_t<double, 3> BexpectedSlice =
                containingSlice->entryToSliceLocal.rotateTo(BexpectedEntry);

        EXPECT_LT(BexpectedSlice(0), 0.0);
        EXPECT_LT(BhostView(0)(0), 0.0);
        EXPECT_NEAR(BhostView(0)(0), BexpectedSlice(0), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(1), BexpectedSlice(1), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(2), BexpectedSlice(2), 1.0e-12);
    }

    TEST_F(BendRepParticleContainerTest,
           TrackerStyleSliceCompositionPreservesFringeFieldForVerticalOffset) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);

        const double curvature = bendAngle / bodyLength;
        const double entryFringe = bend.getEntryFringeSupportLength();
        const double referenceS  = -0.5 * entryFringe;
        const Vector_t<double, 3> localReferenceParticle(0.0, 1.0e-3, 0.0);

        const Vector_t<double, 3> referenceOrigin = fieldLocalToEntryCartesian(
                Vector_t<double, 3>(0.0, 0.0, referenceS), curvature, bodyLength);
        const Vector_t<double, 3> tangentInEntry = fieldLocalVectorToEntryCartesian(
                Vector_t<double, 3>(0.0, 0.0, 1.0), curvature, referenceS,
                bodyLength);
        CoordinateSystemTrafo entryToReference(
                referenceOrigin, getQuaternion(tangentInEntry, Vector_t<double, 3>(0.0, 0.0, 1.0)));
        CoordinateSystemTrafo referenceToEntry = entryToReference.inverted();

        auto pc = makeContainer(
                {localReferenceParticle(0), localReferenceParticle(1), localReferenceParticle(2)});
        pc->setToLabTrafo(referenceToEntry);

        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * pc->getToLabTrafo();
            CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
            pc->transformBunch(refToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(sliceLocalToRef);
        }

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        const Vector_t<double, 3> entryCartesian =
                pc->getToLabTrafo().transformTo(localReferenceParticle);
        const Vector_t<double, 3> fieldLocalPosition =
                entryCartesianToFieldLocal(entryCartesian, curvature, bodyLength);
        Vector_t<double, 3> Eexpected(0.0);
        Vector_t<double, 3> Bexpected(0.0);
        bend.apply(fieldLocalPosition, Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
        const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                Bexpected, curvature, fieldLocalPosition(2), bodyLength);
        const Vector_t<double, 3> BexpectedRef = pc->getToLabTrafo().rotateFrom(BexpectedEntry);

        EXPECT_GT(std::abs(BexpectedRef(0)), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(0), BexpectedRef(0), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(1), BexpectedRef(1), 1.0e-12);
        EXPECT_NEAR(BhostView(0)(2), BexpectedRef(2), 1.0e-12);
    }

    TEST_F(BendRepParticleContainerTest,
           FullSliceSequenceMatchesHostFieldForDeterministicManyParticleSet) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);
        bend.getField().setNormalComponent(1, 0.4);
        bend.getField().setSkewComponent(1, -0.2);
        const double curvature = bendAngle / bodyLength;

        const double entryFringe = bend.getEntryFringeSupportLength();
        const double exitFringe  = bend.getExitFringeSupportLength();
        const std::vector<Vector_t<double, 3>> fieldLocalPositions = {
                {0.0, 0.0, -0.5 * entryFringe},
                {1.0e-3, 0.0, 0.05 * bodyLength},
                {-1.0e-3, 0.0, 0.35 * bodyLength},
                {0.0, 1.0e-3, 0.65 * bodyLength},
                {0.0, -1.0e-3, 0.75 * bodyLength},
                {0.8e-3, -0.6e-3, bodyLength + 0.5 * exitFringe}};
        std::vector<std::array<double, 3>> entryPositions;
        for (const auto& fieldLocalPosition : fieldLocalPositions) {
            const Vector_t<double, 3> entryCartesian = fieldLocalToEntryCartesian(
                    fieldLocalPosition, curvature, bodyLength);
            entryPositions.push_back({entryCartesian(0), entryCartesian(1), entryCartesian(2)});
        }
        auto pc = makeContainer(entryPositions);

        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            pc->transformBunch(slice.entryToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(slice.entryToSliceLocal.inverted());
        }

        auto BhostView = pc->B.getHostMirror();
        Kokkos::deep_copy(BhostView, pc->B.getView());

        for (std::size_t i = 0; i < fieldLocalPositions.size(); ++i) {
            Vector_t<double, 3> Eexpected(0.0);
            Vector_t<double, 3> Bexpected(0.0);
            bend.apply(fieldLocalPositions[i], Vector_t<double, 3>(0.0), 0.0, Eexpected, Bexpected);
            const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                    Bexpected, curvature, fieldLocalPositions[i](2), bodyLength);

            EXPECT_NEAR(BhostView(i)(0), BexpectedEntry(0), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(BhostView(i)(1), BexpectedEntry(1), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(BhostView(i)(2), BexpectedEntry(2), 1.0e-12) << "particle " << i;
        }
    }

    TEST_F(BendRepParticleContainerTest,
           FullSliceSequenceBorisUpdateMatchesScalarHostUpdateForDeterministicManyParticleSet) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.02;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 1.0e-11;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(1.2);
        bend.getField().setNormalComponent(1, 0.4);
        bend.getField().setSkewComponent(1, -0.2);
        const double curvature = bendAngle / bodyLength;

        const double entryFringe = bend.getEntryFringeSupportLength();
        const double exitFringe  = bend.getExitFringeSupportLength();
        const std::vector<Vector_t<double, 3>> fieldLocalPositions = {
                {0.0, 0.0, -0.5 * entryFringe},
                {1.0e-3, 0.0, 0.05 * bodyLength},
                {-1.0e-3, 0.0, 0.35 * bodyLength},
                {0.0, 1.0e-3, 0.65 * bodyLength},
                {0.0, -1.0e-3, 0.75 * bodyLength},
                {0.8e-3, -0.6e-3, bodyLength + 0.5 * exitFringe}};
        const Vector_t<double, 3> initialMomentum(1.0e-3, -2.0e-3, 1.2);
        std::vector<std::array<double, 3>> entryPositions;
        for (const auto& fieldLocalPosition : fieldLocalPositions) {
            const Vector_t<double, 3> entryCartesian = fieldLocalToEntryCartesian(
                    fieldLocalPosition, curvature, bodyLength);
            entryPositions.push_back({entryCartesian(0), entryCartesian(1), entryCartesian(2)});
        }
        auto pc = makeContainer(entryPositions);

        auto Phost  = pc->P.getHostMirror();
        auto dthost = pc->dt.getHostMirror();
        for (std::size_t i = 0; i < fieldLocalPositions.size(); ++i) {
            Phost(i)  = initialMomentum;
            dthost(i) = dt;
        }
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            pc->transformBunch(slice.entryToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(slice.entryToSliceLocal.inverted());
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::kickDeterministicSet", massEV, charge);

        pc->switchToUnitlessPositions();
        pushParticleContainer(*pc, "BendRepParticleContainerTest::pushDeterministicSet", 0.0);
        pc->switchOffUnitlessPositions();

        auto RhostView = pc->R.getHostMirror();
        auto Pafter    = pc->P.getHostMirror();
        Kokkos::deep_copy(RhostView, pc->R.getView());
        Kokkos::deep_copy(Pafter, pc->P.getView());

        BorisPusher pusher;
        for (std::size_t i = 0; i < fieldLocalPositions.size(); ++i) {
            Vector_t<double, 3> Eexpected(0.0);
            Vector_t<double, 3> Bexpected(0.0);
            bend.apply(fieldLocalPositions[i], initialMomentum, 0.0, Eexpected, Bexpected);
            const Vector_t<double, 3> BexpectedEntry = fieldLocalVectorToEntryCartesian(
                    Bexpected, curvature, fieldLocalPositions[i](2), bodyLength);
            const Vector_t<double, 3> initialPositionEntry(
                    entryPositions[i][0], entryPositions[i][1], entryPositions[i][2]);

            Vector_t<double, 3> pExpected = initialMomentum;
            pusher.kick(
                    Vector_t<double, 3>(0.0), pExpected, Eexpected, BexpectedEntry, dt, massEV,
                    charge);

            Vector_t<double, 3> rExpected = initialPositionEntry * (1.0 / (Physics::c * dt));
            pusher.push(rExpected, pExpected, 0.0);
            rExpected *= Physics::c * dt;

            EXPECT_NEAR(Pafter(i)(0), pExpected(0), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(Pafter(i)(1), pExpected(1), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(Pafter(i)(2), pExpected(2), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(RhostView(i)(0), rExpected(0), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(RhostView(i)(1), rExpected(1), 1.0e-12) << "particle " << i;
            EXPECT_NEAR(RhostView(i)(2), rExpected(2), 1.0e-12) << "particle " << i;
        }
    }

    TEST_F(BendRepParticleContainerTest,
           ProtonLikeEntryFringeBorisKickProducesNegativeVerticalMomentumChange) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 5.0e-12;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;
        constexpr double betaGamma      = 1.2857059679563454;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        const double curvature = bendAngle / bodyLength;
        const double entryFringe = bend.getEntryFringeSupportLength();
        const Vector_t<double, 3> fieldLocalPosition(0.0, 1.0e-3, -0.5 * entryFringe);
        const Vector_t<double, 3> entryCartesian =
                fieldLocalToEntryCartesian(fieldLocalPosition, curvature, bodyLength);
        const Vector_t<double, 3> initialMomentum(0.0, 0.0, betaGamma);

        auto pc     = makeContainer({entryCartesian(0), entryCartesian(1), entryCartesian(2)});
        auto Phost  = pc->P.getHostMirror();
        auto dthost = pc->dt.getHostMirror();
        Phost(0)    = initialMomentum;
        dthost(0)   = dt;
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            pc->transformBunch(slice.entryToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(slice.entryToSliceLocal.inverted());
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::protonLikeEntryFringeKick", massEV, charge);

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_LT(Pafter(0)(1), 0.0);
    }

    TEST_F(BendRepParticleContainerTest,
           RuntimeNormalizedEntryFringeSliceKickProducesPositiveVerticalMomentumChange) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = Physics::pi / 4.0;
        constexpr double fringeHalfGap   = 0.015;
        constexpr double fringeIntegral  = 0.5;
        constexpr std::size_t nSlices    = 64;
        constexpr double dt              = 5.0e-12;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr double betaGamma       = 1.2857059679563454;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);

        BMultipoleField normalizedField;
        normalizedField.setNormalComponent(0, bendAngle / bend.getEffectiveFieldLength());
        bend.setNormalizedField(normalizedField);
        const double momentumEV = betaGammaFromKineticEnergy(kineticEnergyEV, massEV) * massEV;
        bend.updatePhysicalFieldFromMomentumEV(momentumEV, charge);

        const double entryFringe = bend.getEntryFringeSupportLength();
        const Vector_t<double, 3> fieldLocalPosition(0.0, 1.0e-3, -0.5 * entryFringe);
        const Vector_t<double, 3> entryCartesian = fieldLocalToEntryCartesian(
                fieldLocalPosition, bendAngle / bodyLength, bodyLength);
        const Vector_t<double, 3> initialMomentum(0.0, 0.0, betaGamma);

        auto pc     = makeContainer({entryCartesian(0), entryCartesian(1), entryCartesian(2)});
        auto Phost  = pc->P.getHostMirror();
        auto dthost = pc->dt.getHostMirror();
        Phost(0)    = initialMomentum;
        dthost(0)   = dt;
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            pc->transformBunch(slice.entryToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(slice.entryToSliceLocal.inverted());
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::runtimeNormalizedEntryFringeKick", massEV,
                charge);

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_GT(Pafter(0)(1), 0.0);
    }

    TEST_F(BendRepParticleContainerTest,
           BenchmarkStyleLineStartFrameKickProducesNegativeVerticalMomentumChange) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 5.0e-12;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;
        constexpr double betaGamma      = 1.2857059679563454;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        const double entryFringe = bend.getEntryFringeSupportLength();
        auto pc                  = makeContainer({0.0, 1.0e-3, 0.0});
        auto Phost               = pc->P.getHostMirror();
        auto dthost              = pc->dt.getHostMirror();
        Phost(0)                 = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        dthost(0)                = dt;
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const CoordinateSystemTrafo lineStartToEntryLocal(
                Vector_t<double, 3>(0.0, 0.0, entryFringe), Quaternion());
        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * lineStartToEntryLocal;
            CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
            pc->transformBunch(refToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(sliceLocalToRef);
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkStyleLineStartKick", massEV, charge);

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_LT(Pafter(0)(1), 0.0);
    }

    TEST_F(BendRepParticleContainerTest,
           BenchmarkStyleHalfStepKickProducesNegativeVerticalMomentumChange) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 5.0e-12;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;
        constexpr double betaGamma      = 1.2857059679563454;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        const double entryFringe = bend.getEntryFringeSupportLength();
        const double beta        = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double halfStepZ   = 0.5 * beta * Physics::c * dt;
        auto pc                  = makeContainer({0.0, 1.0e-3, halfStepZ});
        auto Phost               = pc->P.getHostMirror();
        auto dthost              = pc->dt.getHostMirror();
        Phost(0)                 = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        dthost(0)                = dt;
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const CoordinateSystemTrafo lineStartToEntryLocal(
                Vector_t<double, 3>(0.0, 0.0, entryFringe), Quaternion());
        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * lineStartToEntryLocal;
            CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
            pc->transformBunch(refToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(sliceLocalToRef);
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkStyleHalfStepKick", massEV, charge);

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_LT(Pafter(0)(1), 0.0);
    }

    TEST_F(BendRepParticleContainerTest,
           BenchmarkExplicitPlacementHalfStepKickUsesNegativeVerticalMomentumChange) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 5.0e-12;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;
        constexpr double betaGamma      = 1.2857059679563454;
        constexpr double theta          = -Physics::pi / 4.0;
        constexpr double bodyX          = -2.760036392e-01;
        constexpr double bodyZNominal   = 4.130686370e-01;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        Quaternion rotTheta(std::cos(0.5 * theta), 0.0, std::sin(0.5 * theta), 0.0);
        CoordinateSystemTrafo globalToLocal(
                Vector_t<double, 3>(bodyX, 0.0, bodyZNominal + bend.getEntryFringeSupportLength()),
                rotTheta.conjugate());
        bend.setCSTrafoGlobal2Local(globalToLocal);
        bend.fixPosition();

        const double beta      = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double halfStepZ = 0.5 * beta * Physics::c * dt;
        auto pc                = makeContainer({0.0, 1.0e-3, halfStepZ});
        auto Phost             = pc->P.getHostMirror();
        auto dthost            = pc->dt.getHostMirror();
        Phost(0)               = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        dthost(0)              = dt;
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        const CoordinateSystemTrafo lineStartToEntryLocal =
                bend.getPlacedElement().getNominalEntryTransform();
        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * lineStartToEntryLocal;
            CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
            pc->transformBunch(refToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(sliceLocalToRef);
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkExplicitPlacementHalfStepKick",
                massEV, charge);

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_LT(Pafter(0)(1), 0.0);
    }

    TEST_F(BendRepParticleContainerTest,
           BenchmarkExplicitPlacementFullPicStepUsesNegativeVerticalMomentumChange) {
        constexpr double bodyLength     = 1.0;
        constexpr double bendAngle      = Physics::pi / 4.0;
        constexpr double fringeHalfGap  = 0.015;
        constexpr double fringeIntegral = 0.5;
        constexpr std::size_t nSlices   = 64;
        constexpr double dt             = 5.0e-12;
        constexpr double massEV         = Physics::m_p * Units::GeV2eV;
        constexpr double charge         = 1.0;
        constexpr double betaGamma      = 1.2857059679563454;
        constexpr double theta          = -Physics::pi / 4.0;
        constexpr double bodyX          = -2.760036392e-01;
        constexpr double bodyZNominal   = 4.130686370e-01;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setElementLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(Physics::pi / 8.0);
        bend.setExitAngle(Physics::pi / 8.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);
        bend.setNSlices(nSlices);
        bend.setB(-1.2);

        Quaternion rotTheta(std::cos(0.5 * theta), 0.0, std::sin(0.5 * theta), 0.0);
        CoordinateSystemTrafo globalToLocal(
                Vector_t<double, 3>(bodyX, 0.0, bodyZNominal + bend.getEntryFringeSupportLength()),
                rotTheta.conjugate());
        bend.setCSTrafoGlobal2Local(globalToLocal);
        bend.fixPosition();

        auto pc = makeContainer(
                std::vector<std::array<double, 3>>{{0.0, 0.0, 0.0}, {0.0, 1.0e-3, 0.0}});
        auto Phost  = pc->P.getHostMirror();
        auto dthost = pc->dt.getHostMirror();
        for (std::size_t i = 0; i < 2; ++i) {
            Phost(i)  = Vector_t<double, 3>(0.0, 0.0, betaGamma);
            dthost(i) = dt;
        }
        Kokkos::deep_copy(pc->P.getView(), Phost);
        Kokkos::deep_copy(pc->dt.getView(), dthost);
        Kokkos::fence();

        pc->switchToUnitlessPositions();
        pushParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkExplicitPlacementHalfDrift", 0.0);
        pc->switchOffUnitlessPositions();

        const CoordinateSystemTrafo lineStartToEntryLocal =
                bend.getPlacedElement().getNominalEntryTransform();
        const auto slices = bend.buildTrackingSlices();
        for (const auto& slice : slices) {
            CoordinateSystemTrafo refToSliceLocal = slice.entryToSliceLocal * lineStartToEntryLocal;
            CoordinateSystemTrafo sliceLocalToRef = refToSliceLocal.inverted();
            pc->transformBunch(refToSliceLocal);
            bend.applySlice(pc, slice);
            pc->transformBunch(sliceLocalToRef);
        }

        kickParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkExplicitPlacementKick", massEV,
                charge);

        pc->switchToUnitlessPositions();
        pushParticleContainer(
                *pc, "BendRepParticleContainerTest::benchmarkExplicitPlacementSecondHalfDrift",
                0.0);
        pc->switchOffUnitlessPositions();

        auto Pafter = pc->P.getHostMirror();
        Kokkos::deep_copy(Pafter, pc->P.getView());

        EXPECT_NEAR(Pafter(0)(1), 0.0, 1.0e-15);
        EXPECT_LT(Pafter(1)(1), 0.0);
    }

    TEST(SBendRep, ProtonOnAxisIntegratesToNegativeFortyFiveDegreeDeflection) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double arcLength       = 1.0;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(arcLength, angle / arcLength);
        bend.setElementLength(arcLength);
        bend.setBendAngle(angle);

        const double radius    = arcLength / angle;
        const double betaGamma = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by        = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(by, 0.0);

        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = arcLength / (beta * Physics::c);
        const ReferenceState state   = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps);

        const double deflection = std::atan2(state.p(0), state.p(2));
        EXPECT_NEAR(deflection, -angle, 5.0e-5);
        EXPECT_NEAR(state.r(0), -std::copysign(radius * (1.0 - std::cos(angle)), angle), 1.0e-6);
        EXPECT_NEAR(state.r(2), radius * std::sin(std::abs(angle)), 1.0e-6);
    }

    TEST(SBendRep, RuntimeNormalizedElectronOnAxisIntegratesToPositiveFortyFiveDegreeDeflection) {
        constexpr double angle           = -Physics::pi / 4.0;
        constexpr double arcLength       = 1.0;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;

        SBendRep bend("SBENDNEG");
        bend.getGeometry() = PlanarArcGeometry(arcLength, angle / arcLength);
        bend.setElementLength(arcLength);
        bend.setBendAngle(angle);

        const double radius    = std::abs(arcLength / angle);
        const double betaGamma = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        BMultipoleField normalizedField;
        normalizedField.setNormalComponent(0, angle / arcLength);
        bend.setNormalizedField(normalizedField);
        bend.updatePhysicalFieldFromMomentumEV(betaGamma * massEV, charge);

        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = arcLength / (beta * Physics::c);
        const ReferenceState state   = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps);

        const double deflection = std::atan2(state.p(0), state.p(2));
        EXPECT_NEAR(deflection, -angle, 5.0e-5);
        EXPECT_NEAR(state.r(0), -std::copysign(radius * (1.0 - std::cos(angle)), angle), 1.0e-6);
        EXPECT_NEAR(state.r(2), radius * std::sin(std::abs(angle)), 1.0e-6);
    }

    TEST(SBendRep, RuntimeNormalizationUsesBunchReferenceMomentumWhenDesignEnergyIsUnset) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double length          = 1.0;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(length, angle / length);
        bend.setElementLength(length);
        bend.setBendAngle(angle);
        bend.setFieldAmplitude(angle / length, 0.0);

        BMultipoleField normalizedField;
        normalizedField.setNormalComponent(0, angle / length);
        bend.setNormalizedField(normalizedField);

        const double momentumEV = betaGammaFromKineticEnergy(kineticEnergyEV, massEV) * massEV;
        bend.updatePhysicalFieldFromMomentumEV(momentumEV, 1.0);

        const double expectedBy = momentumEV / Physics::c * (angle / length);
        EXPECT_NEAR(bend.getB(), expectedBy, 1.0e-9 * std::abs(expectedBy));
        EXPECT_GT(bend.getB(), 0.0);
    }

    TEST(SBendRep, RuntimeNormalizedProtonEntryFringeProducesExpectedVerticalKickSign) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double entryAngle      = 0.15;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double yOffset         = 1.0e-4;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(entryAngle);
        bend.setExitAngle(0.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        BMultipoleField normalizedField;
        normalizedField.setNormalComponent(0, bendAngle / bend.getEffectiveFieldLength());
        bend.setNormalizedField(normalizedField);
        const double momentumEV = betaGammaFromKineticEnergy(kineticEnergyEV, massEV) * massEV;
        bend.updatePhysicalFieldFromMomentumEV(momentumEV, charge);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta = momentumEV / std::sqrt(momentumEV * momentumEV + massEV * massEV);
        const double integrationTime = (-fieldBegin) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, fieldBegin);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, fieldBegin);
        yState.p = Vector_t<double, 3>(0.0, 0.0, momentumEV / massEV);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / bend.getEffectiveFieldLength();
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(entryAngle) * std::sin(entryAngle))
                           / std::cos(entryAngle);
        const double expectedVertical = signedCurvature * std::tan(entryAngle - psi);
        const double deltaYP          = (yFocused.p(1) - base.p(1)) / (momentumEV / massEV);

        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
    }

    TEST(RBendRep, RuntimeNormalizationUsesReferenceChargeSign) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double chordLength     = 1.0;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(chordLength);
        bend.getGeometry().setBendAngle(angle);
        bend.setElementLength(chordLength);
        bend.setBendAngle(angle);
        bend.setFieldAmplitude(angle / chordLength, 0.0);

        BMultipoleField normalizedField;
        normalizedField.setNormalComponent(0, angle / chordLength);
        bend.setNormalizedField(normalizedField);

        const double momentumEV = betaGammaFromKineticEnergy(kineticEnergyEV, massEV) * massEV;
        bend.updatePhysicalFieldFromMomentumEV(momentumEV, charge);

        const double expectedBy = momentumEV / (charge * Physics::c) * (angle / chordLength);
        EXPECT_NEAR(bend.getB(), expectedBy, 1.0e-9 * std::abs(expectedBy));
        EXPECT_LT(bend.getB(), 0.0);
    }

    TEST(SBendRep, ReferenceMomentumResolutionPrefersDesignEnergyWhenProvided) {
        constexpr double bunchKineticEnergyEV   = 300.0 * Units::MeV2eV;
        constexpr double designKineticEnergyMeV = 590.0;
        constexpr double designKineticEnergyEV  = designKineticEnergyMeV * Units::MeV2eV;
        constexpr double massEV                 = Physics::m_p * Units::GeV2eV;

        SBendRep bend("SBEND");
        bend.setDesignEnergy(designKineticEnergyMeV, false);

        const double bunchMomentumEV =
                betaGammaFromKineticEnergy(bunchKineticEnergyEV, massEV) * massEV;
        const double designMomentumEV =
                betaGammaFromKineticEnergy(designKineticEnergyEV, massEV) * massEV;
        EXPECT_NEAR(
                bend.resolveReferenceMomentumEV(bunchMomentumEV, massEV), designMomentumEV,
                1.0e-12 * designMomentumEV);
    }

    TEST(RBendRep, ElectronOnAxisIntegratesToNegativeFortyFiveDegreeDeflection) {
        constexpr double angle           = Physics::pi / 4.0;
        constexpr double chordLength     = 1.0;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(chordLength);
        bend.getGeometry().setBendAngle(angle);
        bend.setElementLength(chordLength);
        bend.setBendAngle(angle);

        const double radius    = chordLength / (2.0 * std::sin(angle / 2.0));
        const double arcLength = radius * angle;
        const double betaGamma = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by        = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(by, 0.0);

        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = arcLength / (beta * Physics::c);
        const ReferenceState state   = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps);

        const double deflection = std::atan2(state.p(0), state.p(2));
        EXPECT_NEAR(deflection, -angle, 5.0e-5);
        EXPECT_NEAR(state.r(0), -radius * (1.0 - std::cos(angle)), 1.0e-6);
        EXPECT_NEAR(state.r(2), radius * std::sin(angle), 1.0e-6);
    }

    TEST(SBendRep, FringeSupportProducesExpectedOpalEngeProfile) {
        constexpr double bodyLength = 1.0;
        constexpr double angle      = Physics::pi / 4.0;
        constexpr double halfGap    = 0.02;
        const double edgeScale      = 1.0 / (1.0 + std::exp(0.478959));
        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, angle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(angle);
        bend.setFringeHalfGap(halfGap);
        bend.setFringeIntegral(0.5);
        bend.setB(-1.0);
        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);

        Vector_t<double, 3> E(0.0);
        Vector_t<double, 3> B(0.0);
        EXPECT_NEAR(fieldBegin, -10.0 * halfGap, 1.0e-12);
        EXPECT_NEAR(fieldEnd, bodyLength + 10.0 * halfGap, 1.0e-12);
        EXPECT_NEAR(bend.getEffectiveFieldLength(), bodyLength, 2.0e-6);

        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, -0.21), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), 0.0, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, 0.0), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -edgeScale, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, 0.5), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -1.0, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, bodyLength),
                Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), -edgeScale, 1.0e-12);

        B = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, 1.21),
                Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(1), 0.0, 1.0e-12);

        const double expAtEdge           = std::exp(0.478959);
        const double profileGap          = 2.0 * halfGap;
        const double entryScaleDerivative = 1.911289 * expAtEdge
                                           / (profileGap * std::pow(1.0 + expAtEdge, 2.0));
        const double yOffset = 1.0e-3;
        B                    = Vector_t<double, 3>(0.0);
        EXPECT_FALSE(bend.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, yOffset, 0.0), Vector_t<double, 3>(0.0), 0.0, E, B));
        EXPECT_NEAR(B(2), -entryScaleDerivative * yOffset, 1.0e-12);
    }

    TEST(SBendRep, EntryFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double entryAngle      = 0.15;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(entryAngle);
        bend.setExitAngle(0.0);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (-fieldBegin) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, fieldBegin);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, fieldBegin);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, fieldBegin);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(entryAngle) * std::sin(entryAngle))
                           / std::cos(entryAngle);
        const double expectedHorizontal = -signedCurvature * std::tan(entryAngle);
        const double expectedVertical   = signedCurvature * std::tan(entryAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(SBendRep, ExitFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double exitAngle       = 0.14;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 590.0 * Units::MeV2eV;
        constexpr double massEV          = Physics::m_p * Units::GeV2eV;
        constexpr double charge          = 1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(bodyLength, bendAngle / bodyLength);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(0.0);
        bend.setExitAngle(exitAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (fieldEnd - bodyLength) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, bodyLength);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, bodyLength);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, bodyLength);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(exitAngle) * std::sin(exitAngle))
                           / std::cos(exitAngle);
        const double expectedHorizontal = -signedCurvature * std::tan(exitAngle);
        const double expectedVertical   = signedCurvature * std::tan(exitAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(RBendRep, EntryFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double entryAngle      = 0.12;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(entryAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (-fieldBegin) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, fieldBegin);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, fieldBegin);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, fieldBegin);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(entryAngle) * std::sin(entryAngle))
                           / std::cos(entryAngle);
        const double expectedHorizontal = -signedCurvature * std::tan(entryAngle);
        const double expectedVertical   = signedCurvature * std::tan(entryAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

    TEST(RBendRep, ExitFringeProducesEdgeFocusingKick) {
        constexpr double bodyLength      = 1.0;
        constexpr double bendAngle       = 0.2;
        constexpr double exitAngle       = 0.11;
        constexpr double fringeHalfGap   = 0.02;
        constexpr double fringeIntegral  = 0.5;
        constexpr double kineticEnergyEV = 1.0 * Units::GeV2eV;
        constexpr double massEV          = Physics::m_e * Units::GeV2eV;
        constexpr double charge          = -1.0;
        constexpr std::size_t steps      = 20000;
        constexpr double xOffset         = 1.0e-4;
        constexpr double yOffset         = 1.0e-4;

        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(bodyLength);
        bend.getGeometry().setBendAngle(bendAngle);
        bend.setLength(bodyLength);
        bend.setBendAngle(bendAngle);
        bend.setEntranceAngle(bendAngle - exitAngle);
        bend.setFringeHalfGap(fringeHalfGap);
        bend.setFringeIntegral(fringeIntegral);

        const double effectiveLength = bend.getEffectiveFieldLength();
        const double radius          = effectiveLength / bendAngle;
        const double betaGamma       = betaGammaFromKineticEnergy(kineticEnergyEV, massEV);
        const double by              = -dipoleFieldMagnitude(betaGamma, massEV, charge, radius);
        bend.setB(by);
        bend.setFieldAmplitude(bendAngle / effectiveLength, 0.0);

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        bend.getFieldExtend(fieldBegin, fieldEnd);
        const double beta            = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
        const double integrationTime = (fieldEnd - bodyLength) / (beta * Physics::c);

        const ReferenceState base = integrateReferenceParticle(
                bend, massEV, charge, kineticEnergyEV, integrationTime, steps, bodyLength);

        ReferenceState xState{};
        xState.r = Vector_t<double, 3>(xOffset, 0.0, bodyLength);
        xState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState xFocused =
                integrateReferenceParticle(bend, massEV, charge, xState, integrationTime, steps);

        ReferenceState yState{};
        yState.r = Vector_t<double, 3>(0.0, yOffset, bodyLength);
        yState.p = Vector_t<double, 3>(0.0, 0.0, betaGamma);
        const ReferenceState yFocused =
                integrateReferenceParticle(bend, massEV, charge, yState, integrationTime, steps);

        const double signedCurvature = bendAngle / effectiveLength;
        const double psi             = signedCurvature * fringeHalfGap * fringeIntegral
                           * (1.0 + std::sin(exitAngle) * std::sin(exitAngle))
                           / std::cos(exitAngle);
        const double expectedHorizontal = -signedCurvature * std::tan(exitAngle);
        const double expectedVertical   = signedCurvature * std::tan(exitAngle - psi);
        const double deltaXP            = (xFocused.p(0) - base.p(0)) / betaGamma;
        const double deltaYP            = (yFocused.p(1) - base.p(1)) / betaGamma;

        EXPECT_NEAR(deltaXP / xOffset, expectedHorizontal, 5.0e-4);
        EXPECT_NEAR(deltaYP / yOffset, expectedVertical, 5.0e-4);
        EXPECT_LT(std::abs(expectedVertical), std::abs(expectedHorizontal));
    }

}  // namespace
