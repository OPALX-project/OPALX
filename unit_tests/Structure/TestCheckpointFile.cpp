#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/CheckpointFile.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace {
    class TestFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& type) {
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], type);
        }
        void setBC(const std::string& bc) {
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::BCFFTX], bc);
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::BCFFTY], bc);
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::BCFFTZ], bc);
        }
    };

    using Bunch = PartBunch<double, 3>;

    void expectVectorEqual(const Vector_t<double, 3>& actual, const Vector_t<double, 3>& expected) {
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(actual[d], expected[d]);
        }
    }

    class CheckpointFileTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            gmsg = new Inform(nullptr, -1);
            OpalData::getInstance()->storeInputFn("checkpoint_unit.in");
            Options::enableHDF5 = false;
        }

        static void TearDownTestSuite() {
            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0) {
                std::remove("checkpoint_unit.stat");
                std::remove("checkpoint_unit.lbal");
                std::remove("checkpoint_unit_checkpoint.h5");
                std::remove("checkpoint_unit_checkpoint.h5.tmp");
                std::remove("checkpoint_unit_attributes.h5");
                std::remove("checkpoint_unit_attributes.h5.tmp");
            }
            ippl::Comm->barrier();
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            Options::useQMAttributes = false;
            fieldSolver              = std::make_unique<TestFieldSolverCmd>();
            fieldSolver->setType("NONE");
            fieldSolver->setNX(8);
            fieldSolver->setNY(8);
            fieldSolver->setNZ(8);
            fieldSolver->setBC("PERIODIC");

            beamObject = std::make_unique<Beam>();
            beam       = Beam::find("UNNAMED_BEAM");
            ASSERT_NE(beam, nullptr);
        }

        void TearDown() override {
            Options::useQMAttributes = false;
            if (beam) {
                Attributes::setRealArray(*beam->findAttribute("POLARIZATION"), {});
            }
        }

        std::unique_ptr<Bunch> makeBunch(DataSink& sink) {
            return std::make_unique<Bunch>(
                    std::vector<double>{1.25e-12, -2.5e-12}, std::vector<double>{0.511e-3, 0.938},
                    std::vector<Beam*>{beam, beam}, std::vector<std::size_t>{16, 16}, 1.0, "LF2",
                    fieldSolver.get(), &sink);
        }

        void populate(Bunch& bunch) {
            const int rank = ippl::Comm->rank();
            for (std::size_t ci = 0; ci < bunch.getNumParticleContainers(); ++ci) {
                auto pc = bunch.getParticleContainer(ci);
                pc->createParticles(4);
                auto rHost   = pc->R.getHostMirror();
                auto pHost   = pc->P.getHostMirror();
                auto dtHost  = pc->dt.getHostMirror();
                auto binHost = pc->Bin.getHostMirror();
                for (std::size_t i = 0; i < 4; ++i) {
                    const double base = 100.0 * rank + 10.0 * ci + i;
                    rHost(i)          = Vector_t<double, 3>(base + 0.1, base + 0.2, base + 0.3);
                    pHost(i)          = Vector_t<double, 3>(base + 1.1, base + 1.2, base + 1.3);
                    dtHost(i)         = 2.0e-12 + i * 1.0e-15;
                    binHost(i)        = static_cast<Bunch::binIndex_t>(i + ci);
                }
                Kokkos::deep_copy(pc->R.getView(), rHost);
                Kokkos::deep_copy(pc->P.getView(), pHost);
                Kokkos::deep_copy(pc->dt.getView(), dtHost);
                Kokkos::deep_copy(pc->Bin.getView(), binHost);
                pc->set_sPos(0.5 + ci);
                pc->setRefPartR(Vector_t<double, 3>(ci + 0.1, ci + 0.2, ci + 0.3));
                pc->setRefPartP(Vector_t<double, 3>(ci + 1.1, ci + 1.2, ci + 1.3));
                pc->setToLabTrafo(CoordinateSystemTrafo(
                        Vector_t<double, 3>(ci + 2.1, ci + 2.2, ci + 2.3),
                        Quaternion(1.0, 0.0, 0.0, 0.0)));
            }
            bunch.setT(7.5e-9);
            bunch.setdT(2.0e-12);
            bunch.setGlobalTrackStep(37);
            bunch.resetPcActive();
            bunch.setPcAtSStop(1);
            Kokkos::fence();
        }

        std::unique_ptr<TestFieldSolverCmd> fieldSolver;
        std::unique_ptr<Beam> beamObject;
        Beam* beam = nullptr;
    };

    TEST_F(CheckpointFileTest, RoundTripRestoresMultiContainerStateWithHdf5DiagnosticsDisabled) {
        DataSink sourceSink;
        auto source = makeBunch(sourceSink);
        populate(*source);

        struct SavedLocal {
            std::vector<Vector_t<double, 3>> r;
            std::vector<Vector_t<double, 3>> p;
            std::vector<double> dt;
            std::vector<long long> id;
            std::vector<Bunch::binIndex_t> bin;
        };
        std::vector<SavedLocal> saved(source->getNumParticleContainers());
        for (std::size_t ci = 0; ci < saved.size(); ++ci) {
            auto pc = source->getParticleContainer(ci);
            const auto r =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            const auto p =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            const auto d =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
            const auto ids =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
            const auto bins =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->Bin.getView());
            for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
                saved[ci].r.push_back(r(i));
                saved[ci].p.push_back(p(i));
                saved[ci].dt.push_back(d(i));
                saved[ci].id.push_back(ids(i));
                saved[ci].bin.push_back(bins(i));
            }
        }

        const std::string path = CheckpointFile::defaultPath("checkpoint_unit");
        OpalData::getInstance()->setMaxPhase("TEST_CAVITY", 0.375);
        const int phasesBeforeRead = OpalData::getInstance()->getNumberOfMaxPhases();
        ASSERT_NO_THROW(CheckpointFile::write(path, *source, 2, 7));
        if (ippl::Comm->rank() == 0) {
            EXPECT_TRUE(std::filesystem::exists(path));
            EXPECT_FALSE(std::filesystem::exists(path + ".tmp"));
        }

        const auto inspected = CheckpointFile::inspect(path);
        EXPECT_EQ(inspected.globalTrackStep, 37);
        EXPECT_DOUBLE_EQ(inspected.time, 7.5e-9);
        EXPECT_DOUBLE_EQ(inspected.dt, 2.0e-12);
        EXPECT_EQ(inspected.numContainers, 2u);
        EXPECT_EQ(inspected.stepSizeSegment, 2u);
        EXPECT_EQ(inspected.stepsCompletedInSegment, 7u);

        DataSink targetSink;
        auto target         = makeBunch(targetSink);
        const auto restored = CheckpointFile::read(path, *target);
        ASSERT_EQ(OpalData::getInstance()->getNumberOfMaxPhases(), phasesBeforeRead + 1);
        auto lastPhase = OpalData::getInstance()->getLastMaxPhases();
        --lastPhase;
        EXPECT_EQ(lastPhase->first, "TEST_CAVITY");
        EXPECT_DOUBLE_EQ(lastPhase->second, 0.375);
        EXPECT_EQ(restored.globalTrackStep, 37);
        EXPECT_EQ(target->getGlobalTrackStep(), 37);
        EXPECT_DOUBLE_EQ(target->getT(), 7.5e-9);
        EXPECT_DOUBLE_EQ(target->getdT(), 2.0e-12);
        EXPECT_FALSE(target->pcAtSStop(0));
        EXPECT_TRUE(target->pcAtSStop(1));

        for (std::size_t ci = 0; ci < saved.size(); ++ci) {
            auto pc = target->getParticleContainer(ci);
            ASSERT_EQ(pc->getLocalNum(), saved[ci].r.size());
            EXPECT_DOUBLE_EQ(pc->get_sPos(), 0.5 + ci);
            EXPECT_DOUBLE_EQ(pc->getChargePerParticle(), ci == 0 ? 1.25e-12 : -2.5e-12);
            EXPECT_DOUBLE_EQ(pc->getMassPerParticle(), ci == 0 ? 0.511e-3 : 0.938);
            expectVectorEqual(pc->getRefPartR(), Vector_t<double, 3>(ci + 0.1, ci + 0.2, ci + 0.3));
            expectVectorEqual(pc->getRefPartP(), Vector_t<double, 3>(ci + 1.1, ci + 1.2, ci + 1.3));

            const auto r =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            const auto p =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
            const auto d =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
            const auto ids =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
            const auto bins =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->Bin.getView());
            for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
                expectVectorEqual(r(i), saved[ci].r[i]);
                expectVectorEqual(p(i), saved[ci].p[i]);
                EXPECT_DOUBLE_EQ(d(i), saved[ci].dt[i]);
                EXPECT_EQ(ids(i), saved[ci].id[i]);
                EXPECT_EQ(bins(i), saved[ci].bin[i]);
            }
        }
    }

    TEST_F(CheckpointFileTest, RoundTripRestoresPerParticleChargeMassAndSpin) {
        Options::useQMAttributes = true;
        Attributes::setRealArray(*beam->findAttribute("POLARIZATION"), {0.0, 0.0, 1.0});

        DataSink sourceSink;
        auto source = makeBunch(sourceSink);
        populate(*source);

        for (std::size_t ci = 0; ci < source->getNumParticleContainers(); ++ci) {
            auto pc      = source->getParticleContainer(ci);
            auto qHost   = Kokkos::create_mirror_view(pc->getQView());
            auto mHost   = Kokkos::create_mirror_view(pc->getMView());
            auto polHost = pc->Pol.getHostMirror();
            for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
                qHost(i)   = (ci + 1.0) * 1.0e-12 + i * 1.0e-15;
                mHost(i)   = (ci + 1.0) * 0.1 + i * 1.0e-3;
                polHost(i) = Bunch::ParticleContainer_t::spin_vector_type{
                        static_cast<float>(0.1 * i), static_cast<float>(0.2 * i),
                        static_cast<float>(1.0 - 0.1 * i)};
            }
            Kokkos::deep_copy(pc->getQView(), qHost);
            Kokkos::deep_copy(pc->getMView(), mHost);
            Kokkos::deep_copy(pc->Pol.getView(), polHost);
        }

        const std::string path = "checkpoint_unit_attributes.h5";
        ASSERT_NO_THROW(CheckpointFile::write(path, *source, 0, 5));

        DataSink targetSink;
        auto target = makeBunch(targetSink);
        ASSERT_NO_THROW(CheckpointFile::read(path, *target));

        for (std::size_t ci = 0; ci < target->getNumParticleContainers(); ++ci) {
            auto pc = target->getParticleContainer(ci);
            ASSERT_TRUE(pc->hasSpin());
            const auto q = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
            const auto mass =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getMView());
            const auto pol =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->Pol.getView());
            for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
                EXPECT_DOUBLE_EQ(q(i), (ci + 1.0) * 1.0e-12 + i * 1.0e-15);
                EXPECT_DOUBLE_EQ(mass(i), (ci + 1.0) * 0.1 + i * 1.0e-3);
                EXPECT_FLOAT_EQ(pol(i)[0], static_cast<float>(0.1 * i));
                EXPECT_FLOAT_EQ(pol(i)[1], static_cast<float>(0.2 * i));
                EXPECT_FLOAT_EQ(pol(i)[2], static_cast<float>(1.0 - 0.1 * i));
            }
        }
    }
}  // namespace
