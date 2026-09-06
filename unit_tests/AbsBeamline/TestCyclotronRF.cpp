#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <unistd.h>
#include "Fields/CyclotronRFProfile.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Ippl.h"

class CyclotronRFTest : public ::testing::Test {
protected:
    std::string filename;
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
    void SetUp() override {
        filename = (std::filesystem::temp_directory_path()
                    / ("opalx-rf-" + std::to_string(getpid()) + ".dat")).string();
        std::ofstream out(filename);
        // f(u)=1+u^3, exactly representable by Hermite interpolation.
        out << "3\n0 1 0\n0.5 1.125 0.75\n1 2 3\n";
    }
    void TearDown() override { std::filesystem::remove(filename); }
};

TEST_F(CyclotronRFTest, HermiteValuesDerivativesAndEndpoints) {
    CyclotronRFProfile profile(filename);
    for (double u : {0., .1, .5, .9, 1.}) {
        double f, d;
        ASSERT_TRUE(CyclotronRFProfile::evaluate(profile.host, u, f, d));
        EXPECT_NEAR(f, 1 + u*u*u, 1e-14);
        EXPECT_NEAR(d, 3*u*u, 1e-14);
    }
    double f, d;
    EXPECT_FALSE(CyclotronRFProfile::evaluate(profile.host, -1e-6, f, d));
    EXPECT_EQ(f, 0);
    EXPECT_FALSE(CyclotronRFProfile::evaluate(profile.host, 1.000001, f, d));
}

TEST_F(CyclotronRFTest, RejectsMalformedProfiles) {
    for (const std::string text : {"0", "2\n0 1 0\n", "2\n0 1 0\n0 2 1\n",
                                   "2\n0 1 0\n0.9 2 1\n", "2\n0 1 0\n1 2 1\nextra"}) {
        { std::ofstream out(filename); out << text; }
        EXPECT_THROW(CyclotronRFProfile profile(filename), OpalException);
    }
}

TEST_F(CyclotronRFTest, EnergyGainTransitTimeAndPhase) {
    CyclotronRFProfile profile(filename);
    const double mass = 938272088.16;
    const double gamma = 1 + 72e6/mass;
    const double pz = std::sqrt(gamma*gamma-1);
    for (double width : {0., .3}) {
        for (double phase : {0., Physics::pi, Physics::pi/2}) {
            CyclotronRFKick kick{847000, Physics::two_pi*50.65e6, phase, width, 5.2};
            kick.validate();
            Vector_t<double, 3> p(0, 0, pz);
            ASSERT_TRUE(kick.apply(profile.host, 0., 0., mass, p));
            const double a = kick.omega*width/(2*Physics::c*(pz/gamma));
            const double transit = width == 0 ? 1 : std::sin(a)/a;
            EXPECT_NEAR((std::sqrt(1+dot(p,p))-gamma)*mass,
                        847000*transit*std::cos(phase), 3e-7);
        }
    }
}

TEST_F(CyclotronRFTest, FocusingDeviceParityAndInvalidState) {
    CyclotronRFProfile profile(filename);
    CyclotronRFKick kick{847000, Physics::two_pi*50.65e6, .4, .3, 5.2};
    Vector_t<double, 3> expected(.01, 0, .4);
    ASSERT_TRUE(kick.apply(profile.host, .6, 1e-9, 938272088.16, expected));
    const double mass = 938272088.16;
    const double bg = std::hypot(.01, .4), gamma = std::sqrt(1+bg*bg);
    const double phi = kick.omega*1e-9-kick.phase;
    const double a = kick.omega*.3/(2*Physics::c*(bg/gamma));
    const double g1 = gamma+847000*(1+std::pow(.6,3))*std::sin(a)/a*std::cos(phi)/mass;
    const double angle = -3*.6*.6*847000/5.2*std::sin(phi)*Physics::c
                       /(kick.omega*Physics::two_pi*bg*mass);
    const double tangent = std::sqrt(g1*g1-1-.01*.01);
    EXPECT_NEAR(expected[0], std::cos(angle)*.01+std::sin(angle)*tangent, 1e-14);
    EXPECT_NEAR(expected[2], -std::sin(angle)*.01+std::cos(angle)*tangent, 1e-14);
    Kokkos::View<double*> output("result", 4);
    const auto grid = profile.data;
    Kokkos::parallel_for("rf kick test", 1, KOKKOS_LAMBDA(int) {
        Vector_t<double, 3> p(.01, 0, .4);
        output(3) = kick.apply(grid, .6, 1e-9, 938272088.16, p);
        for (int i = 0; i < 3; ++i) output(i) = p[i];
    });
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
    EXPECT_EQ(host(3), 1);
    for (int i=0; i<3; ++i) EXPECT_NEAR(host(i), expected[i], 1e-14);
    Vector_t<double, 3> unsupported(.01, .001, .4);
    EXPECT_FALSE(kick.apply(profile.host, .6, 0, 938272088.16, unsupported));
    EXPECT_EQ(unsupported[1], .001);
    kick.omega = 0;
    EXPECT_THROW(kick.validate(), OpalException);
}

TEST_F(CyclotronRFTest, ZeroVoltageAndRejectedKickDoNotMutateMomentum) {
    CyclotronRFProfile profile(filename);
    CyclotronRFKick kick{0, Physics::two_pi*50.65e6, 0, .3, 5.2};
    Vector_t<double, 3> p(.01, 0, .4);
    ASSERT_TRUE(kick.apply(profile.host, .5, 0, 938272088.16, p));
    EXPECT_EQ(p[0], .01);
    EXPECT_EQ(p[2], .4);
    kick.voltage = -100e6;
    EXPECT_FALSE(kick.apply(profile.host, .5, 0, 938272088.16, p));
    EXPECT_EQ(p[0], .01);
    EXPECT_EQ(p[2], .4);
}

TEST_F(CyclotronRFTest, GapCloneAndNoContinuousFieldContribution) {
    RFCavityRep cavity("gap");
    cavity.configureCyclotronGap(filename, 1, 2,
        CyclotronRFKick{847000, Physics::two_pi*50.65e6, 0, .3, 1});
    std::unique_ptr<RFCavityRep> clone(static_cast<RFCavityRep*>(cavity.clone()));
    ASSERT_TRUE(clone->isCyclotronGap());
    EXPECT_TRUE(clone->gapSupports(1));
    EXPECT_TRUE(clone->gapSupports(2));
    EXPECT_FALSE(clone->gapSupports(.9));
    Vector_t<double, 3> r(1.5,0,0), p(0,0,.4), e(2), b(3);
    clone->apply(r,p,0,e,b);
    EXPECT_FALSE(clone->isInside(r));
    EXPECT_FALSE(clone->isInsideBody(r));
    for(int i=0;i<3;++i) { EXPECT_EQ(e[i],2); EXPECT_EQ(b[i],3); }
    ASSERT_TRUE(clone->applyGapKick(1.5,0,938272088.16,p));
    EXPECT_GT(p[2], .4);
}
