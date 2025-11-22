#include <gtest/gtest.h>
#include <memory>
#include <Kokkos_Core.hpp>

#include "Gaussian.h"
#include "Distribution.h"
#include "SamplingBase.hpp"

// -----------------------------------------------------------------------------
// Minimal stubs to make Gaussian work inside the test environment
// -----------------------------------------------------------------------------

class TestDistribution : public Distribution {
public:
    TestDistribution() {
        sigmaR_ = {1.0, 1.0, 1.0};
        sigmaP_ = {0.1, 0.1, 0.1};
        cutoffR_ = {3.0, 3.0, 3.0};
        avrgpz_ = 5.0;
    }

    Vector_t<double,3> getSigmaR() const override { return sigmaR_; }
    Vector_t<double,3> getSigmaP() const override { return sigmaP_; }
    Vector_t<double,3> getCutoffR() const override { return cutoffR_; }
    double getAvrgpz() const override { return avrgpz_; }

private:
    Vector_t<double,3> sigmaR_;
    Vector_t<double,3> sigmaP_;
    Vector_t<double,3> cutoffR_;
    double avrgpz_;
};


// -----------------------------------------------------------------------------
// A very small dummy particle container that mimics IPPL storage
// -----------------------------------------------------------------------------

template<typename T, int Dim>
class DummyParticles {
public:
    using view_type = Kokkos::View<T*[Dim]>;

    view_type getView() { return view_; }
    size_t size() const { return n_; }

    void create(size_t n) {
        n_ = n;
        view_ = view_type("particles", n);
    }

private:
    size_t n_ = 0;
    view_type view_;
};

template<typename T, int Dim>
class DummyParticleContainer {
public:
    DummyParticles<T,Dim> R;
    DummyParticles<T,Dim> P;
};


// -----------------------------------------------------------------------------
// Test fixture
// -----------------------------------------------------------------------------

class GaussianSamplerTest : public ::testing::Test {
protected:
    void SetUp() override {
        Kokkos::initialize();

        pc = std::make_shared<ParticleContainer_t>();
        fc = std::make_shared<FieldContainer_t>();   // not used, but required
        dist = std::make_shared<TestDistribution>();

        sampler = std::make_unique<Gaussian>(pc, fc, dist);
    }

    void TearDown() override {
        Kokkos::finalize();
    }

    std::shared_ptr<ParticleContainer_t> pc;
    std::shared_ptr<FieldContainer_t> fc;
    std::shared_ptr<Distribution_t> dist;

    std::unique_ptr<Gaussian> sampler;
};


// -----------------------------------------------------------------------------
// Actual test
// -----------------------------------------------------------------------------

TEST_F(GaussianSamplerTest, GeneratesCorrectParticleCountAndMomentumShift)
{
    size_t N = 1000;
    Vector_t<double,3> nr = {1.0, 1.0, 1.0};

    sampler->generateParticles(N, nr);

    // ---- Verify particle count ----
    auto Rv = pc->R.getView();
    auto Pv = pc->P.getView();
    ASSERT_EQ(Rv.extent(0), N);
    ASSERT_EQ(Pv.extent(0), N);

    // ---- Compute means on host ----
    double meanRx = 0, meanRy = 0, meanRz = 0;
    double meanPz = 0;

    auto R_h = Kokkos::create_mirror_view(Rv);
    auto P_h = Kokkos::create_mirror_view(Pv);

    Kokkos::deep_copy(R_h, Rv);
    Kokkos::deep_copy(P_h, Pv);

    for (size_t i = 0; i < N; i++) {
        meanRx += R_h(i,0);
        meanRy += R_h(i,1);
        meanRz += R_h(i,2);
        meanPz += P_h(i,2);
    }

    meanRx /= N;
    meanRy /= N;
    meanRz /= N;
    meanPz /= N;

    // ---- Check that R mean ≈ 0 after the correction step ----
    EXPECT_NEAR(meanRx, 0.0, 1e-2);
    EXPECT_NEAR(meanRy, 0.0, 1e-2);
    EXPECT_NEAR(meanRz, 0.0, 1e-2);

    // ---- Check that average pz ≈ avrgpz ----
    EXPECT_NEAR(meanPz, dist->getAvrgpz(), 1e-2);
}
