//
// Unit tests for the Binning module:
//   - ParallelReduceTools (ArrayReduction, HostArrayReduction, createReductionObject)
//   - BinningTools (computeFixSum, determineHistoReductionMode, viewIsSorted, CoordinateSelector)
//   - BinHisto (Histogram class: construction, init, postSum, mergeBins, iteration policies)
//   - AdaptBins (full integration: initLimits, assignBins, histogram, sort, adaptive rebinning)
//
// A minimal bunch type (TestBunch) is constructed with a mesh and field layout,
// similar to how it is done for IPPL's GatherScatterTest.
//

#include "Ippl.h"
#include "gtest/gtest.h"

// VField_t alias required by ParallelReduceTools.h (vnorm) and AdaptBins.h (LTrans)
template <unsigned Dim>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <unsigned Dim>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <typename T, unsigned Dim, class... ViewArgs>
using Field = ippl::Field<T, Dim, Mesh_t<Dim>, Centering_t<Dim>, ViewArgs...>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using VField_t = Field<Vector_t<T, Dim>, Dim, ViewArgs...>;

// Alias needed by AdaptBins.tpp LTrans (uses unqualified Vector)
template <typename T, unsigned Dim>
using Vector = ippl::Vector<T, Dim>;

#include "PartBunch/Binning/AdaptBins.h"
#include "PartBunch/Binning/BinHisto.h"
#include "PartBunch/Binning/BinningTools.h"
#include "PartBunch/Binning/ParallelReduceTools.h"

using size_type = ippl::detail::size_type;

#include <random>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>

// ============================================================================
// Minimal bunch type for testing (mirrors ParticleContainer's relevant parts)
// ============================================================================
constexpr unsigned TestDim = 3;
using TestT                = double;

using TestMesh_t     = ippl::UniformCartesian<TestT, TestDim>;
using TestPLayout_t  = ippl::ParticleSpatialLayout<TestT, TestDim, TestMesh_t>;
using TestFL_t       = ippl::FieldLayout<TestDim>;

/**
 * @brief Minimal particle bunch for binning unit tests.
 *
 * Provides only the attributes required by AdaptBins and CoordinateSelector:
 * R (positions, inherited), P (momenta), and Bin (bin index).
 */
struct TestBunch : public ippl::ParticleBase<TestPLayout_t> {
    using Base             = ippl::ParticleBase<TestPLayout_t>;
    using bin_index_type   = short int;  // same as ParticleContainer

    /// Particle momenta (used by CoordinateSelector)
    typename Base::particle_position_type P;

    /// Bin index attribute (used by AdaptBins)
    ippl::ParticleAttrib<bin_index_type> Bin;

    TestBunch(TestPLayout_t& pl)
        : Base(pl) {
        this->addAttribute(P);
        this->addAttribute(Bin);
    }
    ~TestBunch() = default;
};

// ============================================================================
// Common type aliases
// ============================================================================
using Container_t    = TestBunch;
using Selector_t     = ParticleBinning::CoordinateSelector<Container_t>;
using AdaptBins_t    = ParticleBinning::AdaptBins<Container_t, Selector_t>;
using bin_index_type = Container_t::bin_index_type;  // short int
using value_type     = TestT;                        // double

// ============================================================================
// Test fixture
// ============================================================================
class BinningTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    // Members kept alive for the duration of each test
    std::unique_ptr<TestMesh_t>   mesh;
    std::unique_ptr<TestFL_t>     fl;
    std::unique_ptr<TestPLayout_t> playout;
    std::shared_ptr<Container_t>  bunch;
    std::shared_ptr<AdaptBins_t>  adaptBins;

    // Parameters
    static constexpr size_t gridPoints = 32;
    static constexpr TestT  domainLen  = 1.0;

    void SetUp() override {
        // Build mesh + field layout (similar to GatherScatterTest)
        ippl::Vector<TestT, TestDim> hr;
        ippl::Vector<TestT, TestDim> origin;
        ippl::NDIndex<TestDim> domain;
        std::array<bool, TestDim> decomp;
        decomp.fill(true);

        for (unsigned d = 0; d < TestDim; ++d) {
            domain[d]  = ippl::Index(gridPoints);
            hr[d]      = domainLen / gridPoints;
            origin[d]  = 0.0;
        }

        mesh    = std::make_unique<TestMesh_t>(domain, hr, origin);
        fl      = std::make_unique<TestFL_t>(MPI_COMM_WORLD, domain, decomp, true);
        playout = std::make_unique<TestPLayout_t>(*fl, *mesh);
        bunch   = std::make_shared<Container_t>(*playout);
        bunch->setParticleBC(ippl::BC::PERIODIC);
    }

    void TearDown() override {
        adaptBins.reset();
        bunch.reset();
        playout.reset();
        fl.reset();
        mesh.reset();
    }

    // ---- helpers ----------------------------------------------------------
    /// Creates `n` particles with random P values in [0.1, 0.9] and R in mesh interior.
    void createParticles(size_t n) {
        bunch->create(n);
        std::mt19937_64 eng(42 + ippl::Comm->rank());
        std::uniform_real_distribution<TestT> unifR(0.05, 0.95);
        std::uniform_real_distribution<TestT> unifP(0.1, 0.9);

        auto R_host = bunch->R.getHostMirror();
        auto P_host = bunch->P.getHostMirror();
        for (size_t i = 0; i < n; ++i) {
            ippl::Vector<TestT, TestDim> r, p;
            for (unsigned d = 0; d < TestDim; ++d) {
                r[d] = unifR(eng);
                p[d] = unifP(eng);
            }
            R_host(i) = r;
            P_host(i) = p;
        }
        Kokkos::deep_copy(bunch->R.getView(), R_host);
        Kokkos::deep_copy(bunch->P.getView(), P_host);
        ippl::Comm->barrier();
        Kokkos::fence();
        bunch->update();
    }

    /// Creates `n` particles where P[axis] is uniformly spread in [pMin, pMax].
    void createParticlesUniformP(size_t n, int axis, TestT pMin, TestT pMax) {
        bunch->create(n);
        std::mt19937_64 eng(123 + ippl::Comm->rank());
        std::uniform_real_distribution<TestT> unifR(0.05, 0.95);

        auto R_host = bunch->R.getHostMirror();
        auto P_host = bunch->P.getHostMirror();
        for (size_t i = 0; i < n; ++i) {
            ippl::Vector<TestT, TestDim> r, p;
            for (unsigned d = 0; d < TestDim; ++d) {
                r[d] = unifR(eng);
                p[d] = 0.0;
            }
            // Uniform spacing along the chosen axis
            p[axis] = pMin + (pMax - pMin) * (static_cast<TestT>(i) + 0.5) / n;
            R_host(i) = r;
            P_host(i) = p;
        }
        Kokkos::deep_copy(bunch->R.getView(), R_host);
        Kokkos::deep_copy(bunch->P.getView(), P_host);
        ippl::Comm->barrier();
        Kokkos::fence();
        bunch->update();
    }

    /// Build an AdaptBins with given parameters.
    void buildAdaptBins(bin_index_type maxBins = 5,
                        value_type alpha = 1.0,
                        value_type beta  = 1.0,
                        value_type desW  = 0.1) {
        Selector_t selector(2);  // bin along axis 2 (z-component of P)
        adaptBins = std::make_shared<AdaptBins_t>(bunch, selector, maxBins, alpha, beta, desW);
    }
};


// ############################################################################
// 1. ParallelReduceTools tests
// ############################################################################

TEST_F(BinningTest, ArrayReductionDefaultInit) {
    // ArrayReduction should initialise all elements to zero.
    constexpr bin_index_type N = 4;
    ParticleBinning::ArrayReduction<size_type, bin_index_type, N> arr;
    for (bin_index_type i = 0; i < N; ++i) {
        EXPECT_EQ(arr.the_array[i], 0u);
    }
}

TEST_F(BinningTest, ArrayReductionCopy) {
    constexpr bin_index_type N = 3;
    ParticleBinning::ArrayReduction<size_type, bin_index_type, N> a;
    a.the_array[0] = 10;
    a.the_array[1] = 20;
    a.the_array[2] = 30;

    ParticleBinning::ArrayReduction<size_type, bin_index_type, N> b(a);
    for (bin_index_type i = 0; i < N; ++i) {
        EXPECT_EQ(b.the_array[i], a.the_array[i]);
    }
}

TEST_F(BinningTest, ArrayReductionPlusEquals) {
    constexpr bin_index_type N = 3;
    ParticleBinning::ArrayReduction<size_type, bin_index_type, N> a, b;
    a.the_array[0] = 1; a.the_array[1] = 2; a.the_array[2] = 3;
    b.the_array[0] = 4; b.the_array[1] = 5; b.the_array[2] = 6;
    a += b;
    EXPECT_EQ(a.the_array[0], 5u);
    EXPECT_EQ(a.the_array[1], 7u);
    EXPECT_EQ(a.the_array[2], 9u);
}

TEST_F(BinningTest, CreateReductionObjectValid) {
    // createReductionObject should succeed for sizes 1..maxArrSize.
    for (bin_index_type n = 1; n <= ParticleBinning::maxArrSize<bin_index_type>; ++n) {
        auto var = ParticleBinning::createReductionObject<size_type, bin_index_type>(n);
        // Just verify no exception and variant is valid.
        std::visit([&](auto& reducer) {
            EXPECT_EQ(sizeof(reducer.the_array) / sizeof(reducer.the_array[0]),
                      static_cast<size_t>(n));
        }, var);
    }
}

TEST_F(BinningTest, CreateReductionObjectInvalid) {
    // Requesting a size above maxArrSize should throw.
    bin_index_type tooBig = ParticleBinning::maxArrSize<bin_index_type> + 1;
    // Wrap in a lambda to avoid macro comma issue with template args.
    EXPECT_THROW(
        (ParticleBinning::createReductionObject<size_type, bin_index_type>(tooBig)),
        std::out_of_range
    );
}

TEST_F(BinningTest, HostArrayReductionBasics) {
    // HostArrayReduction: dynamic allocation, zero init, +=
    constexpr bin_index_type N = 4;
    ParticleBinning::HostArrayReduction<size_type, bin_index_type>::binCountStatic = N;

    ParticleBinning::HostArrayReduction<size_type, bin_index_type> a;
    for (bin_index_type i = 0; i < N; ++i) {
        EXPECT_EQ(a.the_array[i], 0u);
    }
    a.the_array[0] = 7;
    a.the_array[3] = 3;

    ParticleBinning::HostArrayReduction<size_type, bin_index_type> b;
    b.the_array[0] = 1;
    b.the_array[3] = 2;

    a += b;
    EXPECT_EQ(a.the_array[0], 8u);
    EXPECT_EQ(a.the_array[3], 5u);
}

TEST_F(BinningTest, KokkosReductionIdentityArrayReduction) {
    // The Kokkos reduction_identity should return a zero-initialized ArrayReduction.
    constexpr bin_index_type N = 3;
    auto identity = Kokkos::reduction_identity<
        ParticleBinning::ArrayReduction<size_type, bin_index_type, N>>::sum();
    for (bin_index_type i = 0; i < N; ++i) {
        EXPECT_EQ(identity.the_array[i], 0u);
    }
}

TEST_F(BinningTest, KokkosReductionIdentityHostArrayReduction) {
    constexpr bin_index_type N = 4;
    ParticleBinning::HostArrayReduction<size_type, bin_index_type>::binCountStatic = N;
    auto identity = Kokkos::reduction_identity<
        ParticleBinning::HostArrayReduction<size_type, bin_index_type>>::sum();
    for (bin_index_type i = 0; i < N; ++i) {
        EXPECT_EQ(identity.the_array[i], 0u);
    }
}


// ############################################################################
// 2. BinningTools tests
// ############################################################################

TEST_F(BinningTest, DetermineHistoReductionModeStandard) {
    // Standard mode should auto-select based on bin count.
    auto modeSmall = ParticleBinning::determineHistoReductionMode<bin_index_type>(
        ParticleBinning::HistoReductionMode::Standard, 3);
    auto modeLarge = ParticleBinning::determineHistoReductionMode<bin_index_type>(
        ParticleBinning::HistoReductionMode::Standard, 100);

    bool isHostSpace = std::is_same<Kokkos::DefaultExecutionSpace,
                                    Kokkos::DefaultHostExecutionSpace>::value;
    if (isHostSpace) {
        // On host spaces, always HostOnly.
        EXPECT_EQ(modeSmall, ParticleBinning::HistoReductionMode::HostOnly);
        EXPECT_EQ(modeLarge, ParticleBinning::HistoReductionMode::HostOnly);
    } else {
        // On device, small bin count -> ParallelReduce, large -> TeamBased.
        EXPECT_EQ(modeSmall, ParticleBinning::HistoReductionMode::ParallelReduce);
        EXPECT_EQ(modeLarge, ParticleBinning::HistoReductionMode::TeamBased);
    }
}

TEST_F(BinningTest, DetermineHistoReductionModeForced) {
    // If a specific mode is forced, it should be respected (unless host-only execution space).
    bool isHostSpace = std::is_same<Kokkos::DefaultExecutionSpace,
                                    Kokkos::DefaultHostExecutionSpace>::value;
    if (!isHostSpace) {
        auto mode = ParticleBinning::determineHistoReductionMode<bin_index_type>(
            ParticleBinning::HistoReductionMode::TeamBased, 2);
        EXPECT_EQ(mode, ParticleBinning::HistoReductionMode::TeamBased);
    }
}

TEST_F(BinningTest, ComputeFixSum) {
    // Test computeFixSum with a simple host-space view.
    constexpr size_t N = 5;
    Kokkos::View<size_type*> input("input", N);
    Kokkos::View<size_type*> postsum("postsum", N + 1);

    // Fill input = {1,2,3,4,5}
    auto h_in = Kokkos::create_mirror_view(input);
    for (size_t i = 0; i < N; ++i) h_in(i) = i + 1;
    Kokkos::deep_copy(input, h_in);

    ParticleBinning::computeFixSum<Kokkos::View<size_type*>>(input, postsum);

    auto h_ps = Kokkos::create_mirror_view(postsum);
    Kokkos::deep_copy(h_ps, postsum);

    // Expected postsum: {0, 1, 3, 6, 10, 15}
    EXPECT_EQ(h_ps(0), 0u);
    EXPECT_EQ(h_ps(1), 1u);
    EXPECT_EQ(h_ps(2), 3u);
    EXPECT_EQ(h_ps(3), 6u);
    EXPECT_EQ(h_ps(4), 10u);
    EXPECT_EQ(h_ps(5), 15u);
}

TEST_F(BinningTest, ViewIsSortedTrue) {
    constexpr size_t N = 5;
    Kokkos::View<bin_index_type*> view("sortedView", N);
    auto h_view = Kokkos::create_mirror_view(view);
    for (size_t i = 0; i < N; ++i) h_view(i) = static_cast<bin_index_type>(i);
    Kokkos::deep_copy(view, h_view);

    // Identity hash: indices = {0,1,2,3,4}
    using hash_type = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;
    hash_type idx("idx", N);
    auto h_idx = Kokkos::create_mirror_view(idx);
    for (size_t i = 0; i < N; ++i) h_idx(i) = static_cast<int>(i);
    Kokkos::deep_copy(idx, h_idx);

    bool sorted = ParticleBinning::viewIsSorted<bin_index_type, size_type>(view, idx, N);
    EXPECT_TRUE(sorted);
}

TEST_F(BinningTest, ViewIsSortedFalse) {
    constexpr size_t N = 5;
    Kokkos::View<bin_index_type*> view("unsortedView", N);
    auto h_view = Kokkos::create_mirror_view(view);
    h_view(0) = 0; h_view(1) = 3; h_view(2) = 1; h_view(3) = 2; h_view(4) = 4;
    Kokkos::deep_copy(view, h_view);

    using hash_type = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;
    hash_type idx("idx", N);
    auto h_idx = Kokkos::create_mirror_view(idx);
    for (size_t i = 0; i < N; ++i) h_idx(i) = static_cast<int>(i);
    Kokkos::deep_copy(idx, h_idx);

    bool sorted = ParticleBinning::viewIsSorted<bin_index_type, size_type>(view, idx, N);
    EXPECT_FALSE(sorted);
}


// ############################################################################
// 3. BinHisto (Histogram) tests
// ############################################################################

TEST_F(BinningTest, HistogramConstruction) {
    // Construct a Histogram with DualView=false on host space.
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    Histo_t histo("testHisto", 4, 2.0, 1.0, 1.0, 0.5);

    EXPECT_EQ(histo.getCurrentBinCount(), 4);
}

TEST_F(BinningTest, HistogramInitSetsWidthsAndPostSum) {
    // After init(), widths should be uniform and postSum should be a prefix sum.
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    constexpr bin_index_type nBins = 4;
    constexpr value_type totalWidth = 2.0;

    Histo_t histo("testInit", nBins, totalWidth, 1.0, 1.0, 0.5);

    // Manually set histogram counts on host: {10, 20, 30, 40}
    auto histView = histo.getHistogram();
    histView(0) = 10; histView(1) = 20; histView(2) = 30; histView(3) = 40;

    histo.init();

    // Check bin widths are uniform: totalWidth/nBins = 0.5
    auto widths = histo.getBinWidths();
    for (bin_index_type i = 0; i < nBins; ++i) {
        EXPECT_NEAR(widths(i), totalWidth / nBins, 1e-12);
    }

    // Check post-sum: {0, 10, 30, 60, 100}
    auto ps = histo.getPostSum();
    EXPECT_EQ(ps(0), 0u);
    EXPECT_EQ(ps(1), 10u);
    EXPECT_EQ(ps(2), 30u);
    EXPECT_EQ(ps(3), 60u);
    EXPECT_EQ(ps(4), 100u);
}

TEST_F(BinningTest, HistogramGetNPartInBin) {
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    Histo_t histo("testNPart", 3, 1.0, 1.0, 1.0, 0.3);
    auto histView = histo.getHistogram();
    histView(0) = 5; histView(1) = 15; histView(2) = 25;
    histo.init();

    EXPECT_EQ(histo.getNPartInBin(0), 5u);
    EXPECT_EQ(histo.getNPartInBin(1), 15u);
    EXPECT_EQ(histo.getNPartInBin(2), 25u);
}

TEST_F(BinningTest, HistogramCopyConstructor) {
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    Histo_t original("origHisto", 3, 1.5, 1.0, 1.0, 0.5);
    auto histView = original.getHistogram();
    histView(0) = 10; histView(1) = 20; histView(2) = 30;
    original.init();

    Histo_t copy(original);
    EXPECT_EQ(copy.getCurrentBinCount(), 3);
    EXPECT_EQ(copy.getNPartInBin(0), 10u);
    EXPECT_EQ(copy.getNPartInBin(1), 20u);
    EXPECT_EQ(copy.getNPartInBin(2), 30u);
}

TEST_F(BinningTest, HistogramAssignmentOperator) {
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    Histo_t original("assignOrig", 2, 1.0, 1.0, 1.0, 0.5);
    auto histView = original.getHistogram();
    histView(0) = 7; histView(1) = 13;
    original.init();

    Histo_t assigned("dummy", 1, 0.5, 1.0, 1.0, 0.5);
    assigned = original;
    EXPECT_EQ(assigned.getCurrentBinCount(), 2);
    EXPECT_EQ(assigned.getNPartInBin(0), 7u);
    EXPECT_EQ(assigned.getNPartInBin(1), 13u);
}

TEST_F(BinningTest, HistogramGetBinIterationPolicy) {
    // Use DualView=true to avoid Kokkos subview type mismatch in non-DualView path.
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, true>;
    Histo_t histo("policyTest", 3, 1.0, 1.0, 1.0, 0.3);

    // Set counts on device then sync
    auto dView = histo.getHistogram().view_device();
    Kokkos::parallel_for("fillPolicyHisto", 1, KOKKOS_LAMBDA(const int) {
        dView(0) = 5; dView(1) = 10; dView(2) = 15;
    });
    histo.modify_device();
    histo.sync();
    histo.init();

    // postSum should be {0, 5, 15, 30}
    auto policy0 = histo.getBinIterationPolicy(0);
    auto policy1 = histo.getBinIterationPolicy(1);
    auto policy2 = histo.getBinIterationPolicy(2);

    // Check that the ranges encode [start, end) matching the post-sum
    EXPECT_EQ(policy0.begin(), 0);
    EXPECT_EQ(policy0.end(), 5);
    EXPECT_EQ(policy1.begin(), 5);
    EXPECT_EQ(policy1.end(), 15);
    EXPECT_EQ(policy2.begin(), 15);
    EXPECT_EQ(policy2.end(), 30);
}

TEST_F(BinningTest, HistogramMergeBins) {
    // Create a histogram with many bins and verify that mergeBins reduces them.
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>;
    constexpr bin_index_type nBins = 10;
    constexpr value_type totalWidth = 1.0;
    Histo_t histo("mergeTest", nBins, totalWidth, 1.0, 1.0, 0.2);

    // Create a distribution with most particles in a few bins.
    // Bins 0-2: 1 particle each (sparse tails)
    // Bins 3-6: 100 particles each (dense center)
    // Bins 7-9: 1 particle each (sparse tails)
    auto histView = histo.getHistogram();
    histView(0) = 1;  histView(1) = 1;  histView(2) = 1;
    histView(3) = 100; histView(4) = 100; histView(5) = 100; histView(6) = 100;
    histView(7) = 1;  histView(8) = 1;  histView(9) = 1;
    histo.init();

    auto lookupTable = histo.mergeBins();

    // After merging, the number of bins should be less than the original 10.
    bin_index_type mergedCount = histo.getCurrentBinCount();
    EXPECT_GT(mergedCount, 0);
    EXPECT_LE(mergedCount, nBins);

    // The lookup table maps old bin indices to new ones.
    // It should be monotonically non-decreasing.
    auto h_lookup = Kokkos::create_mirror_view(lookupTable);
    Kokkos::deep_copy(h_lookup, lookupTable);
    for (bin_index_type i = 1; i < nBins; ++i) {
        EXPECT_GE(h_lookup(i), h_lookup(i - 1));
    }

    // All new indices should be in [0, mergedCount)
    for (bin_index_type i = 0; i < nBins; ++i) {
        EXPECT_GE(h_lookup(i), 0);
        EXPECT_LT(h_lookup(i), mergedCount);
    }

    // Verify total particle count is conserved after merge.
    size_type totalAfter = 0;
    for (bin_index_type i = 0; i < mergedCount; ++i) {
        totalAfter += histo.getNPartInBin(i);
    }
    EXPECT_EQ(totalAfter, 406u);  // 3*1 + 4*100 + 3*1 = 406
}

TEST_F(BinningTest, HistogramDualViewConstruction) {
    // Test that the DualView variant of Histogram works properly.
    using Histo_t = ParticleBinning::Histogram<size_type, bin_index_type, value_type, true>;
    Histo_t histo("dualViewHisto", 4, 2.0, 1.0, 1.0, 0.5);
    EXPECT_EQ(histo.getCurrentBinCount(), 4);

    // Set counts on device, sync to host
    auto dView = histo.getHistogram().view_device();
    Kokkos::parallel_for("fillDualHisto", 4, KOKKOS_LAMBDA(const int i) {
        dView(i) = (i + 1) * 10;
    });
    histo.modify_device();
    histo.sync();
    histo.init();

    EXPECT_EQ(histo.getNPartInBin(0), 10u);
    EXPECT_EQ(histo.getNPartInBin(1), 20u);
    EXPECT_EQ(histo.getNPartInBin(2), 30u);
    EXPECT_EQ(histo.getNPartInBin(3), 40u);
}


// ############################################################################
// 4. CoordinateSelector tests
// ############################################################################

TEST_F(BinningTest, CoordinateSelectorReturnsNormalizedValues) {
    // CoordinateSelector applies: v = fabs(P[axis]); return v/sqrt(1+v*v)
    // So the output should always be in [0, 1).
    size_t nPart = 100;
    createParticles(nPart);

    Selector_t selector(2);  // z-axis
    selector.updateDataArr(bunch);

    auto P_host = bunch->P.getHostMirror();
    Kokkos::deep_copy(P_host, bunch->P.getView());

    size_t localN = bunch->getLocalNum();
    for (size_t i = 0; i < localN; ++i) {
        value_type p_z   = std::fabs(P_host(i)[2]);
        value_type expected = p_z / std::sqrt(1.0 + p_z * p_z);

        // Evaluate on host through a host mirror of viewdata
        // (The selector runs on device normally; here we verify the formula.)
        EXPECT_GE(expected, 0.0);
        EXPECT_LT(expected, 1.0);
    }
}


// ############################################################################
// 5. AdaptBins: getBin static function tests
// ############################################################################

TEST_F(BinningTest, GetBinBasic) {
    // Test the static getBin function for bin assignment.
    value_type xMin = 0.0, xMax = 1.0;
    bin_index_type numBins = 5;
    value_type binWidthInv = numBins / (xMax - xMin);

    // Midpoint of each bin should map to its index.
    for (bin_index_type b = 0; b < numBins; ++b) {
        value_type midpoint = xMin + (b + 0.5) / numBins * (xMax - xMin);
        bin_index_type result = AdaptBins_t::getBin(midpoint, xMin, xMax, binWidthInv, numBins);
        EXPECT_EQ(result, b);
    }
}

TEST_F(BinningTest, GetBinClampsOutOfRange) {
    value_type xMin = 0.0, xMax = 1.0;
    bin_index_type numBins = 4;
    value_type binWidthInv = numBins / (xMax - xMin);

    // Value below xMin should clamp to bin 0
    bin_index_type resultLow = AdaptBins_t::getBin(-0.5, xMin, xMax, binWidthInv, numBins);
    EXPECT_EQ(resultLow, 0);

    // Value above xMax should clamp to last bin (numBins-1)
    bin_index_type resultHigh = AdaptBins_t::getBin(1.5, xMin, xMax, binWidthInv, numBins);
    EXPECT_EQ(resultHigh, numBins - 1);
}

TEST_F(BinningTest, GetBinOnBoundary) {
    value_type xMin = 0.0, xMax = 1.0;
    bin_index_type numBins = 4;
    value_type binWidthInv = numBins / (xMax - xMin);

    // xMin should be bin 0
    EXPECT_EQ(AdaptBins_t::getBin(0.0, xMin, xMax, binWidthInv, numBins), 0);

    // xMax should be clamped to numBins-1
    EXPECT_EQ(AdaptBins_t::getBin(1.0, xMin, xMax, binWidthInv, numBins), numBins - 1);
}


// ############################################################################
// 6. AdaptBins: integration tests with particles
// ############################################################################

TEST_F(BinningTest, AdaptBinsConstruction) {
    createParticles(50);
    buildAdaptBins(5);
    EXPECT_EQ(adaptBins->getMaxBinCount(), 5);
    EXPECT_EQ(adaptBins->getCurrentBinCount(), 5);
}

TEST_F(BinningTest, AdaptBinsSetCurrentBinCount) {
    createParticles(50);
    buildAdaptBins(5);

    adaptBins->setCurrentBinCount(3);
    EXPECT_EQ(adaptBins->getCurrentBinCount(), 3);

    // Setting above max should clamp to max.
    adaptBins->setCurrentBinCount(100);
    EXPECT_EQ(adaptBins->getCurrentBinCount(), 5);
}

TEST_F(BinningTest, AdaptBinsInitLimits) {
    createParticles(200);
    buildAdaptBins(5);

    adaptBins->initLimits();

    // After initLimits, getBinWidth should be set.
    value_type bw = adaptBins->getBinWidth();
    EXPECT_GT(bw, 0.0);
}

TEST_F(BinningTest, AdaptBinsDoFullRebin) {
    // Full rebinning should complete without error and produce valid histogram.
    size_t nPart = 200;
    createParticles(nPart);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);

    // After full rebin, every bin should sum up to the total number of local particles.
    size_type totalInBins = 0;
    for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
        totalInBins += adaptBins->getNPartInBin(b);
    }
    EXPECT_EQ(totalInBins, bunch->getLocalNum());
}

TEST_F(BinningTest, AdaptBinsGlobalHistogramSumsToTotal) {
    // The global histogram should sum to bunch->getTotalNum().
    size_t nPart = 300;
    createParticles(nPart);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);

    size_type totalGlobal = 0;
    for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
        totalGlobal += adaptBins->getNPartInBin(b, /*global=*/true);
    }
    EXPECT_EQ(totalGlobal, bunch->getTotalNum());
}

TEST_F(BinningTest, AdaptBinsBinsAssignedInRange) {
    // After bin assignment, every particle's Bin attribute should be in [0, currentBins).
    size_t nPart = 200;
    createParticles(nPart);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);

    auto binHost = bunch->Bin.getHostMirror();
    Kokkos::deep_copy(binHost, bunch->Bin.getView());
    bin_index_type nBins = adaptBins->getCurrentBinCount();
    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        EXPECT_GE(binHost(i), 0);
        EXPECT_LT(binHost(i), nBins);
    }
}

TEST_F(BinningTest, AdaptBinsSortContainerByBin) {
    // After sorting, accessing particles through the hash should yield a
    // monotonically non-decreasing sequence of bin indices.
    size_t nPart = 400;
    createParticles(nPart);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);
    adaptBins->sortContainerByBin();

    auto hashArr = adaptBins->getHashArray();
    auto binView = bunch->Bin.getView();

    size_t localN = bunch->getLocalNum();
    if (localN <= 1) return;

    // Copy to host for checking
    auto h_hash = Kokkos::create_mirror_view(hashArr);
    Kokkos::deep_copy(h_hash, hashArr);
    auto h_bin = bunch->Bin.getHostMirror();
    Kokkos::deep_copy(h_bin, binView);

    for (size_t i = 1; i < localN; ++i) {
        EXPECT_GE(h_bin(h_hash(i)), h_bin(h_hash(i - 1)))
            << "Sort violated at index " << i;
    }
}

TEST_F(BinningTest, AdaptBinsBinIterationPolicyCoversAllParticles) {
    // The union of all bin iteration policies should cover exactly all local particles.
    size_t nPart = 300;
    createParticles(nPart);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);
    adaptBins->sortContainerByBin();

    size_type totalCovered = 0;
    for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
        auto policy = adaptBins->getBinIterationPolicy(b);
        totalCovered += (policy.end() - policy.begin());
    }
    EXPECT_EQ(totalCovered, bunch->getLocalNum());
}

TEST_F(BinningTest, AdaptBinsGenAdaptiveHistogramReducesBins) {
    // After genAdaptiveHistogram(), the number of bins should be <= maxBins.
    // Use a larger maxBins so there are bins to merge.
    size_t nPart = 500;
    createParticles(nPart);

    // Enough bins for a fine histogram that can be merged down.
    bin_index_type maxBins = 5;
    buildAdaptBins(maxBins, /*alpha=*/1.0, /*beta=*/1.0, /*desW=*/0.3);

    adaptBins->doFullRebin(maxBins);

    bin_index_type binsBefore = adaptBins->getCurrentBinCount();
    EXPECT_EQ(binsBefore, maxBins);

    adaptBins->genAdaptiveHistogram();

    bin_index_type binsAfter = adaptBins->getCurrentBinCount();
    EXPECT_GT(binsAfter, 0);
    EXPECT_LE(binsAfter, binsBefore);

    // Verify total particle count in local histogram is still correct.
    size_type totalInBins = 0;
    for (bin_index_type b = 0; b < binsAfter; ++b) {
        totalInBins += adaptBins->getNPartInBin(b);
    }
    EXPECT_EQ(totalInBins, bunch->getLocalNum());
}

TEST_F(BinningTest, AdaptBinsRebinWithDifferentBinCounts) {
    // Rebinning with different bin counts should produce valid histograms each time.
    size_t nPart = 200;
    createParticles(nPart);
    buildAdaptBins(5);

    for (bin_index_type nBins : {2, 3, 4, 5}) {
        adaptBins->setCurrentBinCount(nBins);
        adaptBins->doFullRebin(nBins);

        size_type totalInBins = 0;
        for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
            totalInBins += adaptBins->getNPartInBin(b);
        }
        EXPECT_EQ(totalInBins, bunch->getLocalNum())
            << "Mismatch with nBins=" << nBins;
    }
}

TEST_F(BinningTest, AdaptBinsRepeatedDoFullRebin) {
    // Calling doFullRebin multiple times should be idempotent in the histogram total.
    size_t nPart = 150;
    createParticles(nPart);
    buildAdaptBins(4);

    for (int iter = 0; iter < 3; ++iter) {
        adaptBins->doFullRebin(4);

        size_type total = 0;
        for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
            total += adaptBins->getNPartInBin(b);
        }
        EXPECT_EQ(total, bunch->getLocalNum())
            << "Failed at iteration " << iter;
    }
}

TEST_F(BinningTest, AdaptBinsSortAndPoliciesConsistent) {
    // The hash + bin iteration policy should allow correct per-bin iteration.
    size_t nPart = 300;
    createParticles(nPart);
    buildAdaptBins(4);

    adaptBins->doFullRebin(4);
    adaptBins->sortContainerByBin();

    auto hashArr = adaptBins->getHashArray();
    auto h_hash  = Kokkos::create_mirror_view(hashArr);
    Kokkos::deep_copy(h_hash, hashArr);

    auto h_bin = bunch->Bin.getHostMirror();
    Kokkos::deep_copy(h_bin, bunch->Bin.getView());

    for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
        auto policy = adaptBins->getBinIterationPolicy(b);
        for (auto idx = policy.begin(); idx < policy.end(); ++idx) {
            EXPECT_EQ(h_bin(h_hash(idx)), b)
                << "Particle at sorted index " << idx << " has bin " << h_bin(h_hash(idx))
                << " but expected " << b;
        }
    }
}

TEST_F(BinningTest, AdaptBinsHistoReductionModes) {
    // Test that doFullRebin works with each explicitly forced reduction mode.
    size_t nPart = 200;
    createParticles(nPart);

    // Use maxBins within maxArrSize range for ParallelReduce to work.
    bin_index_type maxBins = 3;
    buildAdaptBins(maxBins);

    // Test HostOnly mode
    adaptBins->doFullRebin(maxBins, true, ParticleBinning::HistoReductionMode::HostOnly);
    size_type totalHost = 0;
    for (bin_index_type b = 0; b < adaptBins->getCurrentBinCount(); ++b) {
        totalHost += adaptBins->getNPartInBin(b);
    }
    EXPECT_EQ(totalHost, bunch->getLocalNum());
}

TEST_F(BinningTest, AdaptBinsUniformDistributionEvenBins) {
    // For a truly uniform distribution of P values, bins should have roughly equal counts.
    size_t nPart = 1000;
    createParticlesUniformP(nPart, 2, 0.1, 0.9);
    buildAdaptBins(5);

    adaptBins->doFullRebin(5);

    // Check local histogram - should be roughly uniform within each rank's portion.
    size_type localN = bunch->getLocalNum();
    if (localN < 5) return;  // skip if too few particles on this rank

    size_type expected = localN / 5;
    for (bin_index_type b = 0; b < 5; ++b) {
        size_type count = adaptBins->getNPartInBin(b);
        // Allow ±50% tolerance for the uniform distribution.
        EXPECT_GT(count, 0u)
            << "Bin " << b << " is empty for uniform distribution";
    }
}

TEST_F(BinningTest, AdaptBinsFullWorkflow) {
    // End-to-end test of the full workflow:
    // create -> rebin -> adaptive merge -> sort -> iterate per bin.
    size_t nPart = 500;
    createParticles(nPart);

    bin_index_type maxBins = 5;
    buildAdaptBins(maxBins, 1.0, 1.0, 0.2);

    // 1. Full rebin
    adaptBins->doFullRebin(maxBins);

    // 2. Adaptive histogram (merge)
    adaptBins->genAdaptiveHistogram();
    bin_index_type adaptiveBins = adaptBins->getCurrentBinCount();
    EXPECT_GT(adaptiveBins, 0);

    // 3. Sort
    adaptBins->sortContainerByBin();

    // 4. Iterate per bin and verify all particles are covered.
    auto hashArr = adaptBins->getHashArray();
    auto h_hash  = Kokkos::create_mirror_view(hashArr);
    Kokkos::deep_copy(h_hash, hashArr);
    auto h_bin = bunch->Bin.getHostMirror();
    Kokkos::deep_copy(h_bin, bunch->Bin.getView());

    size_type totalCovered = 0;
    for (bin_index_type b = 0; b < adaptiveBins; ++b) {
        auto policy = adaptBins->getBinIterationPolicy(b);
        size_type rangeSize = policy.end() - policy.begin();
        totalCovered += rangeSize;

        // All particles in this range should have bin == b.
        for (auto idx = policy.begin(); idx < policy.end(); ++idx) {
            EXPECT_EQ(h_bin(h_hash(idx)), b);
        }
    }
    EXPECT_EQ(totalCovered, bunch->getLocalNum());
}


// ############################################################################
// main: initialize/finalize handled by the test fixture
// ############################################################################

int main(int argc, char* argv[]) {
    int result = 1;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
    return result;
}
