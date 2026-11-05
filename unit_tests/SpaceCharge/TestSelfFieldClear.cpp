#include <gtest/gtest.h>

#include "SelfFieldClearTestOps.h"

#include <array>
#include <cstddef>
#include <memory>

namespace opalx::spacecharge::testing {
    namespace {

        class SelfFieldClearTest : public ::testing::Test {
        protected:
            using Point        = Vector_t<double, 3>;
            using MemorySpace  = typename decltype(SpaceChargeParticleContainer::E)::memory_space;
            using CapacityView = Kokkos::View<Point*, MemorySpace, Kokkos::MemoryUnmanaged>;

            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            void check(ClearOperation localClear, ClearOperation otherClear) {
                for (const std::size_t count : {0u, 3u, 257u}) {
                    SCOPED_TRACE(count);
                    ippl::NDIndex<3> domain;
                    for (unsigned d = 0; d < 3; ++d) {
                        domain[d] = ippl::Index(4);
                    }
                    Mesh_t<3> mesh(domain, Point(0.25), Point(-0.5));
                    FieldLayout_t<3> layout(
                            MPI_COMM_WORLD, domain, std::array<bool, 3>{false, false, false},
                            false);
                    SpaceChargeParticleContainer particles(mesh, layout);
                    particles.createParticles(count);

                    // Deliberately leave inactive capacity. getView() must expose live entries
                    // only.
                    particles.E.reserve(count + 11);
                    particles.B.reserve(count + 11);
                    const auto e = particles.E.getView();
                    const auto b = particles.B.getView();
                    const CapacityView eCapacity(e.data(), particles.E.size());
                    const CapacityView bCapacity(b.data(), particles.B.size());
                    ASSERT_EQ(e.extent(0), count);
                    ASSERT_EQ(b.extent(0), count);
                    ASSERT_GT(eCapacity.extent(0), count);
                    ASSERT_GT(bCapacity.extent(0), count);

                    Kokkos::deep_copy(particles.R.getView(), Point(0.125, -0.25, 0.375));
                    Kokkos::deep_copy(particles.P.getView(), Point(-1.0, 2.0, 3.0));
                    const auto rBefore = Kokkos::create_mirror(particles.R.getView());
                    const auto pBefore = Kokkos::create_mirror(particles.P.getView());
                    Kokkos::deep_copy(rBefore, particles.R.getView());
                    Kokkos::deep_copy(pBefore, particles.P.getView());

                    // Repeat and reverse call order to exercise both independently compiled
                    // callers.
                    for (const auto clear : {localClear, otherClear, otherClear, localClear}) {
                        Kokkos::deep_copy(eCapacity, Point(11.0, 12.0, 13.0));
                        Kokkos::deep_copy(bCapacity, Point(21.0, 22.0, 23.0));
                        clear(particles);
                        // This fence is observation-only; it occurs AFTER the operation under test.
                        Kokkos::fence("SelfFieldClearTest::observe");
                        const auto eHost =
                                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), eCapacity);
                        const auto bHost =
                                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), bCapacity);
                        const auto rHost = Kokkos::create_mirror_view_and_copy(
                                Kokkos::HostSpace(), particles.R.getView());
                        const auto pHost = Kokkos::create_mirror_view_and_copy(
                                Kokkos::HostSpace(), particles.P.getView());
                        for (std::size_t i = 0; i < eHost.extent(0); ++i) {
                            for (unsigned d = 0; d < 3; ++d) {
                                EXPECT_DOUBLE_EQ(eHost(i)[d], i < count ? 0.0 : 11.0 + d);
                            }
                        }
                        for (std::size_t i = 0; i < bHost.extent(0); ++i) {
                            for (unsigned d = 0; d < 3; ++d) {
                                EXPECT_DOUBLE_EQ(bHost(i)[d], i < count ? 0.0 : 21.0 + d);
                            }
                        }
                        for (std::size_t i = 0; i < count; ++i) {
                            for (unsigned d = 0; d < 3; ++d) {
                                EXPECT_DOUBLE_EQ(rHost(i)[d], rBefore(i)[d]);
                                EXPECT_DOUBLE_EQ(pHost(i)[d], pBefore(i)[d]);
                            }
                        }
                        EXPECT_EQ(particles.E.getView().data(), e.data());
                        EXPECT_EQ(particles.B.getView().data(), b.data());
                        EXPECT_EQ(particles.getLocalNum(), count);
                    }
                }
            }
        };

        TEST_F(SelfFieldClearTest, OriginalAssignment) {
            check(assignmentClear, assignmentClearOtherTranslationUnit);
        }

        TEST_F(SelfFieldClearTest, FencedAssignment) {
            check(fencedAssignmentClear, fencedAssignmentClearOtherTranslationUnit);
        }

        TEST_F(SelfFieldClearTest, DeepCopy) {
            check(clearSelfFields, deepCopyClearOtherTranslationUnit);
        }

    }  // namespace
}  // namespace opalx::spacecharge::testing
