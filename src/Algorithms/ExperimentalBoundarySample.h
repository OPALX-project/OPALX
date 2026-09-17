#ifndef OPALX_EXPERIMENTAL_BOUNDARY_SAMPLE_H
#define OPALX_EXPERIMENTAL_BOUNDARY_SAMPLE_H

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <limits>

/** Internal, opt-in experiment for sampling shared boundary-refinement triggers.
 * This does not sample deposition, field application, or Boris advancement.
 * Unselected hard-edge crossings have a timestep-dependent integration error;
 * the caller must keep the reference trigger and validate timestep/sample size.
 * No input-language option or supported production accuracy guarantee is added.
 */
namespace experimental_boundary {

    /** Parse an unsigned power-of-two stride. Null/0/1 selects the full population.
     * A zero return denotes malformed input, including overflow or whitespace.
     * Kept separate from getenv/MPI so validation can be tested without OPALX runs.
     */
    inline std::uint64_t parseStride(const char* text) {
        if (!text) return 1;
        if (!*text) return 0;
        std::uint64_t value = 0;
        for (; *text; ++text) {
            if (*text < '0' || *text > '9') return 0;
            const unsigned digit = *text - '0';
            if (value > (std::numeric_limits<std::uint64_t>::max() - digit) / 10) return 0;
            value = 10 * value + digit;
        }
        if (value == 0) return 1;
        return (value & (value - 1)) == 0 ? value : 0;
    }

    /** SplitMix64 finalizer applied to the stable particle ID, with fixed offset.
     * Unsigned overflow is intentional. The device performs no RNG or data copy.
     */
    KOKKOS_INLINE_FUNCTION std::uint64_t hash(std::uint64_t id) {
        auto value = id + UINT64_C(0x9e3779b97f4a7c15);
        value      = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        value      = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31);
    }

    /** Nested ID subsets with expected population N/stride; stride must be valid.
     * Ownership and array reordering do not change a particle's membership.
     */
    KOKKOS_INLINE_FUNCTION bool selected(std::uint64_t id, std::uint64_t stride) {
        return stride == 1 || (hash(id) & (stride - 1)) == 0;
    }

}  // namespace experimental_boundary
#endif
