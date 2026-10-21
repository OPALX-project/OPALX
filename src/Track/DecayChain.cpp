#include "Track/DecayChain.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <iterator>
#include <limits>

std::vector<std::size_t> resolveDecayRoots(
        const std::vector<std::string>& names, const std::vector<std::string>& daughterNames,
        const std::vector<bool>& hasDecay) {
    const std::size_t n          = names.size();
    constexpr std::size_t noBeam = std::numeric_limits<std::size_t>::max();

    // parent[d] is the beam whose decay daughters go into beam d.
    std::vector<std::size_t> parent(n, noBeam);
    for (std::size_t i = 0; i < n; ++i) {
        const std::string& daughterName = daughterNames[i];
        if (daughterName.empty()) {
            continue;
        }
        if (!hasDecay[i]) {
            throw OpalException(
                    "resolveDecayRoots", "DAUGHTERBEAM=\"" + daughterName + "\" on beam \""
                                                 + names[i]
                                                 + "\" needs GLOBALPROCESSES = {\"DECAY\"}.");
        }
        const auto it = std::find(names.begin(), names.end(), daughterName);
        if (it == names.end()) {
            throw OpalException(
                    "resolveDecayRoots", "DAUGHTERBEAM=\"" + daughterName + "\" on beam \""
                                                 + names[i]
                                                 + "\" does not match any beam in the TRACK.");
        }
        const auto d = static_cast<std::size_t>(std::distance(names.begin(), it));
        if (d == i) {
            throw OpalException(
                    "resolveDecayRoots",
                    "Beam \"" + names[i] + "\" names itself as its DAUGHTERBEAM.");
        }
        if (parent[d] != noBeam) {
            throw OpalException(
                    "resolveDecayRoots", "Beam \"" + daughterName
                                                 + "\" is the DAUGHTERBEAM of both \""
                                                 + names[parent[d]] + "\" and \"" + names[i]
                                                 + "\"; a beam can receive the decay daughters "
                                                   "of only one beam.");
        }
        parent[d] = i;
    }

    std::vector<std::size_t> roots(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t root  = i;
        std::size_t steps = 0;
        while (parent[root] != noBeam) {
            root = parent[root];
            if (++steps > n) {
                throw OpalException(
                        "resolveDecayRoots",
                        "The DAUGHTERBEAM settings form a loop through beam \"" + names[i] + "\".");
            }
        }
        roots[i] = root;
    }
    return roots;
}
