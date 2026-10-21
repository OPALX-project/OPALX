#ifndef OPALX_TRACK_DECAYCHAIN_H
#define OPALX_TRACK_DECAYCHAIN_H

#include <cstddef>
#include <string>
#include <vector>

/**
 * @brief Find, for every beam of a TRACK, the beam whose reference particle it follows.
 *
 * A beam named as DAUGHTERBEAM by another beam only receives that beam's decay daughters. It
 * has no reference particle of its own: it follows the first beam of its decay chain. For
 * pi -> mu -> e both the muon and the electron beam follow the pion beam. Every other beam
 * follows itself, so beams that are not linked by DAUGHTERBEAM stay independent.
 *
 * @param names Beam names in TRACK BEAMS order.
 * @param daughterNames DAUGHTERBEAM of each beam, empty when not set.
 * @param hasDecay Whether each beam has GLOBALPROCESSES = DECAY.
 * @return For each beam, the index of the beam it follows (its own index if it is not a
 *   daughter).
 * @throws OpalException if a DAUGHTERBEAM is set without DECAY, is not in the TRACK, names the
 *   beam itself or a beam that already receives daughters from another beam, or if the
 *   DAUGHTERBEAM settings form a loop.
 */
std::vector<std::size_t> resolveDecayRoots(
        const std::vector<std::string>& names, const std::vector<std::string>& daughterNames,
        const std::vector<bool>& hasDecay);

#endif
