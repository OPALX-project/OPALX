#ifndef OPAL_MATCH_CMD_H
#define OPAL_MATCH_CMD_H

#include <optional>
#include "AbstractObjects/Action.h"
#include "Algorithms/MatchedCovariance.h"

/** Named zero-current transverse MATCH of a successful static-magnetic RING COF.
 * ORBIT names the solve; MODEEMITTANCES_N supplies two normalized RMS mode
 * eigenemittances [m], no pi, ordered by ascending positive-imaginary phase.
 * OUTPUT is optional atomic version-1 JSON output. Reports the full mechanical
 * section covariance, conventions and residuals at information level 4.
 * Each execution consumes an independent by-value historical COF snapshot.
 * It does not follow later lattice edits/COF reruns, create a beam, or claim
 * a 6D/space-charge equilibrium. Only one bunch is in scope for later sampling.
 * Failed reruns and clones cannot expose a previous successful result.
 */
class MatchCmd : public Action {
public:
    struct Result {
        ClosedOrbitOptics source;
        MatchedCovariance::Result transverse;
        std::string json;
    };
    MatchCmd();
    MatchCmd* clone(const std::string&) override;
    void execute() override;
    /// Resolve a completed named MATCH; returns a read-only execution snapshot.
    static const Result& findResult(const std::string&);

private:
    MatchCmd(const std::string&, MatchCmd*);
    std::optional<Result> result_m;
};
#endif
