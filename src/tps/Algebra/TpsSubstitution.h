// Copyright: see LICENSE (original OPAL/CLASSIC sources).
#ifndef TPS_TPS_SUBSTITUTION_H
#define TPS_TPS_SUBSTITUTION_H
namespace tps {
    /** @brief Internal depth-first monomial traversal record for substitution.
     * @details Multiplying the current parent by coordinate variable produces the
     * monomial at coefficient index. skip is the first record after this node's
     * subtree, allowing traversal to omit terms above the requested order.
     */
    struct TpsSubstitution {
        int index    = 0;  ///< Coefficient index in the monomial table.
        int order    = 0;  ///< Total degree of this monomial.
        int variable = 0;  ///< Coordinate multiplied into the parent monomial.
        int skip     = 0;  ///< Index immediately after the subtree.
    };
}  // namespace tps
#endif
