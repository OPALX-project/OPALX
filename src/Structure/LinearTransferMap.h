// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved

#ifndef OPAL_LINEAR_TRANSFER_MAP_H
#define OPAL_LINEAR_TRANSFER_MAP_H

#include "Algorithms/Matrix.h"
#include "OPALTypes.h"

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

/**
 * @brief Reference-orbit state and continuously transported beam frame at a map boundary.
 *
 * Momentum is stored in OPALX units, \f$\vec{\beta}\gamma=\mathbf p/(mc)\f$.
 * The right-handed orthonormal triad \f$(\mathbf e_x,\mathbf e_y,\mathbf e_s)\f$ has
 * \f$\mathbf e_s=\mathbf p_0/|\mathbf p_0|\f$.  Between reference steps the transverse
 * axes are advanced by the minimum rotation taking the old \f$\mathbf e_s\f$ into the new
 * one (Bishop transport).  Its axis and angle are
 * \f[
 * \widehat{\mathbf a}=\frac{\mathbf e_s\times\mathbf e_s'}
 * {|\mathbf e_s\times\mathbf e_s'|},\qquad
 * \theta=\mathrm{atan2}(|\mathbf e_s\times\mathbf e_s'|,
 *                              \mathbf e_s\!\cdot\!\mathbf e_s'),
 * \f]
 * and the same Rodrigues rotation is applied to \f$\mathbf e_x\f$ and
 * \f$\mathbf e_y\f$.  Thus no arbitrary rotation about the reference tangent is introduced.
 */
struct LinearTransferMapReference {
    Vector_t<double, 3> position = Vector_t<double, 3>(0.0);
    Vector_t<double, 3> momentum = Vector_t<double, 3>(0.0);
    Vector_t<double, 3> xAxis{1.0, 0.0, 0.0};
    Vector_t<double, 3> yAxis{0.0, 1.0, 0.0};
    Vector_t<double, 3> sAxis{0.0, 0.0, 1.0};
    double time{0.0};
    double pathLength{0.0};
};

/**
 * @brief First-order external-field transfer map attached to one runtime element occurrence.
 *
 * The coordinates are
 * \f[
 * X=(x,x',y,y',\zeta,\delta),
 * \qquad x'=p_x/p_s,\quad y'=p_y/p_s,
 * \f]
 * \f[
 * \zeta=-\beta_0c(t-t_0),
 * \qquad \delta=|\mathbf p|/p_0-1.
 * \f]
 * All positions are in metres, slopes and \f$\delta\f$ are dimensionless.
 * Each stored matrix covers one unique path segment over which the active-element set is
 * constant. A copy is attached to every element active in that segment. Consequently an element
 * can own several consecutive segment maps, while an overlap segment is calculated only once for
 * ordered beamline composition. `activeElements` identifies that reference set. Shadow rays
 * instead select the active set \f$A(\mathbf r)\f$ at their own positions and see
 * \f[
 *   \mathbf E(\mathbf r,t)=\sum_{i\in A(\mathbf r)}\mathbf E_i(\mathbf r,t),
 *   \qquad
 *   \mathbf B(\mathbf r,t)=\sum_{i\in A(\mathbf r)}\mathbf B_i(\mathbf r,t).
 * \f]
 * Since the segments are disjoint and path ordered, the complete map is
 * \f[
 *   M_{\mathrm{total}}=M_N M_{N-1}\cdots M_1,
 * \f]
 * including field-free segment maps. A shared overlap segment must therefore occur only once in
 * this product, even though a copy is attached to every participating element.
 *
 * For each coordinate \f$j\f$, two private rays are launched with \f$\pm\epsilon_j\f$.  With
 * \f[
 * D_{\rm in}=\left[\frac{X_1^+-X_1^-}{2}\;\cdots\;
 *                         \frac{X_6^+-X_6^-}{2}\right]
 * \f]
 * and the analogous \f$D_{\rm out}\f$, the stored matrix is obtained from
 * \f[
 * M D_{\rm in}=D_{\rm out},
 * \f]
 * using a pivoted solve.  Centered differences have truncation error
 * \f$O(\epsilon^2)\f$.
 *
 * Field-support crossings are resolved independently for each ray by
 * ExternalFieldRayTracker; segment ownership does not freeze its field selection.
 * External fields from all active elements are summed during ray integration; collective fields
 * and RF structures are excluded. @c determinantResidual stores \f$|\det(M)-1|\f$, which checks
 * phase-space volume preservation.
 * This is necessary, but for a \f$6\times6\f$ map not sufficient, for symplecticity. Since the
 * chosen slopes and mechanical momenta are not globally canonical, @c symplecticResidual is a
 * diagnostic rather than a universal invariant. It stores
 * \f$\max_{i,j}|(M^TJM-J)_{ij}|\f$ for the block-diagonal canonical matrix
 * \f$J=\mathrm{diag}(J_2,J_2,J_2)\f$.
 */
struct LinearTransferMap {
    matrix6x6_t matrix{0.0};
    std::array<double, 6> finiteDifferenceSteps{};
    LinearTransferMapReference entrance;
    LinearTransferMapReference exit;
    /// Unique segment ordinal within this OrbitThreader execution.
    std::size_t segment{std::numeric_limits<std::size_t>::max()};
    /// Per-element attachment ordinal (useful for a ring-seam split or repeated occurrence).
    std::size_t pass{0};
    /// Names of all elements active throughout this segment.
    std::vector<std::string> activeElements;
    double inputConditionNumber{0.0};
    double determinantResidual{0.0};
    double symplecticResidual{0.0};
    bool includesOverlappingFields{false};
};

#endif  // OPAL_LINEAR_TRANSFER_MAP_H
