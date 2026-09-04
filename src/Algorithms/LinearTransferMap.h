// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved

#ifndef OPAL_LINEAR_TRANSFER_MAP_H
#define OPAL_LINEAR_TRANSFER_MAP_H

#include "Algorithms/Matrix.h"
#include "OPALTypes.h"

#include <array>
#include <cstddef>

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
    Vector_t<double, 3> position{0.0};
    Vector_t<double, 3> momentum{0.0};
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
 * All positions are in metres, slopes and \f$\delta\f$ are dimensionless.  The entrance and
 * exit are the first and last reference-orbit crossings of the element's field-support region;
 * hence a tabulated or analytic fringe field is included in the element map.
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
 * This map excludes collective fields.  In the first implementation RF structures and
 * overlapping elements are rejected.  Since the chosen slopes and mechanical momenta are not
 * globally canonical, `symplecticResidual` is a diagnostic rather than a universal invariant.
 * It stores the maximum component of \f$M^TJM-J\f$ for the block-diagonal canonical matrix
 * \f$J=\mathrm{diag}(J_2,J_2,J_2)\f$.
 */
struct LinearTransferMap {
    matrix6x6_t matrix{0.0};
    std::array<double, 6> finiteDifferenceSteps{};
    LinearTransferMapReference entrance;
    LinearTransferMapReference exit;
    std::size_t pass{0};
    double inputConditionNumber{0.0};
    double symplecticResidual{0.0};
    bool includesOverlappingFields{false};
};

#endif  // OPAL_LINEAR_TRANSFER_MAP_H
