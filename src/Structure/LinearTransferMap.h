// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved

#ifndef OPAL_LINEAR_TRANSFER_MAP_H
#define OPAL_LINEAR_TRANSFER_MAP_H

#include "Algorithms/Matrix.h"
#include "OPALTypes.h"

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Reference-orbit state and transported beam frame at a map boundary.
 *
 * Momentum is stored in OPALX units, \f$\vec{\beta}\gamma=\mathbf p/(mc)\f$.
 * The right-handed orthonormal triad \f$(\mathbf e_x,\mathbf e_y,\mathbf e_s)\f$ has
 * \f$\mathbf e_s=\mathbf p_0/|\mathbf p_0|\f$ except at a RING return section, where
 * the starting plane and axes are reused even for a mismatched return momentum.
 * Between reference steps the transverse
 * axes are advanced by the minimum rotation taking the old \f$\mathbf e_s\f$ into the new
 * one (a discrete approximation to Bishop transport). Its axis and angle are
 * \f[
 * \widehat{\mathbf a}=\frac{\mathbf e_s\times\mathbf e_s'}
 * {|\mathbf e_s\times\mathbf e_s'|},\qquad
 * \theta=\mathrm{atan2}(|\mathbf e_s\times\mathbf e_s'|,
 *                              \mathbf e_s\!\cdot\!\mathbf e_s'),
 * \f]
 * and the corresponding Rodrigues rotation is
 * \f[
 * R\mathbf v=\mathbf v\cos\theta+(\widehat{\mathbf a}\times\mathbf v)\sin\theta
 *       +\widehat{\mathbf a}(\widehat{\mathbf a}\cdot\mathbf v)(1-\cos\theta).
 * \f]
 * The implementation rotates \f$\mathbf e_x\f$, reprojects it perpendicular to the new
 * tangent, normalizes it and reconstructs \f$\mathbf e_y=\mathbf e_s\times\mathbf e_x\f$.
 * Parallel tangents need no rotation; at an antiparallel reversal the old x axis is
 * retained to resolve the non-unique rotation axis. This introduces no arbitrary
 * tangent roll but is not an exact finite-step solution for a general curved orbit.
 *
 * References: R. L. Bishop, "There is More than One Way to Frame a Curve",
 * American Mathematical Monthly 82(3), 246-251 (1975),
 * [doi:10.2307/2319846](https://doi.org/10.2307/2319846); for discrete rotation
 * transport, W. Wang et al.,
 * "Computation of Rotation Minimizing Frames", ACM TOG 27(1), Article 2 (2008),
 * section 2.1, https://doi.org/10.1145/1330511.1330513.
 * Rodrigues' formula is also given in Lynch and Park, Modern Robotics, section 3.2.3:
 * https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-2-3-exponential-coordinates-of-rotation-part-2-of-2/
 */
struct LinearTransferMapReference {
    /// Reference position in lab coordinates [m]; see positionCorrection for the low part.
    Vector_t<double, 3> position = Vector_t<double, 3>(0.0);
    /// Lab mechanical momentum \f$\mathbf p/(mc)=\vec\beta\gamma\f$, dimensionless.
    Vector_t<double, 3> momentum = Vector_t<double, 3>(0.0);
    /// Horizontal unit vector of the section, expressed in lab coordinates.
    Vector_t<double, 3> xAxis{1.0, 0.0, 0.0};
    /// Vertical unit vector of the section, expressed in lab coordinates.
    Vector_t<double, 3> yAxis{0.0, 1.0, 0.0};
    /// Section normal in lab coordinates; normally the reference momentum direction.
    Vector_t<double, 3> sAxis{0.0, 0.0, 1.0};
    /// Reference arrival time [s], not the integration time step.
    double time{0.0};
    /// Accumulated signed reference-path coordinate [m], not a nominal element length.
    double pathLength{0.0};
    /// Kahan position residual [m]; represented value = position - positionCorrection.
    Vector_t<double, 3> positionCorrection = Vector_t<double, 3>(0.0);
    /// Kahan time residual [s]; represented value = time - timeCorrection.
    double timeCorrection{0.0};
    /// Kahan path residual [m]; represented value = pathLength - pathLengthCorrection.
    double pathLengthCorrection{0.0};
};

/**
 * @brief First-order external-field transfer map of one reference-path segment.
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
 * Here the momentum components are projections onto the boundary's stored frame,
 * \f$p_0\f$ is the reference momentum magnitude and \f$\beta_0\f$ its speed divided by c.
 * All positions are in metres, slopes and \f$\delta\f$ are dimensionless.
 * The matrix is a differential map (no affine offset is stored). An entry has the units
 * of its output coordinate divided by those of its input coordinate; R16 and R26 mean
 * matrix(0,5) [m] and matrix(1,5) [dimensionless], respectively.
 * Each stored matrix covers one unique path segment over which the nominal body-owner set is
 * constant. A copy is attached to every nominal owner in that segment. Consequently an element
 * can own several consecutive segment maps, while an overlap segment is calculated only once for
 * ordered beamline composition. @c activeElements records sampled field-support contributors,
 * and may differ from the owners (see its sampling definition below). A drift map includes any
 * neighbouring fringe fields; no fringe-support length is added to the nominal body. Shadow rays
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
 * including unowned intervals (which can contain field tails). A shared overlap segment must occur only once in
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
 * \f$O(\epsilon^2)\f$. Optional Richardson refinement is configured through
 * LinearTransferMapBuilder::Settings. Each level halves all six amplitudes;
 * extrapolation acts on the segment matrix before ordered composition. It does
 * not increase integration order or guarantee accuracy or symplecticity.
 *
 * Field-support crossings are resolved independently for each ray by
 * ExternalFieldRayTracker; segment ownership does not freeze its field selection.
 * External fields from all active elements are summed during ray integration; collective fields
 * and RF structures are excluded. @c determinantResidual stores \f$|\det(M)-1|\f$, a check
 * of volume preservation in the reported coordinates. Unit determinant is necessary but
 * not sufficient for canonical symplecticity. Since the chosen slopes and mechanical
 * momenta are not globally canonical, neither residual is a universal invariant of
 * physical tracking in these coordinates. @c symplecticResidual stores
 * \f$\max_{i,j}|(M^TJM-J)_{ij}|\f$ for the block-diagonal canonical matrix
 * \f$J=\mathrm{diag}(J_2,J_2,J_2)\f$ with
 * \f$J_2=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\f$.
 */
struct LinearTransferMap {
    /// Segment derivative in the coordinate order above; zero-initialized, not an identity map.
    matrix6x6_t matrix{0.0};
    /// Starting positive perturbation amplitudes (x,x',y,y',zeta,delta), not time steps.
    /// Metres in entries 0,2,4; dimensionless in entries 1,3,5.
    std::array<double, 6> finiteDifferenceSteps{};
    /// Smallest amplitudes used, equal to finiteDifferenceSteps / 2^richardsonLevels.
    std::array<double, 6> finestFiniteDifferenceSteps{};
    /// Number of additional amplitude halvings/extrapolation levels; zero is centered differences.
    unsigned richardsonLevels{0};
    /// Per-column maximum change between the last two tableau diagonals; NOT an error bound.
    /// Absent for level zero. Matrix entries have mixed physical units.
    std::optional<std::array<double, 6>> richardsonCorrection;
    /// Integrator used for reference and shadow rays: BORIS (default), RK4 or DOP853.
    std::string integrationMethod{"BORIS"};
    /// Input section and normalization reference; fields are not required to vanish here.
    LinearTransferMapReference entrance;
    /// Output section and normalization reference; each ray is tracked to this plane.
    LinearTransferMapReference exit;
    /// Zero-based identity within one builder result / OrbitThreader execution.
    /// Shared-owner copies retain this ID for deduplication; max(size_t) means unassigned.
    std::size_t segment{std::numeric_limits<std::size_t>::max()};
    /// Zero-based attachment ordinal within one element's map vector, NOT a turn number.
    /// Initially the segment ordinal in the builder; OrbitThreader assigns it per owner.
    std::size_t pass{0};
    /// Sorted names from support queries at chord midpoints of reference sample brackets
    /// intersecting this interval. Descriptive only: not the owner set or a field-selection mask.
    /// Brackets may straddle a body boundary; repeated occurrence names are not unique IDs.
    std::vector<std::string> activeElements;
    /// Maximum over refinement levels of
    /// \f$\|D_{\rm in}\|_\infty\,\|D_{\rm in}^{-1}\|_\infty\f$.
    /// Uses the unscaled, mixed-unit coordinate matrix; not a condition number of the optics.
    double inputConditionNumber{0.0};
    /// Absolute value |det(matrix)-1|, not the signed determinant difference.
    double determinantResidual{0.0};
    /// Maximum absolute entry of M^T J M - J, not its determinant or a row-sum norm.
    double symplecticResidual{0.0};
    /// True for multiple nominal owners OR a sampled contributor set of size > 1.
    /// Support participation does not require nonzero fields; not a 3D overlap certificate.
    bool includesOverlappingFields{false};
};

#endif  // OPAL_LINEAR_TRANSFER_MAP_H
