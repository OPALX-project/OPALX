#ifndef OPALX_RBend_HH
#define OPALX_RBend_HH

#include "AbsBeamline/BendBase.h"
#include "BeamlineGeometry/RBendGeometry.h"

/**
 * @class RBend
 * @brief Abstract rectangular bend with straight body and curved reference path.
 *
 * A rectangular bend keeps a straight hardware body while the reference orbit
 * is a circular arc tangent to the entrance and exit faces. In the standard
 * convention the exit edge angle is determined by the bend angle and the
 * entrance angle.
 */
class RBend : public BendBase {
public:
    RBend();
    explicit RBend(const std::string& name);
    RBend(const RBend&);
    ~RBend() override;

    void accept(BeamlineVisitor& visitor) const override;
    ElementType getType() const override;

    RBendGeometry& getGeometry() override             = 0;
    const RBendGeometry& getGeometry() const override = 0;

    CoordinateSystemTrafo getEdgeToBegin() const override;
    CoordinateSystemTrafo getEdgeToEnd() const override;

    double getExitAngle() const override;

protected:
    /**
     * @brief Return the OPAL-compatible rectangular-bend reference path length.
     *
     * For an `RBEND` the input length \f$L\f$ is the straight distance between
     * the rectangular pole faces.  OPAL's `Bend2D` computes the design radius
     * from the pole-face angles,
     * \f[
     *   R = \frac{L}{\sin E_1 + \sin E_2}, \qquad E_2 = \theta - E_1,
     * \f]
     * and normalizes the analytic field with
     * \f[
     *   L_\mathrm{ref} = |\theta|\,|R|.
     * \f]
     * The symmetric case \f$E_1=E_2=\theta/2\f$ reduces to the usual
     * chord-to-arc expression \f$L_\mathrm{ref}=L\,\theta/(2\sin(\theta/2))\f$.
     */
    double getReferencePathLength() const override;
};

#endif  // OPALX_RBend_HH
