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
     * @brief Return the rectangular-bend design-path body length.
     *
     * For an `RBEND` the placed hardware body is straight, but the analytic
     * dipole field is normalized on the curved design reference path tangent to
     * the rotated entrance and exit faces. The corresponding path length is the
     * arc-equivalent length reported by `RBendGeometry`.
     */
    double getReferencePathLength() const override;
};

#endif  // OPALX_RBend_HH
