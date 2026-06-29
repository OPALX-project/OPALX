#ifndef OPALX_RBend_HH
#define OPALX_RBend_HH

#include "AbsBeamline/BendBase.h"

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


    double getExitAngle() const override;
};

#endif  // OPALX_RBend_HH
