#ifndef OPALX_SBend_HH
#define OPALX_SBend_HH

#include "AbsBeamline/BendBase.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"

/**
 * @class SBend
 * @brief Abstract sector bend with planar-arc body geometry.
 *
 * A sector bend is represented by a planar circular arc in the local x-z
 * plane. The entrance and exit faces may have independent edge angles, so the
 * canonical body geometry and the edge-face optics data are kept separately.
 */
class SBend : public BendBase {
public:
    SBend();
    explicit SBend(const std::string& name);
    SBend(const SBend&);
    ~SBend() override;

    void accept(BeamlineVisitor& visitor) const override;
    ElementType getType() const override;

    PlanarArcGeometry& getGeometry() override             = 0;
    const PlanarArcGeometry& getGeometry() const override = 0;

    void setExitAngle(double exitAngle) override;
    double getExitAngle() const override;
};

#endif  // OPALX_SBend_HH
