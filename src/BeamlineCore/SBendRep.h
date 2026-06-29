#ifndef OPALX_SBendRep_HH
#define OPALX_SBendRep_HH

#include "AbsBeamline/SBend.h"

/**
 * @class SBendRep
 * @brief Concrete OPALX representation of an analytic sector bend.
 *
 * The representation owns the planar-arc body geometry and the local analytic
 * multipole field expansion used by the first OPALX-native `SBEND` port.
 */
class SBendRep : public SBend {
public:
    SBendRep();
    explicit SBendRep(const std::string& name);
    SBendRep(const SBendRep&);
    ~SBendRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    BGeometryBase& getGeometry() override;
    const BGeometryBase& getGeometry() const override;

private:
    Geometry geometry_m{Geometry::makeArc(0.0, 0.0)};
};

#endif  // OPALX_SBendRep_HH
