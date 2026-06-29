#ifndef OPALX_LaserRep_HH
#define OPALX_LaserRep_HH

#include "AbsBeamline/Laser.h"
#include "BeamlineGeometry/Geometry.h"

class LaserRep : public Laser {
public:
    explicit LaserRep(const std::string& name);
    LaserRep();
    LaserRep(const LaserRep&);
    ~LaserRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    BGeometryBase& getGeometry() override;
    const BGeometryBase& getGeometry() const override;

private:
    Geometry geometry_m;
};

#endif  // OPALX_LaserRep_HH
