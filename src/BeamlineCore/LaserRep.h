#ifndef OPALX_LaserRep_HH
#define OPALX_LaserRep_HH

#include "AbsBeamline/Laser.h"
#include "BeamlineGeometry/StraightGeometry.h"

class LaserRep : public Laser {
public:
    explicit LaserRep(const std::string& name);
    LaserRep();
    LaserRep(const LaserRep&);
    ~LaserRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    StraightGeometry& getGeometry() override;
    const StraightGeometry& getGeometry() const override;

private:
    StraightGeometry geometry_m;
};

#endif  // OPALX_LaserRep_HH
