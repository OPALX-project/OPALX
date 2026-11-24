#ifndef OPALX_LaserRep_HH
#define OPALX_LaserRep_HH

#include "AbsBeamline/Laser.h"

class LaserRep : public Laser {
public:
    explicit LaserRep(const std::string& name);
    LaserRep();
    LaserRep(const LaserRep&);
    ~LaserRep() override;
};

#endif  // OPALX_LaserRep_HH
