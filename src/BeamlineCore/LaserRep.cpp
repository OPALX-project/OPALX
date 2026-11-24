#include "BeamlineCore/LaserRep.h"

LaserRep::LaserRep() : Laser() {}

LaserRep::LaserRep(const LaserRep& right) : Laser(right) {}

LaserRep::LaserRep(const std::string& name) : Laser(name) {}

LaserRep::~LaserRep() {}
