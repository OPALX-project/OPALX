#include "BeamlineCore/LaserRep.h"

LaserRep::LaserRep() : Laser(), geometry_m() {}

LaserRep::LaserRep(const LaserRep& right) : Laser(right), geometry_m(right.geometry_m) {}

LaserRep::LaserRep(const std::string& name) : Laser(name), geometry_m() {}

LaserRep::~LaserRep() {}

ElementBase* LaserRep::clone() const { return new LaserRep(*this); }

Geometry& LaserRep::getGeometry() { return geometry_m; }

const Geometry& LaserRep::getGeometry() const { return geometry_m; }
