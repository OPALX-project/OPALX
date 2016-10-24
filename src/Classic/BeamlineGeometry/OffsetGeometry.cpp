// ------------------------------------------------------------------------
// $RCSfile: OffsetGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OffsetGeometry
//    Represents a geometry which is offset with respect to some global
//    (super) geometry.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/OffsetGeometry.h"


// Class OffsetGeometry.
// ------------------------------------------------------------------------

OffsetGeometry::OffsetGeometry(const BGeometryBase &g, const BGeometryBase &l,
                               const Euclid3D &t):
    global(g), local(l), g2l(t)
{}


OffsetGeometry::OffsetGeometry(const BGeometryBase &g, const Euclid3D &t):
    global(g), local(g), g2l(t)
{}


OffsetGeometry::OffsetGeometry(const OffsetGeometry &og):
    global(og.global), local(og.local), g2l(og.g2l)
{}


OffsetGeometry::~OffsetGeometry()
{}


double OffsetGeometry::getArcLength() const {
    return global.getArcLength();
}


double OffsetGeometry::getElementLength() const {
    return global.getElementLength();
}


Euclid3D OffsetGeometry::getGtoL() const {
    return g2l;
}


void OffsetGeometry::setGtoL(const Euclid3D &t) {
    g2l = t;
}


double OffsetGeometry::getOrigin() const {
    return global.getOrigin();
}


double OffsetGeometry::getEntrance() const {
    return global.getEntrance();
}


double OffsetGeometry::getExit() const {
    return global.getExit();
}


Euclid3D OffsetGeometry::getTransform(double fromS, double toS) const {
    return global.getTransform(fromS, toS);
}


Euclid3D OffsetGeometry::getTotalTransform() const {
    return global.getTotalTransform();
}


Euclid3D OffsetGeometry::getTransform(double s) const {
    return global.getTransform(s);
}


Euclid3D OffsetGeometry::getEntranceFrame() const {
    return global.getEntranceFrame();
}


Euclid3D OffsetGeometry::getExitFrame() const {
    return global.getExitFrame();
}


Euclid3D OffsetGeometry::getEntrancePatch() const {
    return getGlobalToLocalTransform(global.getEntrance(), local.getEntrance());
}


Euclid3D OffsetGeometry::getExitPatch() const {
    return Inverse(getGlobalToLocalTransform(global.getExit(), local.getExit()));
}


Euclid3D OffsetGeometry::getGlobalToLocalTransform
(double gs, double ls) const {
    Euclid3D tl = local.getTransform(ls);
    Euclid3D tgi = Inverse(local.getTransform(gs));

    return tl * g2l * tgi;
}

