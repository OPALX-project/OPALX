#include "MagneticField.h"

MagneticField::MagneticField()
    : Cyclotron()
    , nullGeom_m(NullGeometry())
    , nullField_m(NullField())
{
}

double MagneticField::getSlices() const {
    return -1.0;
}

double MagneticField::getStepsize() const {
    return -1.0;
}

const EMField &MagneticField::getField() const {
    return nullField_m;
}

EMField &MagneticField::getField() {
    return nullField_m;
}

ElementBase* MagneticField::clone() const {
    return nullptr;
}

const BGeometryBase &MagneticField::getGeometry() const {
    return nullGeom_m;
}

BGeometryBase &MagneticField::getGeometry() {
    return nullGeom_m;
}