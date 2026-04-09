#include "Elements/OpalLaser.h"

#include "Attributes/Attributes.h"
#include "BeamlineCore/LaserRep.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <sstream>
#include <vector>

namespace {
    double getRequiredPositiveReal(const std::vector<Attribute>& attrs,
                                   int index,
                                   const std::string& where) {
        if (!attrs[index]) {
            throw OpalException(where, "\"" + attrs[index].getName() + "\" must be set.");
        }

        const double value = Attributes::getReal(attrs[index]);
        if (value <= 0.0) {
            throw OpalException(where, "\"" + attrs[index].getName() + "\" must be greater than 0.");
        }

        return value;
    }

    Vector_t<double, 3> getVector3(const Attribute& attr,
                                   const std::string& where,
                                   const std::string& semanticName) {
        if (!attr) {
            throw OpalException(where, "\"" + attr.getName() + "\" must be set.");
        }

        const std::vector<double> values = Attributes::getRealArray(attr);
        if (values.size() != 3) {
            std::ostringstream msg;
            msg << "\"" << attr.getName() << "\" must contain exactly 3 values for "
                << semanticName << "; got " << values.size() << ".";
            throw OpalException(where, msg.str());
        }

        return Vector_t<double, 3>(values[0], values[1], values[2]);
    }
}

OpalLaser::OpalLaser()
    : OpalElement(SIZE, "LASER",
                  "The \"LASER\" element defines a passive analytic laser pulse for OPALX.") {
    itsAttr[WAVELENGTH] = Attributes::makeReal("WAVELENGTH", "Laser wavelength in m.");
    itsAttr[PULSEENERGY] = Attributes::makeReal("PULSEENERGY", "Total laser pulse energy in J.");
    itsAttr[PULSELENGTH] = Attributes::makeReal("PULSELENGTH", "Laser pulse duration in s.");
    itsAttr[WAISTX] = Attributes::makeReal("WAISTX", "Horizontal laser waist in m.");
    itsAttr[WAISTY] = Attributes::makeReal("WAISTY", "Vertical laser waist in m.");
    itsAttr[DIR] = Attributes::makeRealArray("DIR", "Laser propagation direction vector.");
    itsAttr[STOKES] = Attributes::makeRealArray(
        "STOKES", "Normalized Stokes polarization vector {xi1, xi2, xi3}.");

    registerOwnership();

    setElement(new LaserRep("LASER"));
}

OpalLaser::OpalLaser(const std::string& name, OpalLaser* parent)
    : OpalElement(name, parent) {
    setElement(new LaserRep(name));
}

OpalLaser::~OpalLaser() {
}

OpalLaser* OpalLaser::clone(const std::string& name) {
    return new OpalLaser(name, this);
}

void OpalLaser::update() {
    constexpr const char* where = "OpalLaser::update()";

    OpalElement::update();

    auto* laser = dynamic_cast<LaserRep*>(getElement());
    if (!laser) {
        OpalElement::updateUnknown(getElement());
        return;
    }

    // Required attributes are unset on the prototype instance.
    // OpalData::update() calls update() on every object in the directory,
    // so skip propagation silently when nothing has been configured yet.
    //
    /// \todo Consider a more robust way to handle this, e.g. with default 
    /// values and then putting checks itself into the representation instance! 
    if (!itsAttr[WAVELENGTH]) {
        OpalElement::updateUnknown(laser);
        return;
    }

    const double wavelength = getRequiredPositiveReal(itsAttr, WAVELENGTH, where);
    const double pulseEnergy = getRequiredPositiveReal(itsAttr, PULSEENERGY, where);
    const double pulseLength = getRequiredPositiveReal(itsAttr, PULSELENGTH, where);
    const double waistX = getRequiredPositiveReal(itsAttr, WAISTX, where);
    const double waistY = getRequiredPositiveReal(itsAttr, WAISTY, where);

    Vector_t<double, 3> direction = getVector3(itsAttr[DIR], where, "the propagation direction");
    const double directionNorm = std::sqrt(dot(direction, direction));
    if (directionNorm <= 0.0) {
        throw OpalException(where, "\"DIR\" must be a non-zero vector.");
    }
    direction /= directionNorm;

    Vector_t<double, 3> stokes(0.0);
    if (itsAttr[STOKES]) {
        stokes = getVector3(itsAttr[STOKES], where, "the Stokes vector");
        const double stokesNorm2 = dot(stokes, stokes);
        if (stokesNorm2 > 1.0 + 1.0e-12) {
            throw OpalException(where, "\"STOKES\" must satisfy xi1^2 + xi2^2 + xi3^2 <= 1.");
        }
    }

    const double length = itsAttr[LENGTH] ? Attributes::getReal(itsAttr[LENGTH]) : 0.0;
    if (length < 0.0) {
        throw OpalException(where, "\"L\" must be greater than or equal to 0.");
    }

    laser->setElementLength(length);
    laser->setWavelength(wavelength);
    laser->setPulseEnergy(pulseEnergy);
    laser->setPulseLength(pulseLength);
    laser->setWaistX(waistX);
    laser->setWaistY(waistY);
    laser->setDirection(direction);
    laser->setStokes(stokes);

    OpalElement::updateUnknown(laser);
}
