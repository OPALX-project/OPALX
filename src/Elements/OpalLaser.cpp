#include "Elements/OpalLaser.h"

#include <memory>

#include "Attributes/Attributes.h"
#include "BeamlineCore/LaserRep.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <sstream>
#include <vector>

namespace {
    double getRequiredPositiveReal(
            const ElementValues& values, const std::string& name, const std::string& where) {
        if (!values.hasValue(name)) {
            throw OpalException(where, "\"" + name + "\" must be set.");
        }

        const double value = values.getReal(name);
        if (value <= 0.0) {
            throw OpalException(where, "\"" + name + "\" must be greater than 0.");
        }

        return value;
    }

    Vector_t<double, 3> getVector3(
            const ElementValues& values, const std::string& name, const std::string& where,
            const std::string& semanticName) {
        if (!values.hasValue(name)) {
            throw OpalException(where, "\"" + name + "\" must be set.");
        }

        const std::vector<double> entries = values.getRealArray(name);
        if (entries.size() != 3) {
            std::ostringstream msg;
            msg << "\"" << name << "\" must contain exactly 3 values for " << semanticName
                << "; got " << entries.size() << ".";
            throw OpalException(where, msg.str());
        }

        return Vector_t<double, 3>(entries[0], entries[1], entries[2]);
    }
}  // namespace

OpalLaser::OpalLaser()
    : OpalElement(
              SIZE, "LASER",
              "The \"LASER\" element defines a passive analytic laser pulse for OPALX.") {
    itsAttr[WAVELENGTH]  = Attributes::makeReal("WAVELENGTH", "Laser wavelength in m.");
    itsAttr[PULSEENERGY] = Attributes::makeReal("PULSEENERGY", "Total laser pulse energy in J.");
    itsAttr[PULSELENGTH] = Attributes::makeReal("PULSELENGTH", "Laser pulse duration in s.");
    itsAttr[WAISTX]      = Attributes::makeReal("WAISTX", "Horizontal laser waist in m.");
    itsAttr[WAISTY]      = Attributes::makeReal("WAISTY", "Vertical laser waist in m.");
    itsAttr[DIR]         = Attributes::makeRealArray("DIR", "Laser propagation direction vector.");
    itsAttr[STOKES]      = Attributes::makeRealArray(
            "STOKES", "Normalized Stokes polarization vector {xi1, xi2, xi3}.");

    registerOwnership();

    setElement(new LaserRep("LASER"));
}

OpalLaser::OpalLaser(const std::string& name, OpalLaser* parent) : OpalElement(name, parent) {
    setElement(new LaserRep(name));
}

OpalLaser::~OpalLaser() {}

OpalLaser* OpalLaser::clone(const std::string& name) { return new OpalLaser(name, this); }

std::shared_ptr<ElementBase> OpalLaser::build(const ElementValues& values) {
    auto element = std::make_shared<LaserRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalLaser::update() { fill(*getElement(), getValues()); }

void OpalLaser::fill(ElementBase& element, const ElementValues& values) {
    constexpr const char* where = "OpalLaser::update()";

    setCommonFields(element, values);

    auto* laser = dynamic_cast<LaserRep*>(&element);
    if (!laser) {
        return;
    }

    // A LASER that sets none of the required attributes is left unconfigured, as before.
    const bool anyConfigured = values.hasValue("WAVELENGTH") || values.hasValue("PULSEENERGY")
                               || values.hasValue("PULSELENGTH") || values.hasValue("WAISTX")
                               || values.hasValue("WAISTY") || values.hasValue("DIR");
    if (!anyConfigured) {
        return;
    }

    const double wavelength  = getRequiredPositiveReal(values, "WAVELENGTH", where);
    const double pulseEnergy = getRequiredPositiveReal(values, "PULSEENERGY", where);
    const double pulseLength = getRequiredPositiveReal(values, "PULSELENGTH", where);
    const double waistX      = getRequiredPositiveReal(values, "WAISTX", where);
    const double waistY      = getRequiredPositiveReal(values, "WAISTY", where);

    Vector_t<double, 3> direction = getVector3(values, "DIR", where, "the propagation direction");
    const double directionNorm    = std::sqrt(dot(direction, direction));
    if (directionNorm <= 0.0) {
        throw OpalException(where, "\"DIR\" must be a non-zero vector.");
    }
    direction /= directionNorm;

    Vector_t<double, 3> stokes(0.0);
    if (values.hasValue("STOKES")) {
        stokes                   = getVector3(values, "STOKES", where, "the Stokes vector");
        const double stokesNorm2 = dot(stokes, stokes);
        if (stokesNorm2 > 1.0 + 1.0e-12) {
            throw OpalException(where, "\"STOKES\" must satisfy xi1^2 + xi2^2 + xi3^2 <= 1.");
        }
    }

    const double length = values.hasValue("L") ? values.getReal("L") : 0.0;
    if (length < 0.0) {
        throw OpalException(where, "\"L\" must be greater than or equal to 0.");
    }

    laser->getGeometry().setElementLength(length);
    laser->setWavelength(wavelength);
    laser->setPulseEnergy(pulseEnergy);
    laser->setPulseLength(pulseLength);
    laser->setWaistX(waistX);
    laser->setWaistY(waistY);
    laser->setDirection(direction);
    laser->setStokes(stokes);
}
