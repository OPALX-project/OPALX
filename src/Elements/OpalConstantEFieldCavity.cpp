#include "Elements/OpalConstantEFieldCavity.h"

#include <memory>
#include "AbsBeamline/ConstantEFieldCavity.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

OpalConstantEFieldCavity::OpalConstantEFieldCavity()
    : OpalElement(
              SIZE, "CONSTANTEFIELDCAVITY",
              "The \"CONSTANTEFIELDCAVITY\" element defines a constant accelerating "
              "electric field in z-direction.") {
    itsAttr[EX] = Attributes::makeReal("EX", "Electric field strength Ex in MV/m", 0.0);
    itsAttr[EY] = Attributes::makeReal("EY", "Electric field strength Ey in MV/m", 0.0);
    itsAttr[EZ] = Attributes::makeReal("EZ", "Electric field strength Ez in MV/m", 0.0);

    registerOwnership();
}

OpalConstantEFieldCavity::~OpalConstantEFieldCavity() {}

std::shared_ptr<ElementBase> OpalConstantEFieldCavity::build(
        const ObjectValues& values, double /*p0*/) {
    auto element = std::make_shared<ConstantEFieldCavity>(values.getName());
    fill(*element, values);
    return element;
}

void OpalConstantEFieldCavity::fill(ConstantEFieldCavity& element, const ObjectValues& values) {
    setCommonFields(element, values);

    double length = values.getReal("L");
    double exMV   = values.getReal("EX");
    double eyMV   = values.getReal("EY");
    double ezMV   = values.getReal("EZ");
    element.getGeometry().setElementLength(length);
    element.setEx(exMV * Units::MVpm2Vpm);
    element.setEy(eyMV * Units::MVpm2Vpm);
    element.setEz(ezMV * Units::MVpm2Vpm);
}
