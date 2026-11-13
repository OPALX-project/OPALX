#include "Elements/OpalConstantEFieldCavity.h"

#include <memory>
#include "Attributes/Attributes.h"
#include "BeamlineCore/ConstantEFieldCavityRep.h"
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

    setElement(new ConstantEFieldCavityRep("CONSTANTEFIELDCAVITY"));
}

OpalConstantEFieldCavity::OpalConstantEFieldCavity(
        const std::string& name, OpalConstantEFieldCavity* parent)
    : OpalElement(name, parent) {
    setElement(new ConstantEFieldCavityRep(name));
}

OpalConstantEFieldCavity::~OpalConstantEFieldCavity() {}

OpalConstantEFieldCavity* OpalConstantEFieldCavity::clone(const std::string& name) {
    return new OpalConstantEFieldCavity(name, this);
}

std::shared_ptr<ElementBase> OpalConstantEFieldCavity::build(const ElementValues& values) {
    auto element = std::make_shared<ConstantEFieldCavityRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalConstantEFieldCavity::update() { fill(*getElement(), getValues()); }

void OpalConstantEFieldCavity::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    ConstantEFieldCavityRep* rep = dynamic_cast<ConstantEFieldCavityRep*>(&element);
    if (rep) {
        double length = values.getReal("L");
        double exMV   = values.getReal("EX");
        double eyMV   = values.getReal("EY");
        double ezMV   = values.getReal("EZ");
        rep->getGeometry().setElementLength(length);
        rep->setEx(exMV * Units::MVpm2Vpm);
        rep->setEy(eyMV * Units::MVpm2Vpm);
        rep->setEz(ezMV * Units::MVpm2Vpm);
    }
}
