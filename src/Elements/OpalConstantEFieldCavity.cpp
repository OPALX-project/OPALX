#include "Elements/OpalConstantEFieldCavity.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ConstantEFieldCavityRep.h"
#include "Physics/Units.h"

OpalConstantEFieldCavity::OpalConstantEFieldCavity()
    : OpalElement(SIZE, "CONSTANTEFIELDCAVITY",
                  "The \"CONSTANTEFIELDCAVITY\" element defines a constant accelerating "
                  "electric field in z-direction.") {
    itsAttr[EX] = Attributes::makeReal("EX", "Electric field strength Ex in MV/m", 0.0);
    itsAttr[EY] = Attributes::makeReal("EY", "Electric field strength Ey in MV/m", 0.0);
    itsAttr[EZ] = Attributes::makeReal("EZ", "Electric field strength Ez in MV/m", 0.0);

    registerOwnership();

    setElement(new ConstantEFieldCavityRep("CONSTANTEFIELDCAVITY"));
}

OpalConstantEFieldCavity::OpalConstantEFieldCavity(const std::string& name,
                                                   OpalConstantEFieldCavity* parent)
    : OpalElement(name, parent) {
    setElement(new ConstantEFieldCavityRep(name));
}

OpalConstantEFieldCavity::~OpalConstantEFieldCavity() {}

OpalConstantEFieldCavity* OpalConstantEFieldCavity::clone(const std::string& name) {
    return new OpalConstantEFieldCavity(name, this);
}

void OpalConstantEFieldCavity::update() {
    OpalElement::update();

    ConstantEFieldCavityRep* rep = dynamic_cast<ConstantEFieldCavityRep*>(getElement());
    if (rep) {
        double length = getLength();
        double exMV   = Attributes::getReal(itsAttr[EX]);
        double eyMV   = Attributes::getReal(itsAttr[EY]);
        double ezMV   = Attributes::getReal(itsAttr[EZ]);
        rep->setElementLength(length);
        rep->setEx(exMV * Units::MVpm2Vpm);
        rep->setEy(eyMV * Units::MVpm2Vpm);
        rep->setEz(ezMV * Units::MVpm2Vpm);
    }

    OpalElement::updateUnknown(getElement());
}

