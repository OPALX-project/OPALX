#include "Elements/OpalConstantEz.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ConstantEzRep.h"
#include "Physics/Units.h"

OpalConstantEz::OpalConstantEz()
    : OpalElement(SIZE, "CONSTANTEZ",
                  "The \"CONSTANTEZ\" element defines a constant accelerating "
                  "electric field in z-direction.") {
    itsAttr[EZ] = Attributes::makeReal("EZ", "Electric field strength in MV/m", 0.0);

    registerOwnership();

    setElement(new ConstantEzRep("CONSTANTEZ"));
}

OpalConstantEz::OpalConstantEz(const std::string& name, OpalConstantEz* parent)
    : OpalElement(name, parent) {
    setElement(new ConstantEzRep(name));
}

OpalConstantEz::~OpalConstantEz() {}

OpalConstantEz* OpalConstantEz::clone(const std::string& name) {
    return new OpalConstantEz(name, this);
}

void OpalConstantEz::update() {
    OpalElement::update();

    ConstantEzRep* rep = dynamic_cast<ConstantEzRep*>(getElement());
    if (rep) {
        double length = getLength();
        double ezMV   = Attributes::getReal(itsAttr[EZ]);
        rep->setElementLength(length);
        rep->setEz(ezMV * Units::MVpm2Vpm);
    }

    OpalElement::updateUnknown(getElement());
}
