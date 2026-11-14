#include "Elements/OpalConstantFocusing.h"

#include <memory>

#include "Attributes/Attributes.h"
#include "BeamlineCore/ConstantFocusingRep.h"
#include "Utilities/OpalException.h"

OpalConstantFocusing::OpalConstantFocusing()
    : OpalElement(
              SIZE, "CONSTANTFOCUSING",
              "The \"CONSTANTFOCUSING\" element applies linear 3D focusing scaled from the "
              "initial space-charge field.") {
    itsAttr[STRENGTH] = Attributes::makeReal(
            "STRENGTH", "Dimensionless multiple of the initial mean space-charge field.", 1.0);
    itsAttr[RADIUS] = Attributes::makeReal("RADIUS", "Initial beam radius in m.", 0.0);

    registerOwnership();
}

OpalConstantFocusing::~OpalConstantFocusing() {}

std::shared_ptr<ElementBase> OpalConstantFocusing::build(const ElementValues& values) {
    auto element = std::make_shared<ConstantFocusingRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalConstantFocusing::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    const double radius = values.getReal("RADIUS");
    if (!values.defaultUsed("RADIUS") && radius <= 0.0) {
        throw OpalException("OpalConstantFocusing::build", "CONSTANTFOCUSING requires RADIUS > 0.");
    }

    auto* rep = dynamic_cast<ConstantFocusingRep*>(&element);
    if (rep != nullptr) {
        rep->getGeometry().setElementLength(values.getReal("L"));
        rep->setStrength(values.getReal("STRENGTH"));
        rep->setRadius(radius);
    }
}
