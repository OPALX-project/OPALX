
#include "Elements/OpalTanh.h"
#include "Elements/ElementValues.h"

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform* gmsg;

OpalTanh::OpalTanh()
    : OpalElement(
              SIZE, "TANH",
              "The \"TANH\" element defines a tanh field fall off for plugging"
              "into analytical field models.") {
    itsAttr[LENGTH] =
            Attributes::makeReal("LENGTH", "Length of the central region of the end field.");
    itsAttr[FRINGE_LENGTH] = Attributes::makeReal(
            "FRINGE_LENGTH", "Exponential decay constant for the fringe region.");
    registerOwnership();
}

void OpalTanh::registerModel(const ElementValues& values) {
    double x0     = values.getReal("LENGTH") / 2.0;
    double lambda = values.getReal("FRINGE_LENGTH");
    auto efm      = std::make_shared<endfieldmodel::Tanh>(x0, lambda, 10);

    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    efmMan->setEndFieldModel(values.getName(), efm);
}
