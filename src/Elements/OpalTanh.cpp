
#include "Elements/OpalTanh.h"

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform *gmsg;

OpalTanh::OpalTanh() :
    OpalElement(SIZE, "TANH",
             "The \"TANH\" element defines a tanh field fall off for plugging"
             "into analytical field models.") {
    itsAttr[X0] = Attributes::makeReal("X0",
                        "Length of the central region of the end field.");
    itsAttr[LAMBDA] = Attributes::makeReal("LAMBDA",
                        "Exponential decay constant for the fringe region.");
    registerOwnership();
}

void OpalTanh::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    double x0 = Attributes::getReal(itsAttr[X0]);
    double lambda = Attributes::getReal(itsAttr[LAMBDA]);
    auto efm = std::make_shared<endfieldmodel::Tanh>(x0, lambda, 10);

    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    efmMan->setEndFieldModel(getOpalName(), efm);
}

OpalTanh::OpalTanh(const std::string &name, OpalTanh *parent):
    OpalElement(name, parent) {
}

OpalTanh *OpalTanh::clone(const std::string &name) {
    return new OpalTanh(name, this);
}
