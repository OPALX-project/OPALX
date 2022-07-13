
#include "Elements/OpalEnge.h"

#include "AbsBeamline/EndFieldModel/Enge.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform *gmsg;

OpalEnge::OpalEnge() :
    OpalElement(SIZE, "ENGE",
             "The \"ENGE\" element defines an enge field fall off for plugging"
             "into analytical field models.") {
    itsAttr[X0] = Attributes::makeReal("X0",
                        "Length of the central region of the enge element.");
    itsAttr[LAMBDA] = Attributes::makeReal("LAMBDA",
                        "Scales the end field fall off.");
    itsAttr[COEFFICIENTS] = Attributes::makeRealArray("COEFFICIENTS",
                        "Polynomial coefficients for the Enge function.");
    registerOwnership();
}

void OpalEnge::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    double x0 = Attributes::getReal(itsAttr[X0]);
    double lambda = Attributes::getReal(itsAttr[LAMBDA]);
    std::vector<double> aVec = Attributes::getRealArray(itsAttr[COEFFICIENTS]);

    endfieldmodel::EndFieldModel::setEndFieldModel(getOpalName(),
                    std::make_shared<endfieldmodel::Enge>(aVec, x0, lambda));
}

OpalEnge::OpalEnge(const std::string &name, OpalEnge *parent):
    OpalElement(name, parent) {
}


OpalEnge *OpalEnge::clone(const std::string &name) {
    return new OpalEnge(name, this);
}
