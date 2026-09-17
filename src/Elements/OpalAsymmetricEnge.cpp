
#include "Elements/OpalAsymmetricEnge.h"

#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform *gmsg;

OpalAsymmetricEnge::OpalAsymmetricEnge() :
    OpalElement(SIZE, "ASYMMETRIC_ENGE",
             "The \"ASYMMETRIC_ENGE\" element defines an enge field fall off for"
             "plugging into analytical field models. The Asymmetric version"
             "has different parameters for the start and end of the field.") {
    itsAttr[X0_START] = Attributes::makeReal("X0_START",
      "Offset of the central region of the enge element from the start.");
    itsAttr[LAMBDA_START] = Attributes::makeReal("LAMBDA_START",
      "Scales the field rise at the element entrance.");
    itsAttr[COEFFICIENTS_START] = Attributes::makeRealArray("COEFFICIENTS_START",
      "Polynomial coefficients for the Enge function at the element entrance.");
    itsAttr[X0_END] = Attributes::makeReal("X0_END",
      "Offset of the central region of the enge function element from the end.");
    itsAttr[LAMBDA_END] = Attributes::makeReal("LAMBDA_END",
      "Scales the field rise at the element exit.");
    itsAttr[COEFFICIENTS_END] = Attributes::makeRealArray("COEFFICIENTS_END",
      "Polynomial coefficients for the Enge function at the element exit.");
    registerOwnership();
}

void OpalAsymmetricEnge::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    double x0Start = Attributes::getReal(itsAttr[X0_START]);
    double lambdaStart = Attributes::getReal(itsAttr[LAMBDA_START]);
    std::vector<double> aVecStart =
              Attributes::getRealArray(itsAttr[COEFFICIENTS_START]);
    double x0End = Attributes::getReal(itsAttr[X0_END]);
    double lambdaEnd = Attributes::getReal(itsAttr[LAMBDA_END]);
    std::vector<double> aVecEnd =
              Attributes::getRealArray(itsAttr[COEFFICIENTS_END]);

    endfieldmodel::EndFieldModel::setEndFieldModel(getOpalName(),
             std::make_shared<endfieldmodel::AsymmetricEnge>(aVecStart,
                                                             x0Start,
                                                             lambdaStart,
                                                             aVecEnd,
                                                             x0End,
                                                             lambdaEnd));

}

OpalAsymmetricEnge::OpalAsymmetricEnge(const std::string &name,
                                       OpalAsymmetricEnge *parent):
    OpalElement(name, parent) {
}


OpalAsymmetricEnge *OpalAsymmetricEnge::clone(const std::string &name) {
    return new OpalAsymmetricEnge(name, this);
}

OpalAsymmetricEnge::~OpalAsymmetricEnge() {}

