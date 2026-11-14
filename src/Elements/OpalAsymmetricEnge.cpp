
#include "Elements/OpalAsymmetricEnge.h"
#include "Elements/ElementValues.h"

#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform* gmsg;

OpalAsymmetricEnge::OpalAsymmetricEnge()
    : OpalElement(
              SIZE, "ASYMMETRIC_ENGE",
              "The \"ASYMMETRIC_ENGE\" element defines an enge field fall off for"
              "plugging into analytical field models. The Asymmetric version"
              "has different parameters for the start and end of the field.") {
    itsAttr[START_HALF_LENGTH] = Attributes::makeReal(
            "START_HALF_LENGTH",
            "Offset of the central region of the enge element from the start.");
    itsAttr[START_FRINGE_LENGTH] =
            Attributes::makeReal("START_FRINGE_LENGTH", "E-fold length at the element entrance.");
    itsAttr[START_COEFFICIENTS] = Attributes::makeRealArray(
            "START_COEFFICIENTS",
            "Polynomial coefficients for the Enge function at the element entrance.");
    itsAttr[END_HALF_LENGTH] = Attributes::makeReal(
            "END_HALF_LENGTH",
            "Offset of the central region of the enge function element from the end.");
    itsAttr[END_FRINGE_LENGTH] =
            Attributes::makeReal("END_FRINGE_LENGTH", "E-fold length at the element exit.");
    itsAttr[END_COEFFICIENTS] = Attributes::makeRealArray(
            "END_COEFFICIENTS",
            "Polynomial coefficients for the Enge function at the element exit.");
    registerOwnership();
}

void OpalAsymmetricEnge::registerModel(const ElementValues& values) {
    double x0Start                = values.getReal("START_HALF_LENGTH");
    double lambdaStart            = values.getReal("START_FRINGE_LENGTH");
    std::vector<double> aVecStart = values.getRealArray("START_COEFFICIENTS");
    double x0End                  = values.getReal("END_HALF_LENGTH");
    double lambdaEnd              = values.getReal("END_FRINGE_LENGTH");
    std::vector<double> aVecEnd   = values.getRealArray("END_COEFFICIENTS");

    auto efm = std::make_shared<endfieldmodel::AsymmetricEnge>(
            aVecStart, x0Start, lambdaStart, aVecEnd, x0End, lambdaEnd);
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    efmMan->setEndFieldModel(values.getName(), efm);
}
