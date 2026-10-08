
#include "Elements/OpalEnge.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/EndFieldModel/Enge.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform* gmsg;

OpalEnge::OpalEnge()
    : OpalElement(
              SIZE, "ENGE",
              "The \"ENGE\" element defines an enge field fall off for plugging"
              "into analytical field models.") {
    itsAttr[LENGTH] = Attributes::makeReal("LENGTH", "Length of the central region of the enge element [m].");
    itsAttr[FRINGE_LENGTH]       = Attributes::makeReal("FRINGE_LENGTH", "E-fold length of the fringe field [m].");
    itsAttr[COEFFICIENTS] = Attributes::makeRealArray(
            "COEFFICIENTS", "Polynomial coefficients for the Enge function.");
    registerOwnership();
}

void OpalEnge::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    double x0                = Attributes::getReal(itsAttr[LENGTH])/2.0;
    double lambda            = Attributes::getReal(itsAttr[FRINGE_LENGTH]);
    std::vector<double> aVec = Attributes::getRealArray(itsAttr[COEFFICIENTS]);

    auto efm    = std::make_shared<endfieldmodel::Enge>(aVec, x0, lambda);
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    efmMan->setEndFieldModel(getOpalName(), efm);
}

OpalEnge::OpalEnge(const std::string& name, OpalEnge* parent) : OpalElement(name, parent) {}

OpalEnge* OpalEnge::clone(const std::string& name) { return new OpalEnge(name, this); }
