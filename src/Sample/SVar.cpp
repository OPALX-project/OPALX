#include "Sample/SVar.h"
#include "Attributes/Attributes.h"

namespace {
    enum {
        VARIABLE,
        LOWERBOUND,
        UPPERBOUND,
        TYPE,
        SIZE
    };
}

SVar::SVar():
    Definition(SIZE, "SVAR", "The SVAR statement defines a variable for sampling")
{
    itsAttr[VARIABLE] = Attributes::makeString("VARIABLE",
                                               "Variable name that should be varied during sampling");
    
    itsAttr[LOWERBOUND] = Attributes::makeReal("LOWERBOUND",
                                               "Lower limit of the range of values that the variable should assume");
    
    itsAttr[UPPERBOUND] = Attributes::makeReal("UPPERBOUND",
                                               "Upper limit of the range of values that the variable should assume");
    
    itsAttr[TYPE] = Attributes::makeString("TYPE",
                                           "Sampling type: SEQUENCE, QUADRATURE, UNIFORM_INT, UNIFORM_DOUBLE");
    
    registerOwnership(AttributeHandler::STATEMENT);
}

SVar::SVar(const std::string &name, SVar *parent):
    Definition(name, parent)
{ }

SVar::~SVar()
{ }

void SVar::execute() {

}

std::string SVar::getVariable() const {
    return Attributes::getString(itsAttr[VARIABLE]);
}

double SVar::getLowerBound() const {
    return Attributes::getReal(itsAttr[LOWERBOUND]);
}

double SVar::getUpperBound() const {
    return Attributes::getReal(itsAttr[UPPERBOUND]);
}

std::string SVar::getType() const {
    return Attributes::getString(itsAttr[TYPE]);
}
