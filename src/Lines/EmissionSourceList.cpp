#include "Lines/EmissionSourceList.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "OpalParser/Statement.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"

using namespace Expressions;

EmissionSourceList::EmissionSourceList()
    : Definition(
              0, "EMISSIONSOURCELIST",
              "The EMISSIONSOURCELIST statement defines a list of emission sources.\n"
              "\t<name> : EMISSIONSOURCELIST = (ES1, ES2, ...)") {
    registerOwnership(AttributeHandler::STATEMENT);
}

EmissionSourceList::EmissionSourceList(const std::string& name, EmissionSourceList* parent)
    : Definition(name, parent) {}

EmissionSourceList::~EmissionSourceList() {}

bool EmissionSourceList::canReplaceBy(Object* object) {
    return dynamic_cast<EmissionSourceList*>(object) != nullptr;
}

EmissionSourceList* EmissionSourceList::clone(const std::string& name) {
    return new EmissionSourceList(name, this);
}

void EmissionSourceList::execute() {}

EmissionSourceList* EmissionSourceList::find(const std::string& name) {
    Object* obj = OpalData::getInstance()->find(name);
    auto* esl   = dynamic_cast<EmissionSourceList*>(obj);
    if (esl == nullptr) {
        throw OpalException(
                "EmissionSourceList::find()", "EmissionSourceList \"" + name + "\" not found.");
    }
    return esl;
}

void EmissionSourceList::parse(Statement& stat) {
    parseDelimiter(stat, '=');
    parseDelimiter(stat, '(');
    sources_m.clear();
    parseList(stat);
    parseDelimiter(stat, ')');
}

void EmissionSourceList::parseList(Statement& stat) {
    do {
        std::string name = parseString(stat, "EmissionSource name expected.");
        Object* obj      = OpalData::getInstance()->find(name);
        if (obj == nullptr) {
            throw ParseError(
                    "EmissionSourceList::parseList()",
                    "EmissionSource \"" + name + "\" is undefined.");
        }
        auto* es = dynamic_cast<EmissionSource*>(obj);
        if (es == nullptr) {
            throw ParseError(
                    "EmissionSourceList::parseList()",
                    "Object \"" + name + "\" is not an EmissionSource.");
        }
        sources_m.push_back(es);
    } while (stat.delimiter(','));
}
