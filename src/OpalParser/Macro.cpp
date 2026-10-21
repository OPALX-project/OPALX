// ------------------------------------------------------------------------
// $RCSfile: Macro.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Macro/
//   Abstract base class for OPAL macro-like commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/Macro.h"
#include <algorithm>
#include "AbstractObjects/Expressions.h"
#include "OpalParser/Statement.h"
#include "Utilities/ParseError.h"

// Class Macro
// ------------------------------------------------------------------------

Macro::Macro(int size, const char* name, const char* help)
    : Object(size, name, help), formals(), actuals() {}

Macro::Macro(const std::string& name, Object* parent)
    : Object(name, parent), formals(), actuals() {}

Macro::~Macro() {}

Macro* Macro::clone(const std::string&) {
    throw ParseError(
            "Macro::clone()",
            "You cannot execute the command \"" + getOpalName() + "\" without parameters.");
}

const std::string Macro::getCategory() const { return "MACRO"; }

bool Macro::shouldTrace() const { return false; }

bool Macro::shouldUpdate() const { return false; }

void Macro::parseActuals(Statement& stat) {
    // We start after the opening '('.
    actuals.clear();

    if (!stat.delimiter(')')) {
        std::vector<char> closing{')'};
        std::vector<Token> act;

        while (true) {
            if (stat.atEnd()) {
                throw ParseError("Macro::parseActuals()", "MACRO argument list is not closed.");
            }
            Token token = stat.getCurrent();

            if (token.isDel(',') && closing.size() == 1) {
                if (act.empty()) {
                    throw ParseError("Macro::parseActuals()", "Empty MACRO argument.");
                }
                actuals.push_back(act);
                act.clear();
                continue;
            } else if (token.isDel('(')) {
                closing.push_back(')');
            } else if (token.isDel('[')) {
                closing.push_back(']');
            } else if (token.isDel('{')) {
                closing.push_back('}');
            } else if (token.isDel(')') || token.isDel(']') || token.isDel('}')) {
                if (!token.isDel(closing.back())) {
                    throw ParseError(
                            "Macro::parseActuals()", "Mismatched MACRO argument delimiter.");
                }
                closing.pop_back();
                if (closing.empty()) {
                    if (act.empty()) {
                        throw ParseError("Macro::parseActuals()", "Empty MACRO argument.");
                    }
                    actuals.push_back(act);
                    break;
                }
            }

            // Nested delimiters and their commas belong to the current argument.
            act.push_back(token);
        }
    }
}

void Macro::parseFormals(Statement& stat) {
    // We start after the opening '('.
    formals.clear();
    if (!stat.delimiter(')')) {
        do {
            std::string form = Expressions::parseString(stat, "Expected formal argument name.");
            if (std::find(formals.begin(), formals.end(), form) != formals.end()) {
                throw ParseError(
                        "Macro::parseFormals()", "Duplicate MACRO formal argument: " + form);
            }
            formals.push_back(form);
        } while (stat.delimiter(','));

        Expressions::parseDelimiter(stat, ')');
    }

    Expressions::parseDelimiter(stat, ':');
}
