//
// Class MacroCmd
//
//   This class parses the MACRO command.
//   Encapsulate the buffer for the ``archetypes'' of all macros.
//   The macro is stored as tokens. For execution, first the
//   parameters are replaced, then the resulting stream is sent to the parser.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "OpalParser/MacroCmd.h"

#include "Utility/PAssert.h"

#include <memory>
#include "OpalParser/MacroStream.h"
#include "OpalParser/Parser.h"
#include "OpalParser/Statement.h"
#include "Utilities/ParseError.h"

MacroCmd::MacroCmd()
    : Macro(0u, "MACRO",
            "A \"MACRO\" command defines a subroutine:\n"
            "\t<name>(<arguments>):MACRO{<body>}"),
      body() {}

MacroCmd::MacroCmd(const std::string& name, MacroCmd* parent) : Macro(name, parent), body() {}

MacroCmd::~MacroCmd() {}

void MacroCmd::execute() {
    if (!itsParser) {
        throw ParseError("MacroCmd::execute()", "A MACRO template cannot be executed directly.");
    }
    auto stream = std::make_unique<MacroStream>(getOpalName());
    for (Token token : body)
        stream->append(token);
    stream->start();
    // Parser::run takes ownership; the stored tokens remain reusable.
    itsParser->run(stream.release());
}

Object* MacroCmd::makeInstance(
        const std::string& name, Statement& statement, const Parser* parser) {
    parseActuals(statement);

    // Check for consistency in argument number.
    if (formals.size() != actuals.size()) {
        throw ParseError("MacroCmd::makeInstance()", "Inconsistent number of macro arguments.");
    }

    // Substitute the actual arguments.
    auto macro       = std::make_unique<MacroCmd>(name, this);
    macro->itsParser = parser;
    for (const Token& token : body) {
        bool found = false;

        if (token.isWord()) {
            std::string word = token.getWord();

            for (std::vector<std::string>::size_type i = 0; i < formals.size(); i++) {
                if (word == formals[i]) {
                    macro->body.insert(macro->body.end(), actuals[i].begin(), actuals[i].end());
                    found = true;
                    break;
                }
            }
        }

        if (!found) macro->body.push_back(token);
    }

    return macro.release();
}

Object* MacroCmd::makeTemplate(const std::string& name, TokenStream&, Statement& statement) {
    auto macro = std::make_unique<MacroCmd>(name, this);
    macro->parseFormals(statement);

    // Parse macro body->
    [[maybe_unused]] bool isMacro = statement.keyword("MACRO");
    PAssert(isMacro);
    Token token;

    if (statement.delimiter('{')) {
        int level = 1;
        while (true) {
            if (statement.atEnd()) {
                throw ParseError("MacroCmd::makeTemplate()", "MACRO body is not closed.");
            } else {
                token = statement.getCurrent();

                if (token.isDel('{')) {
                    ++level;
                } else if (token.isDel('}')) {
                    if (--level == 0) break;
                }

                macro->body.push_back(token);
            }
        }
    } else {
        throw ParseError("MacroCmd::makeTemplate()", "Missing MACRO body, should be \"{...}\".");
    }

    return macro.release();
}
