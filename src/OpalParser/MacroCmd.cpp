// ------------------------------------------------------------------------
// $RCSfile: MacroCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MacroCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/MacroCmd.h"
#include "AbstractObjects/OpalData.h"
#include "OpalParser/OpalParser.h"
#include "Parser/Statement.h"
#include "Utilities/ParseError.h"
#include <vector>
#include <cassert>

// Class MacroCmd
// ------------------------------------------------------------------------

MacroCmd::MacroCmd():
    Macro(0u, "MACRO",
          "A \"MACRO\" command defines a subroutine:\n"
          "\t<name>(<arguments>):MACRO{<body>}"),
    body(0)
{}


MacroCmd::MacroCmd(const std::string &name, MacroCmd *parent):
    Macro(name, parent), body() {
    body = new MacroStream(name);
}


MacroCmd::~MacroCmd()
{}


void MacroCmd::execute() {
    body->start();
    itsParser->run(&*body);
}


Object *MacroCmd::makeInstance
(const std::string &name, Statement &statement, const Parser *parser) {
    parseActuals(statement);

    // Check for consistency in argument number.
    if(formals.size() != actuals.size()) {
        throw ParseError("MacroCmd::makeInstance()",
                         "Inconsistent number of macro arguments.");
    }

    // Substitute the actual arguments.
    MacroCmd *macro = new MacroCmd(name, this);
    macro->itsParser = parser;
    body->start();
    Token token = body->readToken();

    while(! token.isEOF()) {
        bool found = false;

        if(token.isWord()) {
            std::string word = token.getWord();

            for(std::vector<std::string>::size_type i = 0;
                i < formals.size(); i++) {
                if(word == formals[i]) {
                    std::vector<Token> act = actuals[i];
                    for(std::vector<Token>::iterator t = act.begin();
                        t != act.end(); t++) {
                        macro->body->append(*t);
                    }
                    found = true;
                    break;
                }
            }
        }

        if(! found) macro->body->append(token);
        token = body->readToken();
    }

    return macro;
}


Object *MacroCmd::makeTemplate
(const std::string &name, TokenStream &, Statement &statement) {
    MacroCmd *macro = new MacroCmd(name, this);
    macro->parseFormals(statement);

    // Parse macro body->
    assert(statement.keyword("MACRO"));
    Token token;

    if(statement.delimiter('{')) {
        int level = 1;
        while(true) {
            if(statement.atEnd()) {
                throw ParseError("MacroCmd::makeTemplate()",
                                 "MACRO body is not closed.");
            } else {
                token = statement.getCurrent();

                if(token.isDel('{')) {
                    ++level;
                } else if(token.isDel('}')) {
                    if(--level == 0) break;
                }

                macro->body->append(token);
            }
        }
    } else {
        throw ParseError("MacroCmd::makeTemplate()",
                         "Missing MACRO body, should be \"{...}\".");
    }

    return macro;
}
