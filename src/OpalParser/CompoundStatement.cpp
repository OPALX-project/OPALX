// ------------------------------------------------------------------------
// $RCSfile: CompoundStatement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CompoundStatement
//   Compound statement of the form "{ <statement>; <statement>; }"
//   (may be used for FOR, IF, or WHILE blocks).
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/CompoundStatement.h"
#include <memory>
#include "OpalParser/MacroStream.h"
#include "OpalParser/Parser.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"

// class CompoundStatement
// ------------------------------------------------------------------------

CompoundStatement::CompoundStatement(TokenStream& is) : Statement("", 0), tokens() {
    Token token = is.readToken();
    buffer_name = token.getFile();
    stat_line   = token.getLine();

    // Skip opening brace.
    int level = 1;
    token     = is.readToken();

    while (!token.isEOF()) {
        if (token.isDel('{')) {
            ++level;
        } else if (token.isDel('}')) {
            if (--level == 0) return;
        }

        tokens.push_back(token);
        token = is.readToken();
    }
    if (Options::parseOnly)
        throw ParseError("CompoundStatement", "Compound statement is missing '}'.");
}

CompoundStatement::~CompoundStatement() {}

void CompoundStatement::execute(const Parser& parser) {
    auto stream = std::make_unique<MacroStream>(buffer_name);
    for (Token token : tokens)
        stream->append(token);
    stream->start();
    parser.run(stream.release());
}
