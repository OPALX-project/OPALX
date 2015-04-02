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
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "MadParser/CompoundStatement.h"
#include "AbstractObjects/MadData.h"
#include "Parser/Parser.h"


// class CompoundStatement
// ------------------------------------------------------------------------

CompoundStatement::CompoundStatement(TokenStream &is):
  Statement("", 0), tokens(0)
{
  Token token = is.readToken();
  buffer_name = token.getFile();
  stat_line = token.getLine();
  tokens = new MacroStream(token.getFile());

  // Skip opening brace.
  int level = 1;
  token = is.readToken();

  while (! token.isEOF()) {
    if (token.isDel('{')) {
      ++level;
    } else if (token.isDel('}')) {
      if (--level == 0) return;
    }

    tokens->append(token);
    token = is.readToken();
  }
}


CompoundStatement::~CompoundStatement()
{}


void CompoundStatement::execute(const Parser &parser)
{
  tokens->start();
  parser.run(&*tokens);
}
