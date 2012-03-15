// ------------------------------------------------------------------------
// $RCSfile: WhileStatement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: WhileStatement
//   Representation for MAD WHILE statement.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "MadParser/WhileStatement.h"
#include "AbstractObjects/MadData.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "MadParser/CompoundStatement.h"
#include "Parser/Parser.h"
#include "Parser/Token.h"
#include "Parser/TokenStream.h"
#include "Utilities/MadException.h"


// class WhileStatement
//   Statement of the form "WHILE ( <condition> ) <statement>".
// ------------------------------------------------------------------------

WhileStatement::WhileStatement(const Parser &parser, TokenStream &is):
  Statement("", 0), while_block(0)
{
  Token key = is.readToken();
  Token token = is.readToken();

  if (key.isKey("WHILE") && token.isDel('(')) {
    int level = 1;
    append(token);
    token = is.readToken();

    while (! token.isEOF()) {
      append(token);

      if (token.isDel('(')) {
	level++;
      } else if (token.isDel(')')) {
	level--;
	if (level == 0) break;
      }

      token = is.readToken();
    }

    while_block = parser.readStatement(&is);
  } else {
    throw MadException("WhileStatement::WhileStatement()",
		       "Invalid \"WHILE\" statement.");
  }
}


WhileStatement::~WhileStatement()
{
  if (while_block != 0) delete while_block;
}


void WhileStatement::execute(const Parser &parser)
{
  curr = tokens.begin();
  keep = ++curr;
  Attribute condition = Attributes::makeBool("WHILE()", "while condition");

  try {
    condition.parse(*this, false);
    MAD.update();

    while (Attributes::getBool(condition)) {
      while_block->execute(parser);
      MAD.update();
    }
  } catch (...) {
    throw MadException("WhileStatement::execute()",
		       "Invalid WHILE condition.");
  }
}
