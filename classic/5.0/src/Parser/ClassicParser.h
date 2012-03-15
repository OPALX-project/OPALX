#ifndef CLASSIC_ClassicParser_HH
#define CLASSIC_ClassicParser_HH

// ------------------------------------------------------------------------
// $RCSfile: ClassicParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ClassicParser
//   
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/AttributeSet.h"
#include "Beamlines/SimpleBeamline.h"
#include "Parser/Parser.h"
#include <string>

using std::string;

class TokenStream;


// Class ClassicParser
// ------------------------------------------------------------------------
/// A parser for the standard input format (SIF) described in
//  {center}
//  D. C. Carey and F. Ch. Iselin:{br}
//  A Standard Input Language for Particle Beam and Accelerator{br}
//  Computer Programs.{br}
//  1984 Snowmass Summer Study.
//  {/center}
//  This is the default parser for the CLASSIC standard input language.
//  The program should declare only one object of class ClassicParser.

class ClassicParser: public Parser {

public:

  ClassicParser();
  virtual ~ClassicParser();

  /// Parse and execute the statement.
  virtual void parse(Statement &stat) const;

  /// Read complete statement from token stream.
  virtual Statement *readStatement(TokenStream *ts) const;

  /// Read statements and parse.
  //  Read one statement at a time on token stream, parse it, and execute it.
  virtual void run(TokenStream *ts) const;

private:

  // Not implemented.
  ClassicParser(const ClassicParser &);
  void operator=(const ClassicParser &);

  // Evaluate built-in function.
  double evaluate(int oper, double args[]) const;

  // Parse routine for beamlines.
  SimpleBeamline *parseLine(Statement &statement) const;

  // Parse routines for expressions.
  double parseExpression(Statement &statement) const;
  double parseFactor(Statement &statement) const;
  double parsePrimary(Statement &statement) const;
  double parseTerm(Statement &statement) const;
};

#endif // CLASSIC_ClassicParser_HH
