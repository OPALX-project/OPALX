// ------------------------------------------------------------------------
// $RCSfile: ClassicParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ClassicParser
//   This is the default parser for the CLASSIC standard input language.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Parser/ClassicParser.h"
#include "Parser/SimpleStatement.h"
#include "Parser/TokenStream.h"
#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/ElementBase.h"
#include "Construction/ElementFactory.h"
#include "Parser/Statement.h"
#include "Utilities/ClassicException.h"
#include "Utilities/DivideError.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include "Utilities/ParseError.h"

#include "BeamlineCore/CorrectorRep.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/MarkerRep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/Octupole.h"
#include "BeamlineCore/Quadrupole.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/RFCavityRep.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "BeamlineCore/SBendRep.h"
#include "BeamlineCore/SeparatorRep.h"
#include "BeamlineCore/Sextupole.h"
#include "BeamlineCore/SkewOctupole.h"
#include "BeamlineCore/SkewQuadrupole.h"
#include "BeamlineCore/SkewSextupole.h"
#include "BeamlineCore/SolenoidRep.h"
#include "BeamlineCore/XCorrectorRep.h"
#include "BeamlineCore/XMonitorRep.h"
#include "BeamlineCore/YCorrectorRep.h"
#include "BeamlineCore/YMonitorRep.h"
#include "Beamlines/ElmPtr.h"

#include <cerrno>
#include <cmath>
#include <exception>
#include <new>
#include <iostream>

#include <Ippl.h>

extern Inform *gmsg;



namespace {

  // The factory object.
  ElementFactory factory;

  // Table of built-in functions.
  struct FunctionEntry {
    const char *name;
    int   arguments;
  };

  static FunctionEntry functionTable[] = {
    // One argument.
    { "SIGN",  1 },  { "FRAC",  1 },  { "INT",   1 },  { "SQRT",  1 },
    { "LOG",   1 },  { "EXP",   1 },  { "SIN",   1 },  { "COS",   1 },
    { "ABS",   1 },  { "TAN",   1 },  { "ASIN",  1 },  { "ACOS",  1 },
    { "ATAN",  1 },
    // Two arguments.
    { "POW",   2 },  { "ATAN2", 2 },  { "MAX",   2 },  { "MIN",   2 }
  };

  // Table sizes.
  int const functionSize = sizeof(functionTable) / sizeof(char*);

  // Argument values to built-in function.
  double args[2];

  // Map of defined constants.
  AttributeSet constants;
}


// Class ClassicParser
// ------------------------------------------------------------------------


ClassicParser::ClassicParser()
{
  factory.define(new CorrectorRep("CORRECTOR"));
  factory.define(new DriftRep("DRIFT"));
  factory.define(new MultipoleRep("MULTIPOLE"));
  factory.define(new MarkerRep("MARKER"));
  factory.define(new MonitorRep("MONITOR"));
  factory.define(new Octupole("OCTUPOLE"));
  factory.define(new Quadrupole("QUADRUPOLE"));
  factory.define(new RBendRep("RBEND"));
  factory.define(new RFCavityRep("RFCAVITY"));
  factory.define(new TravelingWaveRep("TRAVELINGWAVE"));
  factory.define(new SBendRep("SBEND"));
  factory.define(new SeparatorRep("ELSEPARATOR"));
  factory.define(new Sextupole("SEXTUPOLE"));
  factory.define(new SkewOctupole("SKEWOCTUPOLE"));
  factory.define(new SkewQuadrupole("SKEWQUADRUPOLE"));
  factory.define(new SkewSextupole("SKEWSEXTUPOLE"));
  factory.define(new SolenoidRep("SOLENOID"));
  factory.define(new XCorrectorRep("CORRECTOR"));
  factory.define(new XMonitorRep("XMONITOR"));
  factory.define(new YCorrectorRep("CORRECTOR"));
  factory.define(new YMonitorRep("YMONITOR"));
}


ClassicParser::~ClassicParser()
{}


void ClassicParser::parse(Statement &statement) const
{
  try {
    string name;

    if (! statement.word(name)) {
      throw ParseError("ClassicParser::parse()",
		       "Definition should begin with a name.");
    }

    if (statement.delimiter('=')  ||  statement.delimiter(":=")) {
      // Parameter definition: <name> [:]= <expression>
      double value = parseExpression(statement);

      if (constants.hasAttribute(name)) {
	throw ParseError("ClassicParser::parse()",
			 "Parameter already defined.");
      } else {
	constants.setAttribute(name, value);
      }
    } else {
      // Statement must be an element or beamline definition.
      string type;

      if (! statement.delimiter(':')  ||  ! statement.word(type)) {
	throw ParseError("ClassicParser::parse()",
			 "Expected: <name>:<type>{,<keyword>=<expression>}.");
      }

      statement.mark();

      if (type == "LINE") {
	if (statement.delimiter('=')) {
	  // "<name> : LINE = (<list>)"
	  SimpleBeamline *line = parseLine(statement);
	  line->setName(name);
	  factory.define(line);
	} else {
	  throw ParseError("ClassicParser::parse()", "Equals sign missing.");
	}
      } else if (statement.delimiter(',') || statement.atEnd()) {
	// Element definition: "<name> : <class> {, <attributes> }"
	statement.restore();
	AttributeSet set;
	
	while (statement.delimiter(',')) {
	  string attrName;
	
	  if (statement.word(attrName) && statement.delimiter('=')) {
	    double value = parseExpression(statement);

	    if (set.hasAttribute(attrName)) {
	      throw ParseError("ClassicParser::parse()",
			       "Attribute \"" + attrName +
			       "\" is already defined.");
	    } else {
	      set.setAttribute(attrName, value);
	    }
	  } else {
	    throw ParseError("ClassicParser::parse()",
			     "Expected: <attributeName>=<expression>.");
	  }
	}

	if (! factory.makeElement(type, name, set)) {
	  throw ParseError("ClassicParser::parse()",
			   "Unknown element keyword \"" + type + "\".");
	}
	
	if (! statement.atEnd()) {
	  throw ParseError("ClassicParser::parse()",
			   "Could not consume the entire statement.");
	}
      }
    }
  } catch (ParseError &ex) {
    *gmsg << endl << "*** Parse error detected by method \"" << ex.where()
	 << "\"" << endl;
    statement.printWhere(true);
    statement.print();
    *gmsg << ex.what() << endl << endl;
  } catch (ClassicException &ex) {
    *gmsg << endl << "*** Error detected by method \"" << ex.where()
	 << "\"" << endl;
    statement.printWhere(false);
    statement.print();
    *gmsg << ex.what() << endl << endl;
  } catch (std::bad_alloc &) {
    *gmsg << endl << "*** Error:" << endl;
    statement.printWhere(false);
    statement.print();
    *gmsg << "Sorry, virtual memory exhausted." << endl << endl;
  } catch (std::exception &ex) {
    *gmsg << endl << "*** Error:" << endl;
    statement.printWhere(false);
    statement.print();
    *gmsg << "Internal CLASSIC error " << ex.what() << endl << endl;
  } catch (...) {
    *gmsg << endl << "*** Error:" << endl;
    statement.printWhere(false);
    statement.print();
    *gmsg << "Unexpected exception caught." << endl << endl;
    abort();
  }
}


double ClassicParser::evaluate(int oper, double args[]) const
{
  double result;
  errno = 0;

  // Switch on operation code.
  switch (oper) {
  case 0:    // Sign of argument.
    result = (args[0] >= 0.0 ? 1.0 : -1.0);
    break;
  case 1:    // Fractional part of argument.
    result = args[0] - double(int(args[0]));
    break;
  case 2:    // Integer part of argument.
    result = double(int(args[0]));
    break;
  case 3:    // Square root.
    result = sqrt(args[0]);
    break;
  case 4:    // Natural logarithm.
    result = log(args[0]);
    break;
  case 5:    // Exponential.
    result = exp(args[0]);
    break;
  case 6:    // Trigonometric sine.
    result = sin(args[0]);
    break;
  case 7:    // Trigonometric cosine.
    result = cos(args[0]);
    break;
  case 8:    // Absolute value.
    result = abs(args[0]);
    break;
  case 9:    // Trigonometric tangent.
    result = tan(args[0]);
    break;
  case 10:   // Arc sine.
    result = asin(args[0]);
    break;
  case 11:   // Arc cosine.
    result = acos(args[0]);
    break;
  case 12:   // Arc tangent.
    result = atan(args[0]);
    break;
  case 13:   // Power.
    result = pow(args[0], args[1]);
    break;
  case 14:   // Arc tangent of type atan2(x,y).
    result = atan2(args[0], args[1]);
    break;
  case 15:   // Maximum value.
    result = (args[0] > args[1]) ? args[0] : args[1];
    break;
  case 16:   // Minimum value.
    result = (args[0] < args[1]) ? args[0] : args[1];
    break;
  default:   // Illegal.
    return 0.0;
  }

  // Check for correct operation.
  switch (errno) {
  case EDOM:
    throw DomainError("ClassicParser::evaluate()");
  case ERANGE:
    // Ignore underflow.
    if (result == 0.0) return result;
    throw OverflowError("ClassicParser::evaluate()");
  default:
    return result;
  }
}


SimpleBeamline *ClassicParser::parseLine
(Statement &statement) const
{
  SimpleBeamline *line = 0;

  try {
    line = new SimpleBeamline();

    if (statement.delimiter('(')) {
      do {
	ElementBase *element = 0;
	string name;
	statement.mark();
	
	if (statement.word(name)) {
	  element = factory.find(name);
	} else if (statement.delimiter('(')) {
	  statement.restore();
	  element = parseLine(statement);
	} else {
	  throw ParseError("ClassicParser::parseLine()",
			   "Element or line name missing.");
	}
	
	line->append(ElmPtr(element));
      } while (statement.delimiter(','));

      if (! statement.delimiter(')')) {
	throw ParseError("ClassicParser::parseLine()",
			 "Missing comma ',' or right parenthesis ')'");
      }
    }

    return line;
  } catch (...) {
    delete line;
    return 0;
  }
}


double ClassicParser::parseExpression
(Statement &statement) const
{
  double result = 0.0;

  // Unary operator and leading term.
  if (statement.delimiter('+')) {
    result = parseTerm(statement);
  } else if (statement.delimiter('-')) {
    result = - parseTerm(statement);
  } else {
    result = parseTerm(statement);
  }

  // Subsequent add operators and terms.
  while (true) {
    if (statement.delimiter('+')) {
      result += parseTerm(statement);
    } else if (statement.delimiter('-')) {
      result -= parseTerm(statement);
    } else {
      break;
    }
  }

  return result;
}


double ClassicParser::parseFactor(Statement &statement) const
{
  // Leading expression primary.
  double result = parsePrimary(statement);

  if (statement.delimiter('^')) {
    double exponent = parsePrimary(statement);
    result = pow(result, exponent);
  }

  return result;
}


double ClassicParser::parsePrimary(Statement &statement) const
{
  double result = 0.0;
  string name;

  if (statement.delimiter('(')) {
    result = parseExpression(statement);

    if (! statement.delimiter(')')) {
      throw ParseError("ClassicParser::parse()",
		       "Missing right parenthesis.");
    }
  } else if (statement.word(name)) {
    if (statement.delimiter('(')) {
      // "<name>(...)" expected, i. e. function designator.
      for (int oper = 0; oper < functionSize; ++oper) {
      	if (name == functionTable[oper].name) {
	  // "name" is known function; get first argument.
	  args[0] = parseExpression(statement);

	  // Get subsequent arguments.
	  for (int i = 1; i < functionTable[oper].arguments; ++i) {
	    if (! statement.delimiter(',')) {
	      throw ParseError("ClassicParser::parse()",
			       "Missing \",\" in function argument list.");
	    } else {
	      args[i] = parseExpression(statement);
	    }
	  }

	  if (! statement.delimiter(')')) {
	    throw ParseError("ClassicParser::parse()",
			     "Missing \")\" after function argument list.");
	  } else {
	    return evaluate(oper, args);
	  }

	  break;
	}
      }

      // Unknown function: result is still zero.
      throw ParseError("ClassicParser::parsePrimary()",
		       "Unknown function name \"" + name + "\".");
    } else if (! constants.hasAttribute(name)){
      // "name" should refer to a named constant.
      throw ParseError("ClassicParser::parsePrimary()",
		       "Unknown constant name \"" + name + "\".");
    }

    result = constants.getAttribute(name);
  } else if (! statement.real(result)) {
    // Primary should have been a number.
    throw ParseError("ClassicParser::parsePrimary()",
		     "Expression primary expected.");
  }

  return result;
}


double ClassicParser::parseTerm(Statement &statement) const
{
  // First factor.
  double result = parseFactor(statement);

  // Subsequence multiply operators and factors.
  while (true) {
    if (statement.delimiter('*')) {
      result *= parseFactor(statement);
    } else if (statement.delimiter('/')) {
      double arg = parseFactor(statement);

      if (arg == 0.0) {
	throw DivideError("ClassicParser::parseTerm()");
      } else {
	result /= result;
      }
    } else {
      break;
    }
  }

  return result;
}


Statement *ClassicParser::readStatement(TokenStream *is) const
{
  SimpleStatement *statement = 0;
  Token token = is->readToken();

  // Must be ordinary (not STOP) statement.
  if (! token.isEOF()  &&  ! token.isKey("STOP")) {
    try {
      statement = new SimpleStatement(is->getName(), is->getLine());

      while (! token.isDel(';')) {
	statement->append(token);
	token = is->readToken();
	if (token.isEOF()) break;
      }
    } catch (...) {
      // In case of error simulate end of file.
      delete statement;
      statement = 0;
    }
  }

  return statement;
}


void ClassicParser::run(TokenStream *is) const
{
  while (Statement *statement = readStatement(is)) {
    statement->start();
    parse(*statement);
    delete statement;
  }
}
