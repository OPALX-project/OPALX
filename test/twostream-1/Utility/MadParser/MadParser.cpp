// ------------------------------------------------------------------------
// $RCSfile: MadParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class MadParser:
//   This is the default parser for MAD statements.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:17:27 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "MadParser/MadParser.h"
#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/MadData.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Attributes/Attributes.h"
#include "MadParser/CompoundStatement.h"
#include "MadParser/IfStatement.h"
#include "MadParser/WhileStatement.h"
#include "MemoryManagement/Pointer.h"
#include "Parser/SimpleStatement.h"
#include "Parser/Token.h"
#include "Utilities/MadException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <cassert>
#include <ctime>
#include <exception>
#include <iostream>
#include <new>
using namespace std;

using namespace Expressions;
using std::cerr;
using std::endl;


// Class MadParser
// ------------------------------------------------------------------------

std::vector<Pointer<TokenStream> > MadParser::inputStack;


MadParser::MadParser(): stopFlag(false)
{}


MadParser::~MadParser()
{}


void MadParser::parse(Statement &stat) const
{
  if (stat.keyword("SHARED")) {
    // "SHARED ...": Shared object definition.
    parseDefine(stat);
  } else if (stat.keyword("CONSTANT") || stat.keyword("CONST") ||
	     stat.keyword("BOOL") || stat.keyword("REAL") ||
	     stat.keyword("STRING") || stat.keyword("VECTOR")) {
    // Keywords introducing variable definitions.
    parseAssign(stat);
  } else {
    string name = parseString(stat, "Identifier or keyword expected.");

    if (stat.delimiter('?')) {
      // "<class>?": give help for class.
      printHelp(name);
    } else if (stat.delimiter(':')) {
      // "<object>:<class>...": labelled command.
      parseDefine(stat);
    } else if (stat.delimiter('(')) {
      // "<macro>(...)...": macro definition or call.
      // We are positioned just after the '(' of the argument list.
      parseMacro(name, stat);
    } else if (stat.delimiter(',') || stat.delimiter(';') ||
	       stat.atEnd()) {
      // "<class>" or "<class>,<attributes>": Executable command.
      parseAction(stat);
    } else {
      // Assignment beginning with a name.
      parseAssign(stat);
    }
  }
}


void MadParser::execute(Object *object, const string &name) const
{
  // Trace execution.
  if (Options::trace && object->shouldTrace()) {
    double time = double(clock()) / double(CLOCKS_PER_SEC);
    cerr << "\nBegin execution: \"" << name
	 << "\", CPU time = " << time << " seconds.\n" << endl;
  }

  // Force updating of all attributes which might have been changed.
  if (object->shouldUpdate()) {
    MAD.update();
  }

  // Execute or check the command.
  object->execute();

  // Trace execution.
  if (Options::trace && object->shouldTrace()) {
    double time = double(clock()) / double(CLOCKS_PER_SEC);
    cerr << "\nEnd execution:   \"" << name
	 << "\", CPU time = " << time << " seconds.\n" << endl;
  }
}


Object *MadParser::find(const string &name) const
{
  return MAD.find(name);
}


void MadParser::parseAction(Statement &stat) const
{
  stat.start();
  string cmdName = parseString(stat, "Command name expected");

  if (cmdName == "STOP") {
    stopFlag = true;
  } else if (cmdName == "HELP"  &&  stat.delimiter(',')) {
    cmdName = parseString(stat, "Object name expected");
    printHelp(cmdName);
  } else if (Object *object = find(cmdName)) {
    Object* copy = 0;
    try {
      copy = object->clone("");
      copy->parse(stat);
      parseEnd(stat);
      execute(copy, cmdName);
    } catch (...) {
      delete copy;
      throw;
    }
  } else {
    throw ParseError("MadParser::parseAction()",
		     "Command \"" + cmdName + "\" is unknown.");
  }
}


void MadParser::parseAssign(Statement &stat) const
{
  stat.start();

  // Find various model objects.
  static Object *boolConstant   = MAD.find("BOOL_CONSTANT");
  static Object *realConstant   = MAD.find("REAL_CONSTANT");
  static Object *realVariable   = MAD.find("REAL_VARIABLE");
  static Object *realVector     = MAD.find("REAL_VECTOR");
  static Object *stringConstant = MAD.find("STRING_CONSTANT");

  // Gobble up any prefix.
  int code = 0x00;
  while (true) {
    if (stat.keyword("CONSTANT") || stat.keyword("CONST")) {
      code |= 0x01;
    } else if (stat.keyword("BOOL")) {
      code |= 0x02;
    } else if (stat.keyword("REAL")) {
      code |= 0x04;
    } else if (stat.keyword("STRING")) {
      code |= 0x08;
    } else if (stat.keyword("VECTOR")) {
      code |= 0x10;
    } else {
      break;
    }
  }

  string objName = parseString(stat, "Object name expected.");

  // Test for attribute name.
  Object *object = 0;
  string attrName;

  if (stat.delimiter("->")) {
    // Assignment to object attribute.
    attrName = parseString(stat, "Attribute name expected.");

    if (code != 0) {
      throw ParseError("MadParser::parseAssign()",
		       "Invalid type specification for this value.");
    } else if ((object = MAD.find(objName)) == 0) {
      throw ParseError("MadParser::parseAssign()",
		       "The object \"" + objName + "\" is unknown.");
    }
  } else {
    // Assignment to variable-like object.
    if ((object = MAD.find(objName)) == 0) {
      Object *model = 0;
      switch (code) {
      case 0x01:  // CONSTANT
      case 0x05:  // CONSTANT REAL
	model = realConstant;
	break;
      case 0x02:  // BOOL
      case 0x03:  // BOOL CONSTANT
	model = boolConstant;
	break;
      case 0x00:  // empty <type>.
      case 0x04:  // REAL
	model = realVariable;
	break;
      case 0x10:  // VECTOR
      case 0x11:  // CONSTANT VECTOR
      case 0x14:  // REAL VECTOR
      case 0x15:  // CONSTANT REAL VECTOR
	model = realVector;
	break;
      case 0x08:  // STRING
      case 0x09:  // STRING CONSTANT
	model = stringConstant;
	break;
      default:
	break;
      }

      if (model != 0) {
	object = model->clone(objName);
	MAD.define(object);
      } else {
	throw ParseError("MadParser::parseAssign()", "Invalid <type> field.");
      }
    } else if (object->isTreeMember(realConstant)) {
      throw ParseError("MadParser::parseAssign()",
		       "You cannot redefine the constant \""+objName+"\".");
    }

    attrName = "VALUE";
  }

  // Test for index; it is evaluated immediately.
  int index = 0;

  if (stat.delimiter('[')) {
    index = int(Round(parseRealConst(stat)));
    parseDelimiter(stat, ']');

    if (index <= 0) {
      throw ParseError("Expressions::parseReference()",
		       "Index must be positive.");
    }
  }

  if (object != 0) {
    if (Attribute *attr = object->findAttribute(attrName)) {
      if (stat.delimiter('=') || object->isTreeMember(realConstant)) {
	if (index > 0) {
	  attr->parseComponent(stat, true, index);
	} else {
	  attr->parse(stat, true);
	}
      } else if (stat.delimiter(":=")) {
	if (index > 0) {
	  attr->parseComponent(stat, false, index);
	} else {
	  attr->parse(stat, false);
	}
      }
    } else {
      throw ParseError("MadParser::parseAssign()",
		       "Object \"" + objName + "\" has no attribute \"" +
		       attrName + "\".");
    }

    parseEnd(stat);
    MAD.makeDirty(object);
  }
}


void MadParser::parseDefine(Statement &stat) const
{
  stat.start();
  bool isShared = stat.keyword("SHARED");
  string objName = parseString(stat, "Object name expected.");

  if (stat.delimiter(':')) {
    string clsName = parseString(stat, "Class name expected.");
    Object *classObject = find(clsName);

    if (classObject == 0) {
      throw ParseError("MadParser::parseDefine()",
		       "The object \"" + clsName + "\" is unknown.");
    }

    Object* copy = 0;
    try {
      if (stat.delimiter('(')) {
	// Macro-like objects are always classes, instances never.
	// There is no further check required.
	copy = classObject->makeInstance(objName, stat, this);
      } else {
	copy = classObject->clone(objName);
	copy->parse(stat);
	copy->setShared(isShared);
      }

      parseEnd(stat);
      execute(copy, clsName);
      MAD.define(copy);
    } catch (...) {
      delete copy;
      throw;
    }
  } else {
    // Redefine an object to be a class.
    Object *classObject = find(objName);
    Object *copy = classObject->clone(objName);
    copy->parse(stat);
    copy->setShared(isShared);
  }
}


void MadParser::parseEnd(Statement &stat) const
{
  if (! stat.atEnd()  &&  ! stat.delimiter(';')) {
    throw ParseError("MadParser::parseEnd()",
		     "Syntax error (maybe missing comma or semicolon ? )");
  }
}


void MadParser::parseMacro(const string &macName, Statement &stat) const
{
  // Record the position just after the '(' of the argument list.
  stat.mark();

  // Skip argument list.
  int par_level = 1;
  while (true) {
    if (stat.delimiter('(')) {
      ++par_level;
    } else if (stat.delimiter(')')) {
      if (--par_level == 0) break;
    } else {
      stat.getCurrent();
    }
  }

  if (stat.delimiter(':')) {
    // Macro definition.
    string className = parseString(stat, "Class name expected.");

    if (Object *macro = MAD.find(className)) {
      // Backtrack to first argument.
      stat.restore();

      if (Object *copy =
	  macro->makeTemplate(macName, *inputStack.back(), stat)) {
	MAD.define(copy);
      } else {
	throw ParseError("MadParser::parseMacro()", "Command \"" +
			 macName + "\" cannot be defined with arguments.");
      }
    } else {
      throw ParseError("MadParser::parseMacro()",
		       "Object \"" + className + "\" is unknown.");
    }
  } else {
    // Macro call.
    if (Object *macro = MAD.find(macName)) {
      // Backtrack to first argument.
      stat.restore();
      Object* instance = 0;
      try {
	instance = macro->makeInstance(macName, stat, this);
	execute(instance, macName);
      } catch (...) {
	delete instance;
	throw;
      }
    } else {
      throw ParseError("MadParser::parseMacro()",
		       "Macro \"" + macName + "\" is unknown.");
    }
  }
}


void MadParser::printHelp(const string &cmdName) const
{
  Object *object = find(cmdName);

  if (object == 0) {
    cerr << "\nMadParser::printHelp(): Unknown object \""
	 << cmdName << "\".\n" << endl;
  } else {
    object->printHelp(cerr);
  }
}


void MadParser::parseBracketList(char close, Statement &stat)
{
  Token token = readToken();

  while (! token.isEOF()) {
    stat.append(token);

    if (token.isDel('(')) {
      parseBracketList(')', stat);
    } else if (token.isDel('[')) {
      parseBracketList(']', stat);
    } else if (token.isDel('{')) {
      parseBracketList('}', stat);
    } else if (token.isDel(close)) {
      return;
    }

    token = readToken();
  }
}


void MadParser::parseTokenList(Statement &stat)
{
  Token token = readToken();

  while (! token.isEOF()) {
    // End of list if semicolon occurs outside of brackets.
    if (token.isDel(';')) break;
    stat.append(token);

    if (token.isDel('(')) {
      parseBracketList(')', stat);
    } else if (token.isDel('[')) {
      parseBracketList(']', stat);
    } else if (token.isDel('{')) {
      parseBracketList('}', stat);
    }

    token = readToken();
  }
}


Token MadParser::readToken()
{
  if (inputStack.empty()) {
    return Token("", 0, Token::IS_EOF, "End of input");
  } else {
    return inputStack.back()->readToken();
  }
}


Statement *MadParser::readStatement(TokenStream *is) const
{
  Statement *stat = 0;
  Token token = is->readToken();
  string name;

  try {
    if (token.isDel('{')) {
      // Compound statement.
      inputStack.back()->putBack(token);
      stat = new CompoundStatement(*inputStack.back());
    } else if (token.isKey("IF")) {
      // IF statement.
      inputStack.back()->putBack(token);
      stat = new IfStatement(*this, *inputStack.back());
    } else if (token.isKey("WHILE")) {
      // WHILE statement.
      inputStack.back()->putBack(token);
      stat = new WhileStatement(*this, *inputStack.back());
    } else if (token.isWord() || token.isString()) {
      // Simple statement or MACRO statement.
      stat = new SimpleStatement(token.getFile(), token.getLine());
      stat->append(token);
      token = is->readToken();

      if (! token.isEOF()) {
	if (token.isDel('(')) {
	  // Macro statement; statement already contains initial word.
	  stat->append(token);
	  parseBracketList(')', *stat);
	  token = is->readToken();
	
	  if (! token.isEOF() && token.isDel(':')) {
	    // Macro definition.
	    stat->append(token);
	    token = is->readToken();
	
	    if (! token.isEOF()) {
	      stat->append(token);
	      if (token.isKey("MACRO")) {
		token = is->readToken();
		
		if (! token.isEOF() && token.isDel('{')) {
		  stat->append(token);
		  parseBracketList('}', *stat);
		} else {
		  throw ParseError("MadParser::readStatement()",
				   "MACRO definition lacks \"{...}\".");
		}
	      } else {
		parseTokenList(*stat);
	      }
	    }
	  } else if (! token.isDel(';')) {
	    throw ParseError("MadParser::readStatement()",
			     "MACRO call is not terminated by ';'.");
	  }
	} else if (! token.isDel(';')) {
	  stat->append(token);
	  parseTokenList(*stat);
	}
      }
      stat->start();
    } else if (token.isDel(';')) {
      // Skip empty statement.
      stat = readStatement(is);
    } else if (token.isDel('?')) {
      // Give help.
      cerr << "\ntry typing \"HELP\" or \"SHOW\" for help.\n" << endl;
      stat = readStatement(is);
    } else if (! token.isEOF()) {
      stat = new SimpleStatement(token.getFile(), token.getLine());
      stat->append(token);
      parseTokenList(*stat);
      stat->start();
      throw ParseError("MadParser::readStatement()",
		       "Command should begin with a <name>.");
    }
  } catch (ParseError &ex) {
    cerr << "\n*** Parse error detected by function \""
	 << "MadParser::readStatement()" << "\"\n";
    stat->printWhere(cerr, true);
    cerr << "    ";
    stat->print(cerr);
    cerr << "    " << ex.what() << '\n' << endl;
    stat = readStatement(is);
  }

  return stat;
}


void MadParser::run() const
{
  stopFlag = false;

  while (Statement *stat = readStatement(&*inputStack.back())) {
    try {
      // The dispatch via Statement::execute() allows a special
      // treatment of structured statements.
      stat->execute(*this);
    } catch (ParseError &ex) {
      cerr << "\n*** Parse error detected by function \""
	   << ex.where() << "\"\n";
      stat->printWhere(cerr, true);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    " << ex.what() << '\n' << endl;
    } catch (MadException &ex) {
      cerr << "\n*** User error detected by function \""
	   << ex.where() << "\"\n";
      stat->printWhere(cerr, true);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    " << ex.what() << '\n' << endl;
    } catch (ClassicException &ex) {
      cerr << "\n*** User error detected by function \""
	   << ex.where() << "\"\n";
      stat->printWhere(cerr, false);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    " << ex.what() << '\n' << endl;
    } catch (bad_alloc &) {
      cerr << "\n*** Error:\n";
      stat->printWhere(cerr, false);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    Sorry, virtual memory exhausted.\n" << endl;
    } catch (exception &ex) {
      cerr << "\n*** Error:\n";
      stat->printWhere(cerr, false);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    Internal MAD error: " << ex.what()	<< '\n' << endl;
    } catch (...) {
      cerr << "\n*** Error:\n";
      stat->printWhere(cerr, false);
      cerr << "    ";
      stat->print(cerr);
      cerr << "    Unexpected exception caught.\n" << endl;
      abort();
    }

    delete stat;
    if (stopFlag) break;
  }
}


void MadParser::run(TokenStream *is) const
{
  inputStack.push_back(is);
  run();
  inputStack.pop_back();
}


void MadParser::stop() const
{
  stopFlag = true;
}
