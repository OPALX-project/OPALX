// ------------------------------------------------------------------------
// $RCSfile: OpalParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class OpalParser:
//   This is the default parser for OPAL statements.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:17:27 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Attributes/Attributes.h"
#include "OpalParser/CompoundStatement.h"
#include "OpalParser/IfStatement.h"
#include "OpalParser/WhileStatement.h"
#include "MemoryManagement/Pointer.h"
#include "Parser/SimpleStatement.h"
#include "Parser/Token.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <cassert>
#include <ctime>
#include <exception>
#include <iostream>
#include <new>

#include <Ippl.h>
using namespace std;

using namespace Expressions;
using std::cerr;
using std::endl;

extern Inform *gmsg;

// Class OpalParser
// ------------------------------------------------------------------------

std::vector<Pointer<TokenStream> > OpalParser::inputStack;


OpalParser::OpalParser(): stopFlag(false)
{}


OpalParser::~OpalParser()
{}


void OpalParser::parse(Statement &stat) const {
    if(stat.keyword("SHARED")) {
        // "SHARED ...": Shared object definition.
        parseDefine(stat);
    } else if(stat.keyword("CONSTANT") || stat.keyword("CONST") ||
              stat.keyword("BOOL") || stat.keyword("REAL") ||
              stat.keyword("STRING") || stat.keyword("VECTOR")) {
        // Keywords introducing variable definitions.
        parseAssign(stat);
    } else {
        std::string name = parseString(stat, "Identifier or keyword expected.");

        if(stat.delimiter('?')) {
            // "<class>?": give help for class.
            printHelp(name);
        } else if(stat.delimiter(':')) {
            // "<object>:<class>...": labeled command.
            parseDefine(stat);
        } else if(stat.delimiter('(')) {
            // "<macro>(...)...": macro definition or call.
            // We are positioned just after the '(' of the argument list.
            parseMacro(name, stat);
        } else if(stat.delimiter(',') || stat.delimiter(';') ||
                  stat.atEnd()) {
            // "<class>" or "<class>,<attributes>": Executable command.
            parseAction(stat);
        } else {
            // Assignment beginning with a name.
            parseAssign(stat);
        }
    }
}


void OpalParser::execute(Object *object, const std::string &name) const {
    // Trace execution.
    if(Options::mtrace && object->shouldTrace()) {
        double time = double(clock()) / double(CLOCKS_PER_SEC);
        *gmsg << "\nBegin execution: \"" << name
              << "\", CPU time = " << time << " seconds.\n" << endl;
    }

    // Force updating of all attributes which might have been changed.
    if(object->shouldUpdate()) {
        OpalData::getInstance()->update();
    }

    // Execute or check the command.
    object->execute();

    // Trace execution.
    if(Options::mtrace && object->shouldTrace()) {
        double time = double(clock()) / double(CLOCKS_PER_SEC);
        *gmsg << "\nEnd execution:   \"" << name
              << "\", CPU time = " << time << " seconds.\n" << endl;
    }
}


Object *OpalParser::find(const std::string &name) const {
    return OpalData::getInstance()->find(name);
}


void OpalParser::parseAction(Statement &stat) const {
    stat.start();
    std::string cmdName = parseString(stat, "Command name expected");

    if(cmdName == "STOP") {
        stopFlag = true;
    } else if(cmdName == "QUIT") {
        stopFlag = true;
    } else if(cmdName == "HELP"  &&  stat.delimiter(',')) {
        cmdName = parseString(stat, "Object name expected");
        printHelp(cmdName);
    } else if(Object *object = find(cmdName)) {
        Object *copy = 0;
        try {
            copy = object->clone("");
            copy->parse(stat);
            parseEnd(stat);
            execute(copy, cmdName);
            delete copy;
        } catch(...) {
            delete copy;
            throw;
        }
    } else {
        throw ParseError("OpalParser::parseAction()",
                         "Command \"" + cmdName + "\" is unknown.");
    }
}


void OpalParser::parseAssign(Statement &stat) const {
    stat.start();

    // Find various model objects.
    /*static*/
    Object *boolConstant   = OpalData::getInstance()->find("BOOL_CONSTANT");
    /*static*/
    Object *realConstant   = OpalData::getInstance()->find("REAL_CONSTANT");
    /*static*/
    Object *realVariable   = OpalData::getInstance()->find("REAL_VARIABLE");
    /*static*/
    Object *realVector     = OpalData::getInstance()->find("REAL_VECTOR");
    /*static*/
    Object *stringConstant = OpalData::getInstance()->find("STRING_CONSTANT");

    // Gobble up any prefix.
    int code = 0x00;
    while(true) {
        if(stat.keyword("CONSTANT") || stat.keyword("CONST")) {
            code |= 0x01;
        } else if(stat.keyword("BOOL")) {
            code |= 0x02;
        } else if(stat.keyword("REAL")) {
            code |= 0x04;
        } else if(stat.keyword("STRING")) {
            code |= 0x08;
        } else if(stat.keyword("VECTOR")) {
            code |= 0x10;
        } else {
            break;
        }
    }

    std::string objName = parseString(stat, "Object name expected.");

    // Test for attribute name.
    Object *object = 0;
    std::string attrName;

    if(stat.delimiter("->")) {
        // Assignment to object attribute.
        attrName = parseString(stat, "Attribute name expected.");

        if(code != 0) {
            throw ParseError("OpalParser::parseAssign()",
                             "Invalid type specification for this value.");
        } else if((object = OpalData::getInstance()->find(objName)) == 0) {
            throw ParseError("OpalParser::parseAssign()",
                             "The object \"" + objName + "\" is unknown.");
        }
    } else {
        // Assignment to variable-like object.
        if((object = OpalData::getInstance()->find(objName)) == 0) {
            Object *model = 0;
            switch(code) {
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

            if(model != 0) {
                object = model->clone(objName);
                OpalData::getInstance()->define(object);
            } else {
                throw ParseError("OpalParser::parseAssign()", "Invalid <type> field.");
            }
        } else if(object->isTreeMember(realConstant)) {
            throw ParseError("OpalParser::parseAssign()",
                             "You cannot redefine the constant \"" + objName + "\".");
        }

        attrName = "VALUE";
    }

    // Test for index; it is evaluated immediately.
    int index = 0;

    if(stat.delimiter('[')) {
        index = int(Round(parseRealConst(stat)));
        parseDelimiter(stat, ']');

        if(index <= 0) {
            throw ParseError("Expressions::parseReference()",
                             "Index must be positive.");
        }
    }

    if(object != 0) {
        if(Attribute *attr = object->findAttribute(attrName)) {
            if(stat.delimiter('=') || object->isTreeMember(realConstant)) {
                if(index > 0) {
                    attr->parseComponent(stat, true, index);
                } else {
                    attr->parse(stat, true);
                }
            } else if(stat.delimiter(":=")) {
                if(index > 0) {
                    attr->parseComponent(stat, false, index);
                } else {
                    attr->parse(stat, false);
                }
            }
        } else {
            throw ParseError("OpalParser::parseAssign()",
                             "Object \"" + objName + "\" has no attribute \"" +
                             attrName + "\".");
        }

        parseEnd(stat);
        OpalData::getInstance()->makeDirty(object);
    }
}


void OpalParser::parseDefine(Statement &stat) const {
    stat.start();
    bool isShared = stat.keyword("SHARED");
    std::string objName = parseString(stat, "Object name expected.");

    if(stat.delimiter(':')) {
        std::string clsName = parseString(stat, "Class name expected.");
        Object *classObject = find(clsName);

        if(classObject == 0) {
            throw ParseError("OpalParser::parseDefine()",
                             "The object \"" + clsName + "\" is unknown.");
        }

        Object *copy = 0;
        try {
            if(stat.delimiter('(')) {
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
            OpalData::getInstance()->define(copy);
        } catch(...) {
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


void OpalParser::parseEnd(Statement &stat) const {
    if(! stat.atEnd()  &&  ! stat.delimiter(';')) {
        throw ParseError("OpalParser::parseEnd()",
                         "Syntax error (maybe missing comma or semicolon ? )");
    }
}


void OpalParser::parseMacro(const std::string &macName, Statement &stat) const {
    // Record the position just after the '(' of the argument list.
    stat.mark();

    // Skip argument list.
    int par_level = 1;
    while(true) {
        if(stat.delimiter('(')) {
            ++par_level;
        } else if(stat.delimiter(')')) {
            if(--par_level == 0) break;
        } else {
            stat.getCurrent();
        }
    }

    if(stat.delimiter(':')) {
        // Macro definition.
        std::string className = parseString(stat, "Class name expected.");

        if(Object *macro = OpalData::getInstance()->find(className)) {
            // Backtrack to first argument.
            stat.restore();

            if(Object *copy =
                   macro->makeTemplate(macName, *inputStack.back(), stat)) {
                OpalData::getInstance()->define(copy);
            } else {
                throw ParseError("OpalParser::parseMacro()", "Command \"" +
                                 macName + "\" cannot be defined with arguments.");
            }
        } else {
            throw ParseError("OpalParser::parseMacro()",
                             "Object \"" + className + "\" is unknown.");
        }
    } else {
        // Macro call.
        if(Object *macro = OpalData::getInstance()->find(macName)) {
            // Backtrack to first argument.
            stat.restore();
            Object *instance = 0;
            try {
                instance = macro->makeInstance(macName, stat, this);
                execute(instance, macName);
            } catch(...) {
                delete instance;
                throw;
            }
        } else {
            throw ParseError("OpalParser::parseMacro()",
                             "Macro \"" + macName + "\" is unknown.");
        }
    }
}


void OpalParser::printHelp(const std::string &cmdName) const {
    Object *object = find(cmdName);

    if(object == 0) {
        *gmsg << "\nOpalParser::printHelp(): Unknown object \""
              << cmdName << "\".\n" << endl;
    } else {
        object->printHelp(cerr);
    }
}


void OpalParser::parseBracketList(char close, Statement &stat) {
    Token token = readToken();

    while(! token.isEOF()) {
        stat.append(token);

        if(token.isDel('(')) {
            parseBracketList(')', stat);
        } else if(token.isDel('[')) {
            parseBracketList(']', stat);
        } else if(token.isDel('{')) {
            parseBracketList('}', stat);
        } else if(token.isDel(close)) {
            return;
        }

        token = readToken();
    }
}


void OpalParser::parseTokenList(Statement &stat) {
    Token token = readToken();

    while(! token.isEOF()) {
        // End of list if semicolon occurs outside of brackets.
        if(token.isDel(';')) break;
        stat.append(token);

        if(token.isDel('(')) {
            parseBracketList(')', stat);
        } else if(token.isDel('[')) {
            parseBracketList(']', stat);
        } else if(token.isDel('{')) {
            parseBracketList('}', stat);
        }

        token = readToken();
    }
}


Token OpalParser::readToken() {
    if(inputStack.empty()) {
        return Token("", 0, Token::IS_EOF, "End of input");
    } else {
        return inputStack.back()->readToken();
    }
}


Statement *OpalParser::readStatement(TokenStream *is) const {
    Statement *stat = 0;
    Token token = is->readToken();
    std::string name;

    try {
        if(token.isDel('{')) {
            // Compound statement.
            inputStack.back()->putBack(token);
            stat = new CompoundStatement(*inputStack.back());
        } else if(token.isKey("IF")) {
            // IF statement.
            inputStack.back()->putBack(token);
            stat = new IfStatement(*this, *inputStack.back());
        } else if(token.isKey("WHILE")) {
            // WHILE statement.
            inputStack.back()->putBack(token);
            stat = new WhileStatement(*this, *inputStack.back());
        } else if(token.isWord() || token.isString()) {
            // Simple statement or MACRO statement.
            stat = new SimpleStatement(token.getFile(), token.getLine());
            stat->append(token);
            token = is->readToken();

            if(! token.isEOF()) {
                if(token.isDel('(')) {
                    // Macro statement; statement already contains initial word.
                    stat->append(token);
                    parseBracketList(')', *stat);
                    token = is->readToken();

                    if(! token.isEOF() && token.isDel(':')) {
                        // Macro definition.
                        stat->append(token);
                        token = is->readToken();

                        if(! token.isEOF()) {
                            stat->append(token);
                            if(token.isKey("MACRO")) {
                                token = is->readToken();

                                if(! token.isEOF() && token.isDel('{')) {
                                    stat->append(token);
                                    parseBracketList('}', *stat);
                                } else {
                                    throw ParseError("OpalParser::readStatement()",
                                                     "MACRO definition lacks \"{...}\".");
                                }
                            } else {
                                parseTokenList(*stat);
                            }
                        }
                    } else if(! token.isDel(';')) {
                        throw ParseError("OpalParser::readStatement()",
                                         "MACRO call is not terminated by ';'.");
                    }
                } else if(! token.isDel(';')) {
                    stat->append(token);
                    parseTokenList(*stat);
                }
            }
            stat->start();
        } else if(token.isDel(';')) {
            // Skip empty statement.
            stat = readStatement(is);
        } else if(token.isDel('?')) {
            // Give help.
            *gmsg << "\ntry typing \"HELP\" or \"SHOW\" for help.\n" << endl;
            stat = readStatement(is);
        } else if(! token.isEOF()) {
            stat = new SimpleStatement(token.getFile(), token.getLine());
            stat->append(token);
            parseTokenList(*stat);
            stat->start();
            throw ParseError("OpalParser::readStatement()",
                             "Command should begin with a <name>.");
        }
    } catch(ParseError &ex) {
        *gmsg << "\n*** Parse error detected by function \""
              << "OpalParser::readStatement()" << "\"\n";
        stat->printWhere(true);
        *gmsg << "    ";
        stat->print();
        *gmsg << "    " << ex.what() << '\n' << endl;
        stat = readStatement(is);
	exit(1);
    }

    return stat;
}


void OpalParser::run() const {
    stopFlag = false;

    while(Statement *stat = readStatement(&*inputStack.back())) {
        try {
            // The dispatch via Statement::execute() allows a special
            // treatment of structured statements.
            stat->execute(*this);
        } catch(ParseError &ex) {
            *gmsg << "\n*** Parse error detected by function \""
                  << ex.where() << "\"\n";
            stat->printWhere(true);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    " << ex.what() << '\n' << endl;
	    exit(1);
        } catch(OpalException &ex) {
            *gmsg << "\n*** User error detected by function \""
                  << ex.where() << "\"\n";
            stat->printWhere(true);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    " << ex.what() << '\n' << endl;
        } catch(ClassicException &ex) {
            *gmsg << "\n*** User error detected by function \""
                  << ex.where() << "\"\n";
            stat->printWhere(false);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    " << ex.what() << '\n' << endl;
        } catch(bad_alloc &) {
            *gmsg << "\n*** Error:\n";
            stat->printWhere(false);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    Sorry, virtual memory exhausted.\n" << endl;
        } catch(exception &ex) {
            *gmsg << "\n*** Error:\n";
            stat->printWhere(false);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    Internal OPAL error: " << ex.what() << '\n' << endl;
        } catch(...) {
            *gmsg << "\n*** Error:\n";
            stat->printWhere(false);
            *gmsg << "    ";
            stat->print();
            *gmsg << "    Unexpected exception caught.\n" << endl;
            abort();
        }

        delete stat;
        if(stopFlag) break;
    }
}


void OpalParser::run(TokenStream *is) const {
    inputStack.push_back(is);
    run();
    inputStack.pop_back();
}


void OpalParser::stop() const {
    stopFlag = true;
}
