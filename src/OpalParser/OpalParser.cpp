
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
#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Attributes/Attributes.h"
#include "OpalParser/CompoundStatement.h"
#include "OpalParser/IfStatement.h"
#include "OpalParser/Macro.h"
#include "OpalParser/WhileStatement.h"
#include "ValueDefinitions/RealConstant.h"

#include <cmath>
#include <ctime>
#include <exception>
#include <iostream>
#include <memory>
#include <new>
#include "OpalParser/SimpleStatement.h"
#include "OpalParser/Token.h"
#include "OpalParser/TokenStream.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"

#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"

using namespace Expressions;

extern Inform* gmsg;

// Class OpalParser
// ------------------------------------------------------------------------

std::vector<std::shared_ptr<TokenStream> > OpalParser::inputStack;

OpalParser::OpalParser() : stopFlag(false) {}

OpalParser::~OpalParser() {}

void OpalParser::parse(Statement& stat) const {
    if (stat.keyword("SHARED")) {
        // "SHARED ...": Shared object definition.
        parseDefine(stat);
    } else if (
            stat.keyword("CONSTANT") || stat.keyword("CONST") || stat.keyword("BOOL")
            || stat.keyword("REAL") || stat.keyword("STRING") || stat.keyword("VECTOR")) {
        // Keywords introducing variable definitions.
        parseAssign(stat);
    } else {
        std::string name = parseString(stat, "Identifier or keyword expected.");

        if (stat.delimiter('?')) {
            // "<class>?": give help for class.
            printHelp(name);
        } else if (stat.delimiter(':')) {
            // "<object>:<class>...": labeled command.
            parseDefine(stat);
        } else if (stat.delimiter('(')) {
            // "<macro>(...)...": macro definition or call.
            // We are positioned just after the '(' of the argument list.
            parseMacro(name, stat);
        } else if (stat.delimiter(',') || stat.delimiter(';') || stat.atEnd()) {
            // "<class>" or "<class>,<attributes>": Executable command.
            parseAction(stat);
        } else {
            // Assignment beginning with a name.
            stat.mark();
            stat.start();
            std::string objName = parseString(stat, "Object name expected.");
            if (OpalData::getInstance()->find(objName) == 0) {
                Token tok = stat.getCurrent();
                stat.restore();

                std::string name              = tok.getLex();
                std::string hint              = getHint(name);
                unsigned int position         = stat.position();
                std::string positionIndicator = std::string(position, ' ') + "^\n";
                std::ostringstream statStr;
                stat.print(statStr);
                if (!hint.empty()) {
                    throw ParseError(
                            "OpalParser::parse()",
                            statStr.str() + positionIndicator
                                    + "Syntax error, either the keyword REAL is missing or\n"
                                    + hint);
                } else {
                    throw ParseError(
                            "OpalParser::parse()",
                            statStr.str() + positionIndicator
                                    + "Syntax error, the keyword REAL is missing\n");
                }
            }
            parseAssign(stat);
        }
    }
}

void OpalParser::execute(Object* object, const std::string& name) const {
    // Trace execution.
    if (Options::mtrace && object->shouldTrace()) {
        double time = double(clock()) / double(CLOCKS_PER_SEC);
        *gmsg << "\nBegin execution: \"" << name << "\", CPU time = " << time << " seconds.\n"
              << endl;
    }

    // Force updating of all attributes which might have been changed.
    if (object->shouldUpdate()) {
        OpalData::getInstance()->update();
    }

    // Execute or check the command.
    object->execute();

    // Trace execution.
    if (Options::mtrace && object->shouldTrace()) {
        double time = double(clock()) / double(CLOCKS_PER_SEC);
        *gmsg << "\nEnd execution:   \"" << name << "\", CPU time = " << time << " seconds.\n"
              << endl;
    }
}

Object* OpalParser::makeFrom(const std::string& name, const std::string& newName) const {
    OpalData* opal = OpalData::getInstance();
    if (const Object* object = opal->find(name)) {
        if (dynamic_cast<const Macro*>(object) != nullptr) {
            throw ParseError(
                    "OpalParser::makeFrom()",
                    "\"" + name + "\" is a macro; call it with its arguments.");
        }
        return opal->copyObject(*object, newName);
    }
    return opal->makeObject(name, newName);
}

void OpalParser::parseAction(Statement& stat) const {
    stat.start();
    std::string cmdName = parseString(stat, "Command name expected");

    if (cmdName == "STOP") {
        stopFlag = true;
    } else if (cmdName == "QUIT") {
        stopFlag = true;
    } else if (cmdName == "HELP" && stat.delimiter(',')) {
        cmdName = parseString(stat, "Object name expected");
        printHelp(cmdName);
    } else if (std::unique_ptr<Object> copy{makeFrom(cmdName, "")}) {
        copy->parse(stat);
        parseEnd(stat);
        execute(copy.get(), cmdName);
    } else {
        std::string hint = getHint(cmdName, "command");
        if (!hint.empty()) {
            throw ParseError("OpalParser::parseAction()", "Syntax error, " + hint);
        }

        throw ParseError("OpalParser::parseAction()", "Command \"" + cmdName + "\" is unknown.");
    }
}

void OpalParser::parseAssign(Statement& stat) const {
    stat.start();

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

    std::string objName = parseString(stat, "Object name expected.");
    // Test for attribute name.
    Object* object = 0;
    std::string attrName;

    if (stat.delimiter("->")) {
        // Assignment to object attribute.
        attrName = parseString(stat, "Attribute name expected.");

        if (code != 0) {
            throw ParseError(
                    "OpalParser::parseAssign()", "Invalid type specification for this value.");
        } else if ((object = OpalData::getInstance()->find(objName)) == 0) {
            if (OpalData::getInstance()->isKeyword(objName)) {
                throw ParseError(
                        "OpalParser::parseAssign()",
                        "\"" + objName + "\" is a keyword, not a defined object.");
            }
            throw ParseError(
                    "OpalParser::parseAssign()", "The object \"" + objName + "\" is unknown.");
        }
    } else {
        // Assignment to variable-like object.
        if ((object = OpalData::getInstance()->find(objName)) == 0) {
            const char* keyword = nullptr;
            switch (code) {
                case 0x01:  // CONSTANT
                case 0x05:  // CONSTANT REAL
                    keyword = "REAL_CONSTANT";
                    break;
                case 0x02:  // BOOL
                case 0x03:  // BOOL CONSTANT
                    keyword = "BOOL_CONSTANT";
                    break;
                case 0x00:  // empty <type>.
                case 0x04:  // REAL
                    keyword = "REAL_VARIABLE";
                    break;
                case 0x10:  // VECTOR
                case 0x11:  // CONSTANT VECTOR
                case 0x14:  // REAL VECTOR
                case 0x15:  // CONSTANT REAL VECTOR
                    keyword = "REAL_VECTOR";
                    break;
                case 0x08:  // STRING
                case 0x09:  // STRING CONSTANT
                    keyword = "STRING_CONSTANT";
                    break;
                default:
                    break;
            }

            if (keyword != nullptr) {
                object = OpalData::getInstance()->makeObject(keyword, objName);
                OpalData::getInstance()->define(object);
            } else {
                throw ParseError("OpalParser::parseAssign()", "Invalid <type> field.");
            }
        } else if (dynamic_cast<RealConstant*>(object) != nullptr) {
            throw ParseError(
                    "OpalParser::parseAssign()",
                    "You cannot redefine the constant \"" + objName + "\".");
        }

        attrName = "VALUE";
    }

    // Test for index; it is evaluated immediately.
    int index = 0;

    if (stat.delimiter('[')) {
        index = int(std::round(parseRealConst(stat)));
        parseDelimiter(stat, ']');

        if (index <= 0) {
            throw ParseError("Expressions::parseReference()", "Index must be positive.");
        }
    }

    if (object != 0) {
        if (Attribute* attr = object->findAttribute(attrName)) {
            if (stat.delimiter('=') || dynamic_cast<RealConstant*>(object) != nullptr) {
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
            throw ParseError(
                    "OpalParser::parseAssign()",
                    "Object \"" + objName + "\" has no attribute \"" + attrName + "\".");
        }

        parseEnd(stat);
        OpalData::getInstance()->makeDirty(object);
    }
}

void OpalParser::parseDefine(Statement& stat) const {
    stat.start();
    stat.keyword("SHARED");  // accepted; it has no effect
    std::string objName = parseString(stat, "Object name expected.");

    if (stat.delimiter(':')) {
        std::string clsName = parseString(stat, "Class name expected.");

        if (clsName == "SURFACEPHYSICS") {
            throw ParseError(
                    "OpalParser::parseDefine()",
                    "The object \"" + clsName + "\" is changed to \"PARTICLEMATTERINTERACTION\".");
        }

        std::unique_ptr<Object> copy;
        if (stat.delimiter('(')) {
            // Macro-like objects are always classes, instances never.
            // There is no further check required.
            Object* macro = OpalData::getInstance()->find(clsName);
            if (macro == nullptr) {
                throw ParseError(
                        "OpalParser::parseDefine()", "The object \"" + clsName + "\" is unknown.");
            }
            copy.reset(macro->makeInstance(objName, stat, this));
        } else {
            copy.reset(makeFrom(clsName, objName));
            if (copy == nullptr) {
                throw ParseError(
                        "OpalParser::parseDefine()", "The object \"" + clsName + "\" is unknown.");
            }
            copy->parse(stat);
        }

        parseEnd(stat);
        execute(copy.get(), clsName);
        OpalData::getInstance()->define(copy.get());
        copy.release();
    } else {
        // "SHARED NAME, ...;": parsed for errors, has no effect.
        std::unique_ptr<Object> copy(makeFrom(objName, objName));
        if (copy == nullptr) {
            throw ParseError(
                    "OpalParser::parseDefine()", "The object \"" + objName + "\" is unknown.");
        }
        copy->parse(stat);
    }
}

void OpalParser::parseEnd(Statement& stat) const {
    if (!stat.atEnd() && !stat.delimiter(';')) {
        unsigned int position         = stat.position();
        std::string positionIndicator = std::string(position + 1, ' ') + "^\n";
        std::ostringstream statStr;
        stat.print(statStr);

        throw ParseError(
                "OpalParser::parseEnd()",
                statStr.str() + positionIndicator
                        + "Syntax error (maybe missing comma or semicolon ? )");
    }
}

void OpalParser::parseMacro(const std::string& macName, Statement& stat) const {
    // Record the position just after the '(' of the argument list.
    stat.mark();

    // Skip argument list.
    int par_level = 1;
    while (true) {
        if (stat.atEnd()) {
            throw ParseError("OpalParser::parseMacro()", "MACRO argument list is not closed.");
        }
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
        std::string className = parseString(stat, "Class name expected.");

        // A defined macro or line, or a new object of keyword MACRO or LINE.
        Object* macro = OpalData::getInstance()->find(className);
        std::unique_ptr<Object> keywordObject;
        if (macro == nullptr) {
            keywordObject.reset(OpalData::getInstance()->makeObject(className, ""));
            macro = keywordObject.get();
        }
        if (macro != nullptr) {
            // Backtrack to first argument.
            stat.restore();

            if (Object* copy = macro->makeTemplate(macName, *inputStack.back(), stat)) {
                OpalData::getInstance()->define(copy);
            } else {
                throw ParseError(
                        "OpalParser::parseMacro()",
                        "Command \"" + macName + "\" cannot be defined with arguments.");
            }
        } else {
            throw ParseError(
                    "OpalParser::parseMacro()", "Object \"" + className + "\" is unknown.");
        }
    } else {
        // Macro call.
        if (Object* macro = OpalData::getInstance()->find(macName)) {
            // Backtrack to first argument.
            stat.restore();
            std::unique_ptr<Object> instance(macro->makeInstance(macName, stat, this));
            execute(instance.get(), macName);
        } else {
            throw ParseError("OpalParser::parseMacro()", "Macro \"" + macName + "\" is unknown.");
        }
    }
}

void OpalParser::printHelp(const std::string& cmdName) const {
    Object* object = OpalData::getInstance()->find(cmdName);
    std::unique_ptr<Object> keywordObject;
    if (object == nullptr) {
        keywordObject.reset(OpalData::getInstance()->makeObject(cmdName, cmdName));
        object = keywordObject.get();
    }

    if (object == 0) {
        *gmsg << "\nOpalParser::printHelp(): Unknown object \"" << cmdName << "\".\n" << endl;
    } else {
        object->printHelp(std::cerr);
    }
}

void OpalParser::parseBracketList(char close, Statement& stat) {
    Token token = readToken();

    while (!token.isEOF()) {
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

void OpalParser::parseTokenList(Statement& stat) {
    Token token = readToken();

    while (!token.isEOF()) {
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

Token OpalParser::readToken() {
    if (inputStack.empty()) {
        return Token("", 0, Token::IS_EOF, "End of input");
    } else {
        return inputStack.back()->readToken();
    }
}

Statement* OpalParser::readStatement(TokenStream* is) const {
    Statement* stat = 0;
    Token token     = is->readToken();

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

            if (!token.isEOF()) {
                if (token.isDel('(')) {
                    // Macro statement; statement already contains initial word.
                    stat->append(token);
                    parseBracketList(')', *stat);
                    token = is->readToken();

                    if (!token.isEOF() && token.isDel(':')) {
                        // Macro definition.
                        stat->append(token);
                        token = is->readToken();

                        if (!token.isEOF()) {
                            stat->append(token);
                            if (token.isKey("MACRO")) {
                                token = is->readToken();

                                if (!token.isEOF() && token.isDel('{')) {
                                    stat->append(token);
                                    parseBracketList('}', *stat);
                                } else {
                                    throw ParseError(
                                            "OpalParser::readStatement()",
                                            "MACRO definition lacks \"{...}\".");
                                }
                            } else {
                                parseTokenList(*stat);
                            }
                        }
                    } else if (!token.isDel(';')) {
                        throw ParseError(
                                "OpalParser::readStatement()",
                                "MACRO call is not terminated by ';'.");
                    }
                } else if (!token.isDel(';')) {
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
            *gmsg << "\ntry typing \"HELP\" for help.\n" << endl;
            stat = readStatement(is);
        } else if (!token.isEOF()) {
            stat = new SimpleStatement(token.getFile(), token.getLine());
            stat->append(token);
            parseTokenList(*stat);
            stat->start();
            throw ParseError("OpalParser::readStatement()", "Command should begin with a <name>.");
        }
    } catch (ParseError& ex) {
        *ippl::Error << "\n*** Parse error detected by function \""
                     << "OpalParser::readStatement()"
                     << "\"\n";
        /// \todo check this /stat->printWhere(*IpplInfo::Error, true);

        std::string what           = ex.what();
        std::string::size_type pos = 0;
        while ((pos = what.find("\n", pos)) != std::string::npos) {
            what.replace(pos, 1, "\n    ");
            pos += 5;  // length of "\n    "
        }

        *ippl::Error << "     " << *stat << "    a" << what << '\n' << endl;

        stat = readStatement(is);
        exit(1);
    }

    return stat;
}

void OpalParser::run() const {
    stopFlag = false;
    while (std::unique_ptr<Statement> stat{readStatement(&*inputStack.back())}) {
        try {
            // The dispatch via Statement::execute() allows a special
            // treatment of structured statements.
            stat->execute(*this);
        } catch (ParseError& ex) {
            Inform errorMsg("Error", std::cerr);
            errorMsg << "\n*** Parse error detected by function \"" << ex.where() << "\"\n";
            stat->printWhere(errorMsg, true);
            std::string what = ex.what();
            size_t pos       = what.find_first_of('\n');
            do {
                errorMsg << "    " << what.substr(0, pos) << endl;
                what = what.substr(pos + 1, std::string::npos);
                pos  = what.find_first_of('\n');
            } while (pos != std::string::npos);
            errorMsg << "    " << what << endl;

            MPI_Abort(MPI_COMM_WORLD, -100);
        }

        if (stopFlag) break;
    }
}

void OpalParser::run(TokenStream* is) const {
    inputStack.push_back(std::shared_ptr<TokenStream>(is));
    try {
        run();
    } catch (...) {
        inputStack.pop_back();
        throw;
    }
    inputStack.pop_back();
}

void OpalParser::stop() const { stopFlag = true; }

std::string OpalParser::getHint(const std::string& name, const std::string& type) {
    auto owner = AttributeHandler::getOwner(name);
    if (owner.empty()) {
        return std::string();
    }

    std::string hint = "the " + type + " '" + name + "' could belong to\n";
    {
        std::string elements;
        auto its = owner.equal_range(AttributeHandler::ELEMENT);
        if (its.first != its.second) {
            elements = (its.first)->second;
            bool any = (its.first)->second == "Any";
            for (auto it = std::next(its.first); it != its.second && !any; ++it) {
                elements += ", " + it->second;
                any = it->second == "Any";
            }
            if (any) {
                hint += std::string("  - any element\n");
            } else {
                hint += std::string("  - the element")
                        + (std::distance(its.first, its.second) > 1 ? "s " : " ") + elements + "\n";
            }
        }
    }
    {
        std::string commands;
        auto its = owner.equal_range(AttributeHandler::COMMAND);
        if (its.first != its.second) {
            commands = (its.first)->second;
            for (auto it = std::next(its.first); it != its.second; ++it) {
                commands += ", " + it->second;
            }
            hint += std::string("  - the command")
                    + (std::distance(its.first, its.second) > 1 ? "s " : " ") + commands + "\n";
        }
    }
    {
        std::string sub_commands;
        auto its = owner.equal_range(AttributeHandler::SUB_COMMAND);
        if (its.first != its.second) {
            sub_commands = (its.first)->second;
            for (auto it = std::next(its.first); it != its.second; ++it) {
                sub_commands += ", " + it->second;
            }
            hint += std::string("  - the sub-command")
                    + (std::distance(its.first, its.second) > 1 ? "s " : " ") + sub_commands + "\n";
        }
    }
    {
        std::string statements;
        auto its = owner.equal_range(AttributeHandler::STATEMENT);
        if (its.first != its.second) {
            statements = (its.first)->second;
            for (auto it = std::next(its.first); it != its.second; ++it) {
                statements += ", " + it->second;
            }
            hint += std::string("  - the statement")
                    + (std::distance(its.first, its.second) > 1 ? "s " : " ") + statements + "\n";
        }
    }

    hint += "but it's not present!";
    return hint;
}
