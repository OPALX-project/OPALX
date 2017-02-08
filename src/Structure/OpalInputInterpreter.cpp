//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/OpalInputInterpreter.h"
// #include "AbstractObjects/OpalData.h"
// #include "OpalConfigure/Configure.h"
// #include "OpalParser/OpalParser.h"
// #include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Expressions/Expressions.h"
// #include "ValueDefinitions/StringConstant.h"

#include <boost/regex.hpp>
// #include <boost/algorithm/string/trim.hpp>

#include <string>
#include <iostream>
#include <fstream>

OpalInputInterpreter::OpalInputInterpreter(const std::string & fileName) {
    sourceParts_m.push_back(ast::CALL);
    sourceParts_m.push_back(fileName);

    while (processCalls()) { }
}

OpalInputInterpreter::OpalInputInterpreter() { }

bool OpalInputInterpreter::processCalls() {
    bool processedCalls = false;
    for (auto it = sourceParts_m.begin(); it != sourceParts_m.end(); ++ it) {
        if (getCall(it)) {
            std::string fileName = boost::get<std::string>(*it);
            std::ifstream ifh(fileName.c_str());
            it = sourceParts_m.erase(it); // remove fileName;
            processedCalls = true;

            parse(ifh, it);

            ifh.close();
            break;
        }
    }
    return processedCalls;
}

void OpalInputInterpreter::parse(std::istream & in,
                                 const OpalInputInterpreter::iterator_t &position) {
    std::string source("");
    std::string str;
    const std::string commentFormat("");
    const boost::regex cppCommentExpr("//.*");
    const boost::regex cCommentExpr("/\\*.*?\\*/");

    while (std::getline(in, str)) {
        str = boost::regex_replace(str, cppCommentExpr, commentFormat);
        source += str + '\n';

        if (!in.good() || in.eof())
            break;
    }

    iterator_t it = sourceParts_m.insert(position, ast::TEXT);
    sourceParts_m.insert(position, boost::regex_replace(source, cCommentExpr, commentFormat));

    parseTGauss(it, position);
    parseGauss(it, position);
    parseCall(it, position);
}

void OpalInputInterpreter::parseTGauss(const OpalInputInterpreter::iterator_t & begin,
                                       const OpalInputInterpreter::iterator_t & end) {
    const boost::regex tgaussExpr("TGAUSS\\s*\\(\\s*([+-]?\\s*\\d*\\.?\\d*E?[+-]?\\d*)\\s*\\)", boost::regex::icase);
    boost::smatch what;
    boost::match_flag_type flags = boost::match_default;

    for (iterator_t it = begin; it != end; ++ it) {
        if (ignoreCode(it)) continue;

        std::string text = boost::get<std::string>(*it);
        std::string::const_iterator start = text.begin();
        std::string::const_iterator end = text.end();

        while (boost::regex_search(start, end, what, tgaussExpr, flags)) {
            std::string part(start, what[0].first);
            std::string attribute_str(what[1].first, what[1].second);
            double attribute = std::stod(attribute_str);

            sourceParts_m.insert(it, part);
            sourceParts_m.insert(it, ast::TGAUSS);
            sourceParts_m.insert(it, attribute);
            sourceParts_m.insert(it, ast::TEXT);

            start = what[0].second;
        }

        std::string rest(start, end);
        *it = rest;
    }
}

void OpalInputInterpreter::parseGauss(const OpalInputInterpreter::iterator_t & begin,
                                      const OpalInputInterpreter::iterator_t & end) {
    const boost::regex gaussExpr("GAUSS\\s*\\(\\s*\\)", boost::regex::icase);
    boost::smatch what;
    boost::match_flag_type flags = boost::match_default;

    for (iterator_t it = begin; it != end; ++ it) {
        if (ignoreCode(it)) continue;

        std::string text = boost::get<std::string>(*it);
        std::string::const_iterator start = text.begin();
        std::string::const_iterator end = text.end();

        while (boost::regex_search(start, end, what, gaussExpr, flags)) {
            std::string part(start, what[0].first);

            sourceParts_m.insert(it, part);
            sourceParts_m.insert(it, ast::GAUSS);
            sourceParts_m.insert(it, ast::TEXT);

            start = what[0].second;
        }

        std::string rest(start, end);
        *it = rest;
    }
}

void OpalInputInterpreter::parseCall(const OpalInputInterpreter::iterator_t & begin,
                                     const OpalInputInterpreter::iterator_t & end) {
    const boost::regex callExpr("CALL\\s*,\\s*FILE\\s*=\\s*\"([a-zA-Z0-9_\\-/\\.]*)\"\\s*;", boost::regex::icase);
    boost::smatch what;
    boost::match_flag_type flags = boost::match_default;

    for (iterator_t it = begin; it != end; ++ it) {
        if (ignoreCode(it)) continue;

        std::string text = boost::get<std::string>(*it);
        std::string::const_iterator start = text.begin();
        std::string::const_iterator end = text.end();

        while (boost::regex_search(start, end, what, callExpr, flags)) {
            std::string part(start, what[0].first);
            std::string filename(what[1].first, what[1].second);

            sourceParts_m.insert(it, part);
            sourceParts_m.insert(it, ast::CALL);
            sourceParts_m.insert(it, filename);
            sourceParts_m.insert(it, ast::TEXT);

            start = what[0].second;
        }

        std::string rest(start, end);
        *it = rest;
    }

}

bool OpalInputInterpreter::getCall(OpalInputInterpreter::iterator_t & it) {
    switch (boost::get<int>(*it)) {
    case ast::TEXT:
    case ast::TGAUSS:
        ++ it;
        break;
    case ast::GAUSS:
        break;
    case ast::CALL:
        it = sourceParts_m.erase(it); // remove ast::CALL
        return true;
    default:
        throw OpalException("OpalInputInterpreter::getCall(iterator_t &)",
                            "ast code not working properly");
    }

    return false;
}

bool OpalInputInterpreter::ignoreCode(OpalInputInterpreter::iterator_t & it) const {
    switch (boost::get<int>(*it)) {
    case ast::TEXT:
        ++ it;
        break;
    case ast::TGAUSS:
    case ast::CALL:
        ++ it;
        return true;
    case ast::GAUSS:
        return true;
    default:
        throw OpalException("OpalInputInterpreter::ingoreCode()",
                            "ast code not working properly");
    }

    return false;
}

std::string OpalInputInterpreter::processAST() const {
    std::string source("");

    for (auto it = sourceParts_m.begin(); it != sourceParts_m.end(); ++ it) {
        switch (boost::get<int>(*it)) {
        case ast::TEXT:
            ++ it;
            source += boost::get<std::string>(*it);
            break;
        case ast::TGAUSS:
            ++ it;
            source += ("(" + std::to_string(Expressions::Tgauss(boost::get<double>(*it))) + ")");
            break;
        case ast::GAUSS:
            source += ("(" + std::to_string(Expressions::gauss()) + ")");
            break;
        default:
            throw OpalException("OpalInputInterpreter::processAST()",
                                "ast code not working properly");
        }
    }

    return source;
}

std::string OpalInputInterpreter::processASTReference() const {
    std::string source("");

    for (auto it = sourceParts_m.begin(); it != sourceParts_m.end(); ++ it) {
        switch (boost::get<int>(*it)) {
        case ast::TEXT:
            ++ it;
            source += boost::get<std::string>(*it);
            break;
        case ast::TGAUSS:
            ++ it;
            source += "0.0";
            break;
        case ast::GAUSS:
            source += "0.0";
            break;
        default:
            throw OpalException("OpalInputInterpreter::processAST()",
                                "ast code not working properly");
        }
    }

    return source;
}

void OpalInputInterpreter::replaceString(const std::string &search,
                                         std::string replace) {
    const boost::regex expr(search, boost::regex::icase);
    boost::smatch what;
    boost::match_flag_type flags = boost::match_default;

    for (iterator_t it = sourceParts_m.begin(); it != sourceParts_m.end(); ++ it) {
        if (ignoreCode(it)) continue;

        std::string text = boost::get<std::string>(*it);
        std::string::const_iterator start = text.begin();
        std::string::const_iterator end = text.end();

        while (boost::regex_search(start, end, what, expr, flags)) {
            std::string part(start, what[0].first);

            for (unsigned int i = 1; i < what.size(); ++ i) {
                const std::string whatstr(what[i].first, what[i].second);
                const std::string matchReplace = "$1" + whatstr + "$2";
                const std::string matchString = "^(.*?)\\$\\{" + std::to_string(i) + "\\}(.*?)$";
                const boost::regex matchRe(matchString);

                replace = boost::regex_replace(replace, matchRe, matchReplace);
            }

            sourceParts_m.insert(it, part);
            sourceParts_m.insert(it, ast::TEXT);
            sourceParts_m.insert(it, replace);
            sourceParts_m.insert(it, ast::TEXT);

            start = what[0].second;
        }

        std::string rest(start, end);
        *it = rest;
    }
}
