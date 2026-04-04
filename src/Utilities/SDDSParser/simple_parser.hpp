//
// Simple recursive descent parser for SDDS format
// Replaces boost::spirit::qi parser
//
// Copyright (c) 2024
// All rights reserved
//
#ifndef SIMPLE_PARSER_HPP_
#define SIMPLE_PARSER_HPP_

#include "ast.hpp"
#include "version.hpp"
#include "description.hpp"
#include "file.hpp"
#include "parameter.hpp"
#include "column.hpp"
#include "data.hpp"
#include "associate.hpp"
#include "array.hpp"
#include "include.hpp"

#include <string>
#include <sstream>
#include <cctype>
#include <stdexcept>

namespace SDDS {
    namespace parser {

        class SimpleParser {
        public:
            SimpleParser(const std::string& input) 
                : input_(input), pos_(0) {}

            file parse() {
                file result;
                skipWhitespace();
                
                // Parse version
                result.sddsVersion_m = parseVersion();
                skipWhitespace();
                
                // Parse optional description
                if (match("&description")) {
                    result.sddsDescription_m = parseDescription();
                    skipWhitespace();
                }
                
                // Parse parameters, columns, associates, arrays, includes
                while (pos_ < input_.length()) {
                    skipWhitespace();
                    if (match("&parameter")) {
                        result.sddsParameters_m.push_back(parseParameter());
                    } else if (match("&column")) {
                        result.sddsColumns_m.push_back(parseColumn());
                    } else if (match("&associate")) {
                        result.sddsAssociates_m.push_back(parseAssociate());
                    } else if (match("&array")) {
                        result.sddsArrays_m.push_back(parseArray());
                    } else if (match("&include")) {
                        result.sddsIncludes_m.push_back(parseInclude());
                    } else if (match("&data")) {
                        result.sddsData_m = parseData();
                        break;
                    } else {
                        break;
                    }
                    skipWhitespace();
                }
                
                return result;
            }

        private:
            std::string input_;
            size_t pos_;

            void skipWhitespace() {
                while (pos_ < input_.length()) {
                    if (std::isspace(static_cast<unsigned char>(input_[pos_]))) {
                        pos_++;
                    } else if (input_[pos_] == '!') {
                        // Skip comments
                        while (pos_ < input_.length() && input_[pos_] != '\n') {
                            pos_++;
                        }
                    } else {
                        break;
                    }
                }
            }

            bool match(const std::string& str) {
                skipWhitespace();
                if (input_.compare(pos_, str.length(), str) == 0) {
                    pos_ += str.length();
                    return true;
                }
                return false;
            }

            bool match(char c) {
                skipWhitespace();
                if (pos_ < input_.length() && input_[pos_] == c) {
                    pos_++;
                    return true;
                }
                return false;
            }

            void expect(const std::string& str) {
                if (!match(str)) {
                    throw std::runtime_error("Expected: " + str);
                }
            }

            void expect(char c) {
                if (!match(c)) {
                    throw std::runtime_error(std::string("Expected: ") + c);
                }
            }

            std::string parseQuotedString() {
                if (pos_ >= input_.length() || input_[pos_] != '"') {
                    return "";
                }
                pos_++; // skip opening quote
                std::string result;
                while (pos_ < input_.length() && input_[pos_] != '"') {
                    if (input_[pos_] == '\\' && pos_ + 1 < input_.length()) {
                        pos_++;
                        result += input_[pos_];
                    } else {
                        result += input_[pos_];
                    }
                    pos_++;
                }
                if (pos_ < input_.length() && input_[pos_] == '"') {
                    pos_++; // skip closing quote
                }
                return result;
            }

            std::string parseIdentifier() {
                skipWhitespace();
                std::string result;
                while (pos_ < input_.length()) {
                    char c = input_[pos_];
                    if (std::isalnum(static_cast<unsigned char>(c)) || 
                        c == '@' || c == '#' || c == ':' || c == '+' || 
                        c == '-' || c == '%' || c == '.' || c == '_' || 
                        c == '$' || c == '&' || c == '/') {
                        result += c;
                        pos_++;
                    } else {
                        break;
                    }
                }
                return result;
            }

            std::string parseString() {
                skipWhitespace();
                if (pos_ < input_.length() && input_[pos_] == '"') {
                    return parseQuotedString();
                }
                return parseIdentifier();
            }

            version parseVersion() {
                version v;
                expect("SDDS");
                skipWhitespace();
                std::string numStr;
                while (pos_ < input_.length() && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                    numStr += input_[pos_];
                    pos_++;
                }
                if (numStr.empty()) {
                    throw std::runtime_error("Expected version number");
                }
                v.layoutVersion_m = static_cast<short>(std::stoi(numStr));
                return v;
            }

            description parseDescription() {
                description desc;
                skipWhitespace();
                
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    if (match("text")) {
                        expect('=');
                        desc.text_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("contents")) {
                        expect('=');
                        desc.content_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else {
                        break;
                    }
                }
                expect("&end");
                return desc;
            }

            parameter parseParameter() {
                parameter param;
                skipWhitespace();
                
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    if (match("name")) {
                        expect('=');
                        param.name_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("type")) {
                        expect('=');
                        std::string typeStr = parseString();
                        if (typeStr == "float") param.type_m = ast::FLOAT;
                        else if (typeStr == "double") param.type_m = ast::DOUBLE;
                        else if (typeStr == "short") param.type_m = ast::SHORT;
                        else if (typeStr == "long") param.type_m = ast::LONG;
                        else if (typeStr == "character") param.type_m = ast::CHARACTER;
                        else if (typeStr == "string") param.type_m = ast::STRING;
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("units")) {
                        expect('=');
                        param.units_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("description")) {
                        expect('=');
                        param.description_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else {
                        // Skip unsupported fields
                        parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    }
                }
                expect("&end");
                return param;
            }

            column parseColumn() {
                column col;
                skipWhitespace();
                
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    if (match("name")) {
                        expect('=');
                        col.name_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("type")) {
                        expect('=');
                        std::string typeStr = parseString();
                        if (typeStr == "float") col.type_m = ast::FLOAT;
                        else if (typeStr == "double") col.type_m = ast::DOUBLE;
                        else if (typeStr == "short") col.type_m = ast::SHORT;
                        else if (typeStr == "long") col.type_m = ast::LONG;
                        else if (typeStr == "character") col.type_m = ast::CHARACTER;
                        else if (typeStr == "string") col.type_m = ast::STRING;
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("units")) {
                        expect('=');
                        col.units_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else if (match("description")) {
                        expect('=');
                        col.description_m = parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    } else {
                        // Skip unsupported fields
                        parseString();
                        skipWhitespace();
                        if (match(',')) skipWhitespace();
                    }
                }
                expect("&end");
                return col;
            }

            associate parseAssociate() {
                associate assoc;
                skipWhitespace();
                // Basic implementation - can be extended
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    parseString(); // Skip content for now
                    skipWhitespace();
                    if (match(',')) skipWhitespace();
                }
                expect("&end");
                return assoc;
            }

            array parseArray() {
                array arr;
                skipWhitespace();
                // Basic implementation - can be extended
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    parseString(); // Skip content for now
                    skipWhitespace();
                    if (match(',')) skipWhitespace();
                }
                expect("&end");
                return arr;
            }

            include parseInclude() {
                include inc;
                skipWhitespace();
                // Basic implementation - can be extended
                while (pos_ < input_.length() && !match("&end")) {
                    skipWhitespace();
                    parseString(); // Skip content for now
                    skipWhitespace();
                    if (match(',')) skipWhitespace();
                }
                expect("&end");
                return inc;
            }

            data parseData() {
                data d;
                skipWhitespace();
                // Parse data mode
                if (match("mode")) {
                    expect('=');
                    std::string modeStr = parseString();
                    if (modeStr == "ascii") d.mode_m = ast::ASCII;
                    else if (modeStr == "binary") d.mode_m = ast::BINARY;
                    skipWhitespace();
                    if (match(',')) skipWhitespace();
                }
                if (match("no_row_counts")) {
                    expect('=');
                    std::string val = parseString();
                    // Parse boolean value
                    skipWhitespace();
                    if (match(',')) skipWhitespace();
                }
                expect("&end");
                skipWhitespace();
                
                // Parse data rows (simplified - actual parsing would be more complex)
                // This is a placeholder - the actual data parsing is done separately
                // in the original code via parameter::parse() and column::parse()
                
                return d;
            }
        };

    } // namespace parser
} // namespace SDDS

#endif /* SIMPLE_PARSER_HPP_ */

