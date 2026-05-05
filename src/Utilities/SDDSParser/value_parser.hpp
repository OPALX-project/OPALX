//
// Simple value parser for SDDS data values
// Replaces boost::spirit::qi parsing for parameter and column values
//
#ifndef VALUE_PARSER_HPP_
#define VALUE_PARSER_HPP_

#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include "ast.hpp"

namespace SDDS {
    namespace parser {

        class ValueParser {
        public:
            ValueParser(const std::string& input, size_t start_pos = 0)
                : input_(input), pos_(start_pos) {}

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

            bool parseFloat(float& value) {
                skipWhitespace();
                size_t start = pos_;
                bool has_dot = false;
                bool has_exp = false;

                if (pos_ < input_.length() && (input_[pos_] == '+' || input_[pos_] == '-')) {
                    pos_++;
                }

                while (pos_ < input_.length()) {
                    char c = input_[pos_];
                    if (std::isdigit(static_cast<unsigned char>(c))) {
                        pos_++;
                    } else if (c == '.' && !has_dot) {
                        has_dot = true;
                        pos_++;
                    } else if ((c == 'e' || c == 'E') && !has_exp) {
                        has_exp = true;
                        pos_++;
                        if (pos_ < input_.length()
                            && (input_[pos_] == '+' || input_[pos_] == '-')) {
                            pos_++;
                        }
                    } else {
                        break;
                    }
                }

                if (pos_ > start) {
                    try {
                        value = std::stof(input_.substr(start, pos_ - start));
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
                return false;
            }

            bool parseDouble(double& value) {
                skipWhitespace();
                size_t start = pos_;
                bool has_dot = false;
                bool has_exp = false;

                if (pos_ < input_.length() && (input_[pos_] == '+' || input_[pos_] == '-')) {
                    pos_++;
                }

                while (pos_ < input_.length()) {
                    char c = input_[pos_];
                    if (std::isdigit(static_cast<unsigned char>(c))) {
                        pos_++;
                    } else if (c == '.' && !has_dot) {
                        has_dot = true;
                        pos_++;
                    } else if ((c == 'e' || c == 'E') && !has_exp) {
                        has_exp = true;
                        pos_++;
                        if (pos_ < input_.length()
                            && (input_[pos_] == '+' || input_[pos_] == '-')) {
                            pos_++;
                        }
                    } else {
                        break;
                    }
                }

                if (pos_ > start) {
                    try {
                        value = std::stod(input_.substr(start, pos_ - start));
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
                return false;
            }

            bool parseShort(short& value) {
                skipWhitespace();
                size_t start = pos_;

                if (pos_ < input_.length() && (input_[pos_] == '-' || input_[pos_] == '+')) {
                    pos_++;
                }

                while (pos_ < input_.length()
                       && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                    pos_++;
                }

                if (pos_ > start) {
                    try {
                        long long_val = std::stol(input_.substr(start, pos_ - start));
                        value         = static_cast<short>(long_val);
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
                return false;
            }

            bool parseLong(long& value) {
                skipWhitespace();
                size_t start = pos_;

                if (pos_ < input_.length() && (input_[pos_] == '-' || input_[pos_] == '+')) {
                    pos_++;
                }

                while (pos_ < input_.length()
                       && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                    pos_++;
                }

                if (pos_ > start) {
                    try {
                        value = std::stol(input_.substr(start, pos_ - start));
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
                return false;
            }

            bool parseChar(char& value) {
                skipWhitespace();
                if (pos_ < input_.length()) {
                    value = input_[pos_];
                    pos_++;
                    return true;
                }
                return false;
            }

            bool parseQuotedString(std::string& value) {
                skipWhitespace();
                if (pos_ >= input_.length() || input_[pos_] != '"') {
                    return false;
                }
                pos_++;  // skip opening quote
                value.clear();

                while (pos_ < input_.length() && input_[pos_] != '"') {
                    if (input_[pos_] == '\\' && pos_ + 1 < input_.length()) {
                        pos_++;
                        value += input_[pos_];
                    } else {
                        value += input_[pos_];
                    }
                    pos_++;
                }

                if (pos_ < input_.length() && input_[pos_] == '"') {
                    pos_++;  // skip closing quote
                    return true;
                }
                return false;
            }

            bool parseString(std::string& value) {
                skipWhitespace();
                if (pos_ < input_.length() && input_[pos_] == '"') {
                    return parseQuotedString(value);
                }

                // Parse unquoted string (identifier)
                size_t start = pos_;
                while (pos_ < input_.length()) {
                    char c = input_[pos_];
                    if (std::isalnum(static_cast<unsigned char>(c)) || c == '@' || c == '#'
                        || c == ':' || c == '+' || c == '-' || c == '%' || c == '.' || c == '_'
                        || c == '$' || c == '&' || c == '/') {
                        pos_++;
                    } else {
                        break;
                    }
                }

                if (pos_ > start) {
                    value = input_.substr(start, pos_ - start);
                    return true;
                }
                return false;
            }

            size_t getPosition() const { return pos_; }
            void setPosition(size_t pos) { pos_ = pos; }

        private:
            const std::string& input_;
            size_t pos_;
        };

    }  // namespace parser
}  // namespace SDDS

#endif /* VALUE_PARSER_HPP_ */
