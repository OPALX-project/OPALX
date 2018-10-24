//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/file.hpp"
#include "Util/SDDSParser/skipper.hpp"
#include "Util/SDDSParser/array.hpp"
#include "Util/SDDSParser/associate.hpp"
#include "Util/SDDSParser/column.hpp"
#include "Util/SDDSParser/data.hpp"
#include "Util/SDDSParser/description.hpp"
#include "Util/SDDSParser/include.hpp"
#include "Util/SDDSParser/parameter.hpp"
//#include "Util/SDDSParser/variant.h"

#include "Util/SDDSParser/SDDSParserException.h"


#include "boost/optional/optional_io.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#ifndef __SDDSPARSER_H__
#define __SDDSPARSER_H__

namespace SDDS {

    class SDDSParser {
    private:
        std::string readFile();
        static void fixCaseSensitivity(std::string &for_string);
        static std::string fixCaseSensitivity(const std::string &for_string) {
            std::string retval(for_string);
            fixCaseSensitivity(retval);
            return retval;
        }
        std::string sddsFileName_m;

        /// mapping from parameter name to offset in params_m
        std::map<std::string, int> paramNameToID_m;
        /// mapping from column name to ID in columns_m
        std::map<std::string, int> columnNameToID_m;

        SDDS::file sddsData_m;

    public:
        SDDSParser();
        SDDSParser(const std::string &input);
        void setInput(const std::string &input);
        file run();

        file getData();
        ast::columnData_t getColumnData(const std::string &columnName);

        ast::datatype getColumnType(const std::string &col_name) {
            int index = getColumnIndex(col_name);
            return *sddsData_m.sddsColumns_m[index].type_m;
        }

        /**
         *  Converts the string value of a parameter at timestep t to a value of
         *  type T.
         *
         *  @param t timestep (beginning at 1, -1 means last)
         *  @param param parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getValue(int t, std::string column_name, T& nval) {

            fixCaseSensitivity(column_name);

            int col_idx = getColumnIndex(column_name);

            // round timestep to last if not in range
            size_t row_idx = 0;
            size_t num_rows = sddsData_m.sddsColumns_m[col_idx].values_m.size();
            if(t <= 0 || static_cast<size_t>(t) > num_rows)
                row_idx = num_rows - 1;
            else
                row_idx = static_cast<size_t>(t) - 1;

            ast::variant_t val = sddsData_m.sddsColumns_m[col_idx].values_m[row_idx];
            nval = boost::get<T>(val);
        }


        /**
         *  Converts the string value of a parameter at a position spos to a value
         *  of type T.
         *
         *  @param spos interpolate value at spos
         *  @param param parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getInterpolatedValue(double spos, std::string col_name, T& nval) {

            T value_before = 0;
            T value_after  = 0;
            double value_before_spos = 0;
            double value_after_spos  = 0;


            size_t col_idx_spos = getColumnIndex("s");
            ast::columnData_t &spos_values = sddsData_m.sddsColumns_m[col_idx_spos].values_m;
            int index = getColumnIndex(col_name);
            ast::columnData_t &col_values = sddsData_m.sddsColumns_m[index].values_m;

            size_t this_row = 0;
            size_t num_rows = spos_values.size();
            for(this_row = 0; this_row < num_rows; this_row++) {
                value_after_spos = boost::get<double>(spos_values[this_row]);

                if(spos < value_after_spos) {

                    size_t prev_row = 0;
                    if(this_row > 0) prev_row = this_row - 1;

                    try {
                        switch (*sddsData_m.sddsColumns_m[index].type_m) {
                        case ast::FLOAT:
                            value_before = boost::get<float>(col_values[prev_row]);
                            value_after  = boost::get<float>(col_values[this_row]);
                            break;
                        case ast::DOUBLE:
                            value_before = boost::get<double>(col_values[prev_row]);
                            value_after  = boost::get<double>(col_values[this_row]);
                            break;
                        case ast::SHORT:
                            value_before = boost::get<short>(col_values[prev_row]);
                            value_after  = boost::get<short>(col_values[this_row]);
                            break;
                        case ast::LONG:
                            value_before = boost::get<long>(col_values[prev_row]);
                            value_after  = boost::get<long>(col_values[this_row]);
                            break;
                        default:
                            throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                                      "can't convert value to double");
                        }
                    } catch (...) {
                        throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                                  "can't convert value");
                    }

                    value_before_spos = boost::get<double>(spos_values[prev_row]);
                    value_after_spos  = boost::get<double>(spos_values[this_row]);

                    break;
                }
            }

            if(this_row == num_rows)
                throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                          "all values < specified spos");

            // simple linear interpolation
            if(spos - value_before_spos < 1e-8)
                nval = value_before;
            else
                nval = value_before + (spos - value_before_spos)
                    * (value_after - value_before)
                    / (value_after_spos - value_before_spos);
        }

        template <typename T>
        void getParameterValue(std::string parameter_name, T& nval) {
            fixCaseSensitivity(parameter_name);

            if (paramNameToID_m.count(parameter_name) > 0) {
                size_t id = paramNameToID_m[parameter_name];
                auto value = sddsData_m.sddsParameters_m[id].value_m;
                nval = boost::get<T>(value);
            } else {
                throw SDDSParserException("SDDSParser::getParameterValue",
                                        "unknown parameter name: '" + parameter_name + "'!");
            }
        }

    private:

        int getColumnIndex(std::string col_name) const;
    };

    inline
    file SDDSParser::getData() {
        return sddsData_m;
    }
}

#endif