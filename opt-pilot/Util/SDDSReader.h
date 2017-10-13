#ifndef __SDDSREADER_H__
#define __SDDSREADER_H__

#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <locale>

#include "Util/OptPilotException.h"

/**
 *  \class SDDSReader
 *  \brief Implements a parser and value extractor for SDDS files
 *
 *  Statistical Opal simulation results are stored in SDDS files. This class
 *  grants the optimizer access to results of the forward solve to evaluate
 *  expressions (objectives, constraints).
 */
class SDDSReader {

public:

    SDDSReader(std::string filename);
    ~SDDSReader();

    void parseFile();

    /**
     *  Converts the string value of a parameter at timestep t to a value of
     *  type T.
     *
     *  @param t timestep (beginning at 1, -1 means last)
     *  @param param parameter name
     *  @param nval store result of type T in nval
     */
    template <typename T>
    void getValue(int t, std::string param_name, T& nval) {

        fixCaseSensitivity(param_name);

        // round timestep to last if not in range
        size_t col_idx = 0;
        if(t <= 0 || static_cast<size_t>(t) > columns_m.size())
            col_idx = columns_m.size() - 1;
        else
            col_idx = static_cast<size_t>(t) - 1;

        if(columnNameToID_m.count(param_name) > 0) {
            int index = columnNameToID_m[param_name];
            std::string val = (columns_m[col_idx])[index];
            nval = stringToNum<T>(val);
        } else {
            throw OptPilotException("SDDSReader::getValue",
                                    "unknown column name: '" + param_name + "'!");
        }
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
    void getInterpolatedValue(double spos, std::string param_name, T& nval) {

        fixCaseSensitivity(param_name);

        T value_before = 0;
        T value_after  = 0;
        double value_before_spos = 0;
        double value_after_spos  = 0;

        size_t index_spos = 1;

        if(columnNameToID_m.count(param_name) > 0) {
            int index = columnNameToID_m[param_name];

            size_t this_col = 0;
            for(this_col = 0; this_col < columns_m.size(); this_col++) {

                value_after_spos = stringToNum<double>(
                                        (columns_m[this_col])[index_spos]);

                if(spos < value_after_spos) {

                    size_t prev_col = 0;
                    if(this_col > 0) prev_col = this_col - 1;

                    value_before = stringToNum<T>(
                                        (columns_m[prev_col])[index]);
                    value_after  = stringToNum<T>(
                                        (columns_m[this_col])[index]);

                    value_before_spos = stringToNum<double>(
                                            (columns_m[prev_col])[index_spos]);
                    value_after_spos  = stringToNum<double>(
                                            (columns_m[this_col])[index_spos]);

                    break;
                }
            }

            if(this_col == columns_m.size())
                throw OptPilotException("SDDSReader::getInterpolatedValue",
                                        "all values < specified spos");

        } else {
            throw OptPilotException("SDDSReader::getInterpolatedValue",
                                    "unknown column name: '" + param_name + "'!");
        }

        // simple linear interpolation
        if(spos - value_before_spos < 1e-8)
            nval = value_before;
        else
            nval = value_before + (spos - value_before_spos)
                   * (value_after - value_before)
                   / (value_after_spos - value_before_spos);
    }

private:

    /// SDDS filename
    std::string filename_m;
    /// number of parameters found in SDDS header
    int numParams_m;
    /// number of columns in the SDDS file
    int numColumns_m;
    /// mapping from parameter name to offset in params_m
    std::map<std::string, int> paramNameToID_m;
    /// mapping from column name to ID in columns_m
    std::map<std::string, int> columnNameToID_m;

    /// vector holding all vectors of columns
    std::vector< std::vector<std::string> > columns_m;
    /// vector of all parameters
    std::vector< std::string > params_m;

    void fixCaseSensitivity(std::string &for_string);
    void addParameter(std::string param_name);
    void addColumn(std::string column_name);


    /**
     *  Split a 'delimiter' separated string into list
     *  @param str string to split
     *  @param delimiter separating items in string
     *  @return vector containing all string items
     */
    std::vector<std::string> split(std::string str, std::string delimiter);

    template <typename T>
    T stringToNum( const std::string& s ) {
        std::istringstream i(s);
        T x;
        if (!(i >> x))
            return 0;
        return x;
    }

};

#endif
