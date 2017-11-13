#include "Structure/SDDSParser.h"

#include <boost/algorithm/string.hpp>

SDDS::SDDSParser::SDDSParser():
    sddsFileName_m("")
{ }

SDDS::SDDSParser::SDDSParser(const std::string &input):
    sddsFileName_m(input)
{ }

void SDDS::SDDSParser::setInput(const std::string &input) {
    sddsFileName_m = input;
}

SDDS::file SDDS::SDDSParser::run() {
    typedef std::string::const_iterator iterator_t;
    typedef SDDS::parser::file_parser<iterator_t> file_parser_t;
    typedef SDDS::parser::skipper<iterator_t> skipper_t;
    typedef SDDS::error_handler<iterator_t> error_handler_t;

    skipper_t skipper;
    std::string contents = readFile();
    iterator_t contentsIter = contents.begin();
    iterator_t contentsEnd = contents.end();
    error_handler_t error_handler(contentsIter, contentsEnd);
    file_parser_t parser(error_handler);

    bool success = phrase_parse(contentsIter, contentsEnd, parser, skipper, sddsData_m);
    {
        SDDS::parameterList::iterator piter = sddsData_m.sddsParameters_m.begin();
        SDDS::parameterList::iterator pend = sddsData_m.sddsParameters_m.end();
        for (; piter != pend && success; ++ piter) {
            success = piter->parse(contentsIter, contentsEnd, skipper);
        }
        while (success && contentsIter != contentsEnd) {
            SDDS::columnList::iterator citer = sddsData_m.sddsColumns_m.begin();
            SDDS::columnList::iterator cend = sddsData_m.sddsColumns_m.end();
            for (; citer != cend && success; ++ citer) {
                success = citer->parse(contentsIter, contentsEnd, skipper);
            }
        }
    }

    if (!success || contentsIter != contentsEnd)
        {
            throw OpalException("StatisticalErrors::parseSDDSFile",
                                "could not parse SDDS file");
        }

    for (const SDDS::parameter &param: sddsData_m.sddsParameters_m) {
        std::string name = *param.name_m;
        fixCaseSensitivity(name);
        paramNameToID_m.insert(std::make_pair(name,
                                              param.order_m));
    }

    for (const SDDS::column &col: sddsData_m.sddsColumns_m) {
        std::string name = *col.name_m;
        fixCaseSensitivity(name);
        columnNameToID_m.insert(std::make_pair(name,
                                               col.order_m));
    }

    std::cout << sddsData_m.sddsParameters_m.size() << std::endl;
    return sddsData_m;
}

std::string SDDS::SDDSParser::readFile() {
    std::ifstream in(sddsFileName_m.c_str());

    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);

        in.read(&contents[0], contents.size());

        in.close();

        return contents;
    }

    throw OpalException("StatisticalErrors::readSDDSFile",
                        "could not open file '" + sddsFileName_m + "'");

    return std::string("");
}


SDDS::ast::columnData_t SDDS::SDDSParser::getColumnData(const std::string &columnName) {
    for (const SDDS::column &col: sddsData_m.sddsColumns_m) {
        if (*col.name_m == columnName) {
            return col.values_m;
        }
    }
    throw OpalException("StatisticalErrors::getColumnData",
                        "could not find column '" + columnName + "'");
}

//XXX use either all upper, or all lower case chars
void SDDS::SDDSParser::fixCaseSensitivity(std::string &for_string) {

    boost::to_lower(for_string);
}