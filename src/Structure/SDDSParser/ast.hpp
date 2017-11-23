//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef AST_HPP_
#define AST_HPP_

#include <boost/variant.hpp>
#include <boost/spirit/include/qi.hpp>

#include <string>
#include <vector>

namespace SDDS {
    namespace ast {
        enum datatype { FLOAT
                      , DOUBLE
                      , SHORT
                      , LONG
                      , CHARACTER
                      , STRING };

        enum datamode { ASCII
                      , BINARY};

        enum endianess { BIGENDIAN
                       , LITTLEENDIAN};

        struct nil {};

        typedef boost::variant<float,
                               double,
                               short,
                               long,
                               char,
                               std::string> variant_t;

        typedef std::vector<variant_t> columnData_t;

        inline
        std::string getDataTypeString(datatype type) {
            switch(type) {
            case FLOAT:
                return "float";
            case DOUBLE:
                return "double";
            case SHORT:
                return "short";
            case LONG:
                return "long";
            case CHARACTER:
                return "char";
            case STRING:
                return "string";
            default:
                return "unknown";
            }
        }

    }

    namespace parser {
        namespace qi = boost::spirit::qi;

        template <typename Iterator, typename Skipper>
        struct string: qi::grammar<Iterator, std::string(), Skipper >
        {
            string();

            boost::spirit::qi::rule<Iterator, std::string(), Skipper> start;
        };

        template <typename Iterator, typename Skipper>
        struct qstring: qi::grammar<Iterator, std::string(), Skipper >
        {
            qstring();

            boost::spirit::qi::rule<Iterator, std::string(), Skipper> start;
        };
}
}

#endif // AST_HPP_