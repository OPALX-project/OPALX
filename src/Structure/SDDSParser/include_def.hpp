//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef INCLUDE_DEF_HPP_
#define INCLUDE_DEF_HPP_

#include "include.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    include_parser<Iterator>::include_parser(error_handler<Iterator> & _error_handler):
        include_parser::base_type(start)
    {
        using qi::on_error;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        qi::_3_type _3;
        qi::_4_type _4;
        qi::lexeme_type lexeme;
        qi::char_type char_;
        qi::lit_type lit;
        qi::alpha_type alpha;
        qi::alnum_type alnum;
        qi::_pass_type _pass;

        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
        string %= quoted_string
                | lexeme[(alpha | '_') >> *(alnum | '_')];

        include_filename = lit("filename") > '=' > string;

        auto complainInclude = phx::bind(&include::complainUnsupported<include::INCLUDE>::apply);

        start =
                lit("&include")[_pass = complainInclude]
                > lit("&end");

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                "Error! Expecting ", _4, _3));
    }
}}
#endif /* INCLUDE_DEF_HPP_ */
