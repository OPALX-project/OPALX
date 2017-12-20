//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef SKIPPER_HPP_
#define SKIPPER_HPP_

#include <boost/spirit/include/qi.hpp>

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    ///////////////////////////////////////////////////////////////////////////////
    //  The skipper grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct skipper : qi::grammar<Iterator>
    {
        skipper() : skipper::base_type(start)
        {
        	qi::eol_type eol;
        	qi::eoi_type eoi;
        	qi::char_type char_;
            ascii::space_type space;

            start =
                    space
				|   "!" >> *(char_ - eol) >> (eol|eoi)       // comments
                ;
        }

        qi::rule<Iterator> start;
    };

    template <typename Iterator>
    struct listskipper : qi::grammar<Iterator>
    {
        listskipper() : listskipper::base_type(start)
        {
            qi::char_type char_;
            ascii::space_type space;

            start =
                    space
                |   char_(',')
                ;
        }

        qi::rule<Iterator> start;
    };

}}

#endif


