//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef FILE_DEF_HPP_
#define FILE_DEF_HPP_

#include "file.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    file_parser<Iterator>::file_parser(error_handler<Iterator> & _error_handler):
        file_parser::base_type(start),
        version_m(_error_handler),
        description_m(_error_handler),
        parameter_m(_error_handler),
        column_m(_error_handler),
        data_m(_error_handler),
        associate_m(_error_handler),
        array_m(_error_handler),
        include_m(_error_handler)
    {
        using qi::on_error;
        using qi::on_success;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        qi::_val_type _val;
        qi::skip_type skip;

        auto push_backParameter = phx::push_back(phx::at_c<2>(_val), _1);
        auto push_backColumn = phx::push_back(phx::at_c<3>(_val), _1);
        auto push_backAssociate = phx::push_back(phx::at_c<5>(_val), _1);
        auto push_backArray = phx::push_back(phx::at_c<6>(_val), _1);
        auto push_backInclude = phx::push_back(phx::at_c<7>(_val), _1);

        start =
                version_m[phx::at_c<0>(_val) = _1]
              > -description_m[phx::at_c<1>(_val) = _1]
              > *(parameter_m[push_backParameter]
                | column_m[push_backColumn]
                | associate_m[push_backAssociate]
                | array_m[push_backArray]
                | include_m[push_backInclude]
                )
              > data_m[phx::at_c<4>(_val) = _1]
              ;

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                "Error! Expecting ", _4, _3));

    }
}}
#endif /* FILE_DEF_HPP_ */
