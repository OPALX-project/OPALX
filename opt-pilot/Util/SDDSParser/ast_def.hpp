//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef AST_DEF_HPP_
#define AST_DEF_HPP_

#include "ast.hpp"

namespace SDDS { namespace parser
{

    template <typename Iterator, typename Skipper>
    string<Iterator, Skipper>::string():
        string::base_type(start)
    {
        qi::char_type char_;
        qi::eol_type eol;
        qi::lexeme_type lexeme;

        start %= lexeme[+(char_-eol)];
    }


    template <typename Iterator, typename Skipper>
    qstring<Iterator, Skipper>::qstring():
        qstring::base_type(start)
    {
        qi::char_type char_;
        qi::eol_type eol;
        qi::lexeme_type lexeme;

        start %= lexeme['"' >> +(char_ - (eol|'"')) >> '"']
               | lexeme[+(char_ - (eol|' '))];
    }
}}

#endif // AST_DEF_HPP_
