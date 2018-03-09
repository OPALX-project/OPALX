//
//  Copyright & License: See Copyright.readme in src directory
//

#include "ast_def.hpp"
#include "skipper.hpp"

typedef std::string::const_iterator iterator_t;
typedef SDDS::parser::skipper<iterator_t> skipper_t;
template struct SDDS::parser::string<iterator_t, skipper_t>;
template struct SDDS::parser::qstring<iterator_t, skipper_t>;
