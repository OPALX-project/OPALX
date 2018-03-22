//
//  Copyright & License: See Copyright.readme in src directory
//

#include "data.hpp"
#include "data_def.hpp"

typedef std::string::const_iterator iterator_t;
template struct SDDS::parser::data_parser<iterator_t>;
