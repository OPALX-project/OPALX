//
//  Copyright & License: See Copyright.readme in src directory
//

#include "parameter_def.hpp"

unsigned int SDDS::parameter::count_m = 0;

typedef std::string::const_iterator iterator_t;
template struct SDDS::parser::parameter_parser<iterator_t>;
