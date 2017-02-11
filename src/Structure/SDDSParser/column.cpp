//
//  Copyright & License: See Copyright.readme in src directory
//

#include "column_def.hpp"

unsigned int SDDS::column::count_m = 0;

typedef std::string::const_iterator iterator_t;
template struct SDDS::parser::column_parser<iterator_t>;
