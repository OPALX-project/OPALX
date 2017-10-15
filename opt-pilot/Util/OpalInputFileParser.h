#ifndef __OPAL_INPUT_FILE_PARSER_H__
#define __OPAL_INPUT_FILE_PARSER_H__

#include <string>
#include <sstream>
#include <set>

#include "Util/Types.h"
#include "Util/CommentAnnotatedInputFileParser.h"

/**
 *  \class OpalInputFileParser
 *  \brief Implements a parser for OPAL input files
 *
 *  This class extracts the optimization problem (objectives, constraints,
 *  design variables) from an Opal input file. All statements are prefixed
 *  with "//", e.g.,
 *
 *  \verbatim
  d1: DVAR, ELEMENT="", VARIABLE="SIGX";
  obj1: OBJECTIVE, EXPR="(energy*energy)";
  objs: OBJECTIVES = (obj1);
  dvars: DVARS = (d1);
  constrs: CONSTRAINTS = ();
  opt: OPTIMIZE, OBJECTIVES=objs, DVARS=dvars, CONSTRAINTS=constrs;
\endverbatim
 */
class OpalInputFileParser : public CommentAnnotatedInputFileParser {

public:
    OpalInputFileParser(std::string filename,
                        functionDictionary_t known_expr_funcs)
        : CommentAnnotatedInputFileParser(filename, "", known_expr_funcs)
    {}

    ~OpalInputFileParser()
    {}

};

#endif
