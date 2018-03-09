#ifndef __PLAIN_INPUTFILE_PARSER_H__
#define __PLAIN_INPUTFILE_PARSER_H__

#include "Util/Types.h"
#include "Util/CommentAnnotatedInputFileParser.h"

/**
 *  \class PlainInputFileParser
 *  \brief Implements a parser for plain input files
 *
 *  This class extracts the optimization problem (objectives, constraints,
 *  design variables) from a plain input file, e.g.,
 *
 *  \verbatim
  d1:    DVAR, ELEMENT="", VARIABLE="SIGX";
  obj1:  OBJECTIVE, EXPR="(energy*energy)";
  objs:  OBJECTIVES = (obj1);
  dvars: DVARS = (d1);
  cons:  CONSTRAINTS = ();
  opt:   OPTIMIZE, OBJECTIVES=objs, DVARS=dvars, CONSTRAINTS=cons;
\endverbatim
 */
class PlainInputFileParser : public CommentAnnotatedInputFileParser {

public:
    PlainInputFileParser(std::string filename,
                         functionDictionary_t known_expr_funcs)
        : CommentAnnotatedInputFileParser(filename, "", known_expr_funcs)
    {}

    ~PlainInputFileParser()
    {}

};

#endif
