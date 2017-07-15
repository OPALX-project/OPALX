#ifndef __COMMENT_ANNOTATED_INPUTFILE_PARSER_H__
#define __COMMENT_ANNOTATED_INPUTFILE_PARSER_H__

#include <string>
#include <sstream>
#include <set>

#include "Types.h"
#include "InputFileParser.h"

/**
 *  \class CommentAnnotatedInputFileParser
 *  \brief Implements a parser for CommentAnnotated input files
 *
 *  This class extracts the optimization problem (objectives, constraints,
 *  design variables) from an comment annotated input file. All statements are
 *  prefixed with a comment symbol, e.g. "\\":
 *
 *  \verbatim
  //d1: DVAR, ELEMENT="", VARIABLE="SIGX";
  //obj1: OBJECTIVE, EXPR="(energy*energy)";
  //objs: OBJECTIVES = (obj1);
  //dvars: DVARS = (d1);
  //constrs: CONSTRAINTS = ();
  //opt: OPTIMIZE, OBJECTIVES=objs, DVARS=dvars, CONSTRAINTS=constrs;
\endverbatim
 */
class CommentAnnotatedInputFileParser : public InputFileParser {

public:
    CommentAnnotatedInputFileParser(std::string filename,
                                    std::string comment_symbol,
                                    functionDictionary_t known_expr_funcs);
    virtual ~CommentAnnotatedInputFileParser();

    /**
     *  Extracts optimizer information from CommentAnnotated input file.
     */
    void doParse();

    /**
     *  Access (for optimizer) to objectives constraints and design variables.
     *  @see Types.h for type definitions
     *  @param objectives container to store objectives
     *  @param constraints container to store constraints
     *  @param dvars container to store design variables
     */
    void getProblem(Expressions::Named_t &objectives,
                    Expressions::Named_t &constraints,
                    DVarContainer_t &dvars);

private:

    /// local container for objectives
    Expressions::Named_t nobjectives_;
    /// local container for constraints
    Expressions::Named_t nconstraints_;
    /// local container for design variables
    DVarContainer_t ndvars_;

    /// name of the simulation/optimization
    std::string optname_;
    /// string representation of list of objectives
    std::string listObjectives_;
    /// string representation of list of constraints
    std::string listConstraints_;
    /// string representation of list of design variables
    std::string listDVars_;

    /// set of names of all objectives
    std::set<std::string> objectives_;
    /// set of names of all constraints
    std::set<std::string> constraints_;
    /// set of names of all design variables
    std::set<std::string> dvars_;

    std::string comment_symbol_;

    /**
     *  Get a specific attribute in string
     *  @param str string to be search
     *  @param attribute attribute to search
     *  @return string containing trimmed attribute (if any)
     */
    std::string getAttribute(std::string str, std::string attribute);

    /**
     *  Get the list of all attributes
     *  @param str to get attributes from
     *  @return set of all attributes found in str
     */
    std::set<std::string> getListAttribute(std::string str);
};

#endif
