#ifndef __INPUTFILEPARSER_H__
#define __INPUTFILEPARSER_H__

#include <string>
#include <vector>
#include <map>

#include "Expression/Parser/function.hpp"

/// general data structures for storing objectives, constraints and design variables
typedef std::vector<std::string> StringList_t;

enum AttributeType_t {STRING, VALUE, EXPRESSION};
typedef std::pair<AttributeType_t, std::string> Attribute_t;
typedef std::vector<Attribute_t> Attributes_t;
typedef std::map<std::string, Attributes_t> AttributeMap_t;
typedef std::map<std::string, AttributeMap_t> OptInfo_t;

/**
 *  \class InputFileParser
 *  \brief Base class for all input file parsers
 *
 *  Every input file parser must implement doParse and getProblem.
 */
class InputFileParser {

public:
    InputFileParser(std::string filename,
                    functionDictionary_t known_expr_funcs)
        : filename_(filename)
        , known_expr_funcs_(known_expr_funcs)
        {}

    ~InputFileParser() {}

    /**
     *  Parse the input file and set objectives, constraints and design
     *  variables.
     */
    virtual void doParse() = 0;

    /**
     *  Access (for optimizer) to objectives constraints and design variables
     *  if implementation is using a different representation.
     *  @see Types.h for type definitions
     *  @param objectives container to store objectives
     *  @param constraints container to store constraints
     *  @param dvars container to store design variables
     */
    virtual void getProblem(Expressions::Named_t &objectives,
                            Expressions::Named_t &constraints,
                            DVarContainer_t &dvars) = 0;

    /// Getter for objectives
    OptInfo_t getObjectives()  { return objectives_;  }
    /// Getter for constraints
    OptInfo_t getConstraints() { return constraints_; }
    /// Getter for design variables
    OptInfo_t getDesignVars()  { return designVars_;  }

protected:
    /// filename of the input file
    std::string filename_;
    functionDictionary_t known_expr_funcs_;
    /// container for objectives
    OptInfo_t objectives_;
    /// container for constraints
    OptInfo_t constraints_;
    /// container for design variables
    OptInfo_t designVars_;
};

#endif
