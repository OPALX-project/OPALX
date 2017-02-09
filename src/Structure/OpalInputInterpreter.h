#ifndef OPAL_INPUTINTERPRETER_H
#define OPAL_INPUTINTERPRETER_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

// #include "Algorithms/Tracker.h"
// #include "Algorithms/ParallelTTracker.h"

#include <boost/variant.hpp>
#include <iostream>
#include <list>

namespace ast {
    enum splitType { TEXT
                   , GAUSS
                   , CALL
                   , TGAUSS
                   };

    typedef boost::variant<int,
                           double,
                           std::string> variant_t;
}

struct OpalInputInterpreter {
    typedef std::list<ast::variant_t>::iterator iterator_t;

    OpalInputInterpreter(const std::string & fileName);
    OpalInputInterpreter();

    void parse(std::istream & in, const iterator_t &position);
    void parseTGauss(const iterator_t & begin, const iterator_t & end);
    void parseGauss(const iterator_t & begin, const iterator_t & end);
    void parseCall(const iterator_t & begin, const iterator_t & end);
    bool processCalls();
    bool getCall(iterator_t & it);
    bool ignoreCode(iterator_t & it) const;

    std::string processAST() const;
    std::string processASTReference() const;

    void replaceString(const std::string & search,
                       std::string replace);

    std::list<ast::variant_t> sourceParts_m;
};

#endif // OPAL_STATISTICALERRORS_H