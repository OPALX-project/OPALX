#ifndef OPAL_OptimizeCmd_HH
#define OPAL_OptimizeCmd_HH

#include "AbstractObjects/Action.h"

#include <string>

// Class OptimizeCmd
// ------------------------------------------------------------------------
/// The RUN command.

class OptimizeCmd: public Action {

public:

    /// Exemplar constructor.
    OptimizeCmd();

    virtual ~OptimizeCmd();

    /// Make clone.
    virtual OptimizeCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    OptimizeCmd(const OptimizeCmd &);
    void operator=(const OptimizeCmd &);

    // Clone constructor.
    OptimizeCmd(const std::string &name, OptimizeCmd *parent);

    void stashEnvironment();
    void popEnvironment();
};

#endif