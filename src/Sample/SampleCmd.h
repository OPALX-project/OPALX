#ifndef OPAL_SampleCmd_HH
#define OPAL_SampleCmd_HH

#include "AbstractObjects/Action.h"

#include <string>

// Class SampleCmd
// ------------------------------------------------------------------------
/// The RUN command.

class SampleCmd: public Action {

public:

    /// Exemplar constructor.
    SampleCmd();

    virtual ~SampleCmd();

    /// Make clone.
    virtual SampleCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    SampleCmd(const SampleCmd &);
    void operator=(const SampleCmd &);

    // Clone constructor.
    SampleCmd(const std::string &name, SampleCmd *parent);

    void stashEnvironment();
    void popEnvironment();
};

#endif