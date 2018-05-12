#ifndef OPAL_SequenceSampleCmd_HH
#define OPAL_SequenceSampleCmd_HH

#include "AbstractObjects/Action.h"

#include <string>

// Class SequenceSampleCmd
// ------------------------------------------------------------------------
/// The RUN command.

class SequenceSampleCmd: public Action {

public:

    /// Exemplar constructor.
    SequenceSampleCmd();

    virtual ~SequenceSampleCmd();

    /// Make clone.
    virtual SequenceSampleCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    SequenceSampleCmd(const SequenceSampleCmd &);
    void operator=(const SequenceSampleCmd &);

    // Clone constructor.
    SequenceSampleCmd(const std::string &name, SequenceSampleCmd *parent);

    void stashEnvironment();
    void popEnvironment();
};

#endif