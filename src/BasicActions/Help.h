#ifndef OPAL_Help_HH
#define OPAL_Help_HH

// ------------------------------------------------------------------------
// $RCSfile: Help.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Help
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Help
// ------------------------------------------------------------------------
/// The HELP commands.

class Help: public Action {

public:

    /// Exemplar constructor.
    Help();

    virtual ~Help();

    /// Make clone.
    virtual Help *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    Help(const Help &);
    void operator=(const Help &);

    // Clone constructor.
    Help(const std::string &name, Help *parent);
};

#endif // OPAL_Help_HH
