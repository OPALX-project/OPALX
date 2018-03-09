#ifndef OPAL_Quit_HH
#define OPAL_Quit_HH

// ------------------------------------------------------------------------
// $RCSfile: Quit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Quit
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Quit
// ------------------------------------------------------------------------
/// The QUIT command.

class Quit: public Action {

public:

    /// Exemplar constructor.
    Quit();

    virtual ~Quit();

    /// Make clone.
    virtual Quit *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Quit(const Quit &);
    void operator=(const Quit &);

    // Clone constructor.
    Quit(const std::string &name, Quit *parent);
};

#endif // OPAL_Quit_HH
