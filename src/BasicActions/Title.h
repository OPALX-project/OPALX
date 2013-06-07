#ifndef OPAL_Title_HH
#define OPAL_Title_HH

// ------------------------------------------------------------------------
// $RCSfile: Title.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Title
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Title
// ------------------------------------------------------------------------
/// The TITLE command.

class Title: public Action {

public:

    /// Exemplar constructor.
    Title();

    virtual ~Title();

    /// Make clone.
    virtual Title *clone(const string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement &);

private:

    // Not implemented.
    Title(const Title &);
    void operator=(const Title &);

    // Clone constructor.
    Title(const string &name, Title *parent);
};

#endif // OPAL_Title_HH

