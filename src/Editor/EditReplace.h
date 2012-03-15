#ifndef OPAL_EditReplace_HH
#define OPAL_EditReplace_HH

// ------------------------------------------------------------------------
// $RCSfile: EditReplace.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditReplace
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditReplace
// ------------------------------------------------------------------------
/// The sequence editor REPLACE command.

class EditReplace: public Editor {

public:

    /// Exemplar constructor.
    EditReplace();

    virtual ~EditReplace();

    /// Make clone.
    virtual EditReplace *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    EditReplace(const EditReplace &);
    void operator=(const EditReplace &);

    // Clone constructor.
    EditReplace(const string &name, EditReplace *parent);
};

#endif // OPAL_EditReplace_HH
