#ifndef OPAL_EditEnd_HH
#define OPAL_EditEnd_HH

// ------------------------------------------------------------------------
// $RCSfile: EditEnd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditEnd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditEnd
// ------------------------------------------------------------------------
/// The sequence editor ENDEDIT command.

class EditEnd: public Editor {

public:

    /// Exemplar constructor.
    EditEnd();

    virtual ~EditEnd();

    /// Make clone.
    virtual EditEnd *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    EditEnd(const EditEnd &);
    void operator=(const EditEnd &);

    // Clone constructor.
    EditEnd(const string &name, EditEnd *parent);
};

#endif // OPAL_EditEnd_HH
