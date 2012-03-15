#ifndef OPAL_EditFlatten_HH
#define OPAL_EditFlatten_HH

// ------------------------------------------------------------------------
// $RCSfile: EditFlatten.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditFlatten
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditFlatten
// ------------------------------------------------------------------------
/// The sequence editor FLATTEN command.

class EditFlatten: public Editor {

public:

    /// Exemplar constructor.
    EditFlatten();

    virtual ~EditFlatten();

    /// Make clone.
    virtual EditFlatten *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    EditFlatten(const EditFlatten &);
    void operator=(const EditFlatten &);

    // Clone constructor.
    EditFlatten(const string &name, EditFlatten *parent);
};

#endif // OPAL_EditFlatten_HH
