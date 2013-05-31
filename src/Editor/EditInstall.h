#ifndef OPAL_EditInstall_HH
#define OPAL_EditInstall_HH

// ------------------------------------------------------------------------
// $RCSfile: EditInstall.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditInstall
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"
#include "AbstractObjects/Element.h"
#include "MemoryManagement/Pointer.h"


// Class EditInstall
// ------------------------------------------------------------------------
/// The sequence editor INSTALL command.

class EditInstall: public Editor {

public:

    /// Exemplar constructor.
    EditInstall();

    virtual ~EditInstall();

    /// Make clone.
    virtual EditInstall *clone(const string &name);

    /// Execute the command.
    virtual void execute();

    /// Parse the command.
    //  Special format for this sub-command.
    virtual void parse(Statement &statement);

private:

    // Not implemented.
    EditInstall(const EditInstall &);
    void operator=(const EditInstall &);

    // Clone constructor.
    EditInstall(const string &name, EditInstall *parent);

    // A pointer to the new element to be installed.
    Pointer<Element> newElement;
};

#endif // OPAL_EditInstall_HH
