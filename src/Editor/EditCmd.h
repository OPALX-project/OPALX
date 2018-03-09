#ifndef OPAL_EditCmd_HH
#define OPAL_EditCmd_HH

// ------------------------------------------------------------------------
// $RCSfile: EditCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class EditCmd
// ------------------------------------------------------------------------
/// The sequence editor SEQEDIT command.

class EditCmd: public Action {

public:

    /// Exemplar constructor.
    EditCmd();

    virtual ~EditCmd();

    /// Make clone.
    virtual EditCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    EditCmd(const EditCmd &);
    void operator=(const EditCmd &);

    // Clone constructor.
    EditCmd(const std::string &name, EditCmd *parent);
};

#endif // OPAL_EditCmd_HH
