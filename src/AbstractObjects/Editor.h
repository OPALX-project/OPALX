#ifndef OPAL_Editor_HH
#define OPAL_Editor_HH

// ------------------------------------------------------------------------
// $RCSfile: Editor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Editor
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"


// Class Editor
// ------------------------------------------------------------------------
/// The base class for all OPAL sequence editor commands.
//  It implements the common behaviour of editor commands, it can also
//  be used via dynamic casting to determine whether an object represents
//  an editor command.

class Editor: public Object {

public:

    virtual ~Editor();

    /// Return the object category as a string.
    //  Return the string "EDITOR".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always true for editor commands.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always false for editor commands.
    virtual bool shouldUpdate() const;

protected:

    /// Constructor for exemplars.
    Editor(int size, const char *name, const char *help);

    /// Constructor for cloning.
    Editor(const std::string &name, Editor *parent);

private:

    // Not implemented.
    Editor();
    Editor(const Editor &);
    void operator=(const Editor &);
};

#endif // OPAL_Editor_HH
