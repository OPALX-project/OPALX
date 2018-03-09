#ifndef OPAL_MPWriter_HH
#define OPAL_MPWriter_HH

// ------------------------------------------------------------------------
// $RCSfile: MPWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FieldReader
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPBase.h"


// Class MPWriter
// ------------------------------------------------------------------------
/// A DOOM writer for writing field errors.
//  Retrieves field errors from a beam line and uses an AlignHandler to
//  store them in the doom data base.

class MPWriter: public MPBase {

public:

    /// Constructor.
    //  Store [b]name[/b] in the DOOM environment.
    explicit MPWriter(const std::string &name);

    virtual ~MPWriter();

    /// Handle error field.
    virtual void fieldError(const std::string &name, int occur,
                            const BMultipoleField &designField,
                            BMultipoleField &errorField);

private:

    // Not implemented.
    MPWriter();
    MPWriter(const MPWriter &);
    void operator=(const MPWriter &);
};

#endif // OPAL_MPWriter_HH
