#ifndef OPAL_MPReader_HH
#define OPAL_MPReader_HH

// ------------------------------------------------------------------------
// $RCSfile: MPReader.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPReader
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPBase.h"


// Class MPReader
// ------------------------------------------------------------------------
/// A DOOM reader for reading field errors.
//  Retrieves field errors from the DOOM data base and uses a MPHandler
//  to store them in the beam line.

class MPReader: public MPBase {

public:

    /// Constructor.
    //  Store [b]name[/b] in the DOOM environment.
    MPReader(const std::string &name);

    virtual ~MPReader();

    /// Handle field error.
    virtual void fieldError(const std::string &name, int occur,
                            const BMultipoleField &mainField,
                            BMultipoleField &errorField);

private:

    // Not implemented.
    MPReader();
    MPReader(const MPReader &);
    void operator=(const MPReader &);
};

#endif // OPAL_MPReader_HH
