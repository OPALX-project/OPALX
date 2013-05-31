#ifndef OPAL_MPRemover_HH
#define OPAL_MPRemover_HH

// ------------------------------------------------------------------------
// $RCSfile: MPRemover.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPRemover
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPBase.h"


// Class MPRemover
// ------------------------------------------------------------------------
/// A visitor used to remove all random field errors from a line.
//  Uses a MPHandler to access all elements which have a field error,
//  and sets all field errors to zero.

class MPRemover: public MPBase {

public:

    MPRemover();
    virtual ~MPRemover();

    /// Remove error field.
    virtual void fieldError(const std::string &name,
                            int occur,
                            const BMultipoleField &mainField,
                            BMultipoleField &errorField);

private:

    // Not implemented.
    MPRemover(const MPRemover &);
    void operator=(const MPRemover &);
};

#endif // OPAL_MPRemover_HH


