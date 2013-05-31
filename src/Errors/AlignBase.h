#ifndef OPAL_AlignBase_HH
#define OPAL_AlignBase_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignBase
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <string>

using std::string;

class AlignWrapper;


// Class AlignBase
// ------------------------------------------------------------------------
/// Handle alignment errors.
//  An abstract mixin class, used for all classes which must have access
//  to misalignment errors for setting or retrieving.

class AlignBase {

public:

    AlignBase();
    virtual ~AlignBase();

    /// Deal with misalignment.
    //  This is the interface method ``mixed in''.
    virtual void misalignment(const AlignWrapper &, int occur) = 0;

private:

    // Not implemented.
    AlignBase(const AlignBase &);
    void operator=(const AlignBase &);
};

#endif // OPAL_AlignBase_HH
