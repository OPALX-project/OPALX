#ifndef OPAL_MPBase_HH
#define OPAL_MPBase_HH

// ------------------------------------------------------------------------
// $RCSfile: MPBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPBase
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/14 07:08:27 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include <string>
using namespace std;

class BMultipoleField;


// Class MPBase
// ------------------------------------------------------------------------
/// Handle multipole errors.
//  An abstract mixin class, used for all classes which must have access
//  to field errors for setting or retrieving.

class MPBase {

public:

    MPBase();
    virtual ~MPBase();

    /// Deal with error field.
    //  This is the interface method ``mixed in''.
    virtual void fieldError(const string &name, int occur,
                            const BMultipoleField &designField,
                            BMultipoleField &errorField) = 0;

    /// Scale factor for absolute error components.
    //  Return $(p0/c)/i!$.
    static double absFactor(int i);

    /// Scale factor for relative error components.
    //  Return $r^(o-c)$.
    static double relFactor(int i, int o, double r);

private:

    // Not implemented.
    MPBase(const MPBase &);
    void operator=(const MPBase &);
};

#endif // OPAL_MPBase_HH









