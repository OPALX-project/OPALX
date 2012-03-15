#ifndef OPAL_OpalSolenoid_HH
#define OPAL_OpalSolenoid_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSolenoid.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSolenoid
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalSolenoid
// ------------------------------------------------------------------------
/// The SOLENOID element.

class OpalSolenoid: public OpalElement {

public:

    /// The attributes of class OpalSolenoid.
    enum {
        KS = COMMON,  // The longitudinal magnetic field.
        FMAPFN,          // The Field filename
        FAST,           // Faster but less accurate
        DX,             // Misalignment: translation in x direction
        DY,             // Misalignment: translation in y direction
        DZ,             // Misalignment: translation in z direction
        SIZE
    };

    /// Exemplar constructor.
    OpalSolenoid();

    virtual ~OpalSolenoid();

    /// Make clone.
    virtual OpalSolenoid *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC solenoid.
    virtual void update();

private:

    // Not implemented.
    OpalSolenoid(const OpalSolenoid &);
    void operator=(const OpalSolenoid &);

    // Clone constructor.
    OpalSolenoid(const string &name, OpalSolenoid *parent);
};

#endif // OPAL_OpalSolenoid_HH

