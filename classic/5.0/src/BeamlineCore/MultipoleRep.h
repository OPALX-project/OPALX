#ifndef CLASSIC_MultipoleRep_HH
#define CLASSIC_MultipoleRep_HH

// ------------------------------------------------------------------------
// $RCSfile: MultipoleRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MultipoleRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Multipole.h"
#include "BeamlineGeometry/StraightGeometry.h"


// Class MultipoleRep
// ------------------------------------------------------------------------
/// Representation for a general multipole.

class MultipoleRep: public Multipole {

public:

    /// Constructor with given name.
    explicit MultipoleRep(const string &name);

    MultipoleRep();
    MultipoleRep(const MultipoleRep &);
    virtual ~MultipoleRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b]r and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Set mulitpole field.
    virtual void setField(const BMultipoleField &field);

    /// Allow field errors.
    //  Build a FieldWrapper pointing to the multipole and return a pointer to
    //  that wrapper.
    virtual ElementBase *makeFieldWrapper();

private:

    /// Multipole geometry.
    StraightGeometry geometry;

    /// Multipole field.
    BMultipoleField field;

    // Not implemented.
    void operator=(const MultipoleRep &);
};

#endif // CLASSIC_MultipoleRep_HH
