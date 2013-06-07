#ifndef CLASSIC_CyclotronValleyRep_HH
#define CLASSIC_CyclotronValleyRep_HH

// ------------------------------------------------------------------------
// $RCSfile: CyclotronValleyRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronValleyRep
//   Defines a representation for a Cyclotron Valley.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2010/12/8 14:57:34 $
// $Author: Chuan Wang $
//
// ------------------------------------------------------------------------
#include "AbsBeamline/CyclotronValley.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/AcceleratingField.h"


// Class CyclotronValleyRep
// ------------------------------------------------------------------------
/// Representation for a Cyclotron Valley.

class CyclotronValleyRep: public CyclotronValley {

public:

    /// Constructor with given name.
    explicit CyclotronValleyRep(const string &name);

    CyclotronValleyRep();
    CyclotronValleyRep(const CyclotronValleyRep &);
    virtual ~CyclotronValleyRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual AcceleratingField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const AcceleratingField &getField() const;

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

 private:

    // Not implemented.
    void operator=(const CyclotronValleyRep &);

    /// The cavity's geometry.
    StraightGeometry geometry;

    /// The cavity's field.
    AcceleratingField field;

    /// Cavities are ignored (amplitude = 0) when this switch is set.
    static bool ignoreCavities;
};

#endif // CLASSIC_CyclotronValleyRep_HH
