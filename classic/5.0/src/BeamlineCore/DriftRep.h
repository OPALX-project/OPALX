#ifndef CLASSIC_DriftRep_HH
#define CLASSIC_DriftRep_HH

// ------------------------------------------------------------------------
// $RCSfile: DriftRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DriftRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Drift.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"


// Class DriftRep
// ------------------------------------------------------------------------
/// Representation for a drift space.

class DriftRep: public Drift {

public:

    /// Constructor with given name.
    explicit DriftRep(const string &name);

    DriftRep();
    DriftRep(const DriftRep &);
    virtual ~DriftRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey);

    /// Get field.
    //  Version for non-constant object.
    virtual NullField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const NullField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

private:

    // Not implemented.
    void operator=(const DriftRep &);

    /// The zero magnetic field.
    NullField field;

    /// The geometry.
    StraightGeometry geometry;
};

#endif // CLASSIC_DriftRep_HH
