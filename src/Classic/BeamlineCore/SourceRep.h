#ifndef CLASSIC_SOURCEREP_HH
#define CLASSIC_SOURCEREP_HH

// ------------------------------------------------------------------------
// $RCSfile: SourceRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SourceRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Source.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"


// Class SourceRep
// ------------------------------------------------------------------------
/// Tepresentation for a solenoid magnet.

class SourceRep: public Source {

public:

    /// Constructor with given name.
    explicit SourceRep(const std::string &name);

    SourceRep();
    SourceRep(const SourceRep &);
    virtual ~SourceRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

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
    void operator=(const SourceRep &);

    // The zero magnetic field.
    NullField field;

    /// The solenoid geometry.
    StraightGeometry geometry;
};

#endif // CLASSIC_SOURCEREP_HH