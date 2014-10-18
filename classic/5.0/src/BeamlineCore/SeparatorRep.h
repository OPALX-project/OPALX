#ifndef CLASSIC_SeparatorRep_HH
#define CLASSIC_SeparatorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: SeparatorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SeparatorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Separator.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/EDipoleField.h"


// Class SeparatorRep
// ------------------------------------------------------------------------
/// Representation for an electrostatic separator.

class SeparatorRep: public Separator {

public:

    /// Constructor with given name.
    explicit SeparatorRep(const std::string &name);

    SeparatorRep();
    SeparatorRep(const SeparatorRep &);
    virtual ~SeparatorRep();

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
    virtual EDipoleField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const EDipoleField &getField() const;

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

    /// Get component.
    //  Return the x-component of the electric field in A/m.
    virtual double getEx() const;

    /// Get component.
    //  Return the y-component of the electric field in A/m.
    virtual double getEy() const;

    /// Set component.
    //  Assign the x-component of the electric field in A/m.
    virtual void setEx(double Ex);

    /// Set component.
    //  Assign the y-component of the electric field in A/m.
    virtual void setEy(double Ey);

private:

    // Not implemented.
    void operator=(const SeparatorRep &);

    /// The separator geometry.
    StraightGeometry geometry;

    /// The separator field.
    EDipoleField field;
};

#endif // CLASSIC_SeparatorRep_HH
