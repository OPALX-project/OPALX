#ifndef CLASSIC_MultipoleWrapper_HH
#define CLASSIC_MultipoleWrapper_HH

// ------------------------------------------------------------------------
// $RCSfile: MultipoleWrapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MultipoleWrapper
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Multipole.h"
#include "Fields/BMultipoleField.h"
#include "MemoryManagement/Pointer.h"


// Class MultipoleWrapper
// ------------------------------------------------------------------------
/// Representation of a perturbed multipole.
//  A MultipoleWrapper represents a unique instance of a multipole magnet
//  in the accelerator model. It defines imperfections of the field,
//  related to an ``ideal'' magnet contained in the wrapper.

class MultipoleWrapper: public Multipole {

public:

    /// Constructor.
    //  Constructs a wrapper for its argument.
    //  The wrapper is not sharable by default.
    explicit MultipoleWrapper(Multipole *);

    MultipoleWrapper(const MultipoleWrapper &);
    virtual ~MultipoleWrapper();


    /// Apply visitor to modified multipole.
    virtual void accept(BeamlineVisitor &) const;

    /// Make clone.
    virtual ElementBase *clone() const;

    /// Make structural copy.
    virtual ElementBase *copyStructure();

    /// Get multipole field error.
    //  This method can be used to get or set the error field. The error field
    //  offset is mutable, so as to allow changing it in a constant structure.
    virtual BMultipoleField &errorField() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Get field.
    //  Return the perturbed field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Return the perturbed field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get element type string.
    virtual ElementBase::ElementType getType() const;

    /// Get design corrector.
    //  Version for constant object.
    virtual const Multipole &getDesign() const;

    /// Get design corrector.
    //  Version for non-constant object.
    virtual Multipole &getDesign();

    /// Make wrapper for this multipole.
    //  Return [b]this[/b], since this is already a field wrapper.
    virtual ElementBase *makeFieldWrapper();

    /// Remove field wrapper.
    virtual ElementBase *removeFieldWrapper();

    /// Remove field wrapper.
    virtual const ElementBase *removeFieldWrapper() const;

    /// Remove all wrappers.
    virtual ElementBase *removeWrappers();

    /// Remove all wrappers.
    virtual const ElementBase *removeWrappers() const;

private:

    // Not implemented.
    MultipoleWrapper();
    void operator=(const MultipoleWrapper &);

    // The pointer to the ``design'' magnet.
    Pointer<Multipole> itsDesign;

    // The error field.
    mutable BMultipoleField itsError;

    // The magnetic field after application of all modifiers.
    mutable BMultipoleField tempField;
};

#endif // CLASSIC_MultipoleWrapper_HH