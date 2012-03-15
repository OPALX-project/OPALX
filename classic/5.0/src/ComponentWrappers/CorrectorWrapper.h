#ifndef CLASSIC_CorrectorWrapper_HH
#define CLASSIC_CorrectorWrapper_HH

// ------------------------------------------------------------------------
// $RCSfile: CorrectorWrapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CorrectorWrapper
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Corrector.h"
#include "Fields/BDipoleField.h"
#include "MemoryManagement/Pointer.h"

class BDipoleField;


// Class CorrectorWrapper
// ------------------------------------------------------------------------
/// Representation for a perturbed closed orbit corrector.
//  A CorrectorWrapper represents a unique instance of a corrector magnet
//  in the accelerator model. It defines imperfections of the field,
//  related to an ``ideal'' magnet contained in the wrapper.

class CorrectorWrapper: public Corrector {

public:

    /// Constructor.
    //  Construct a wrapper for the argument.
    //  The wrapper is not sharable by default.
    CorrectorWrapper(Corrector *);

    CorrectorWrapper(const CorrectorWrapper &);
    virtual ~CorrectorWrapper();


    /// Apply visitor to modified corrector.
    virtual void accept(BeamlineVisitor &) const;

    /// Make clone.
    virtual ElementBase *clone() const;

    /// Make structural copy.
    virtual ElementBase *copyStructure();

    /// Get corrector field error.
    //  This method can be used to get or set the error field. The error field
    //  offset is mutable, so as to allow changing it in a constant structure.
    virtual BDipoleField &errorField() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Get corrector field.
    //  Version for non-constant object.
    virtual BDipoleField &getField();

    /// Get corrector field.
    //  Version for constant object.
    virtual const BDipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get plane of action.
    virtual Plane getPlane() const;

    /// Get element type string.
    virtual const string &getType() const;

    /// Get design corrector.
    //  Version for constant object.
    virtual const Corrector &getDesign() const;

    /// Get design corrector.
    //  Version for non-constant object.
    virtual Corrector &getDesign();

    /// Make wrapper for this corrector.
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
    CorrectorWrapper();
    void operator=(const CorrectorWrapper &);

    // The pointer to the ``design'' magnet.
    Pointer<Corrector> itsDesign;

    /// The field modifiers.
    mutable BDipoleField itsError;

    // The magnetic field after application of error.
    mutable BDipoleField tempField;
};

#endif // CLASSIC_CorrectorWrapper_HH
