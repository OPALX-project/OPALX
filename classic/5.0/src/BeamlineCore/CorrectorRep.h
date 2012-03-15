#ifndef CLASSIC_CorrectorRep_HH
#define CLASSIC_CorrectorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: CorrectorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CorrectorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Corrector.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BDipoleField.h"


// Class CorrectorRep
// ------------------------------------------------------------------------
/// Representation of a closed orbit corrector.
//  The base class acts on both planes.

class CorrectorRep: public Corrector {

public:

    /// Constructor with given name.
    explicit CorrectorRep(const string &name);

    CorrectorRep();
    CorrectorRep(const CorrectorRep &right);
    virtual ~CorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

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

    /// Get plane(s) of action.
    virtual Plane getPlane() const;

    /// Get horizontal field component in Teslas.
    virtual double getBx() const;

    /// Get vertical field component in Teslas.
    virtual double getBy() const;

    /// Get corrector field.
    virtual BDipoleField &getField();

    /// Get corrector field. Version for const corrector.
    virtual const BDipoleField &getField() const;

    /// Set horizontal field component in Teslas.
    virtual void setBx(double);

    /// Set vertical field component in Teslas.
    virtual void setBy(double);

    /// Allow field errors.
    //  Build a FieldWrapper pointing to the corrector and return a pointer to
    //  that wrapper.
    virtual ElementBase *makeFieldWrapper();

    /// Set active flag.
    //  If [b]flag[/b] is true, the corrector is activated,
    //  otherwise deactivated.
    virtual void setActive(bool flag = true);

protected:

    /// The corrector geometry.
    StraightGeometry geometry;

    /// The corrector strengths.
    BDipoleField field;

    /// The active/inactive flag.
    bool active;

private:

    // Not implemented.
    void operator=(const CorrectorRep &);
};

#endif // CLASSIC_CorrectorRep_HH
