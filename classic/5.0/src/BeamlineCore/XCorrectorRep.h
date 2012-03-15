#ifndef CLASSIC_XCorrectorRep_HH
#define CLASSIC_XCorrectorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: XCorrectorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: XCorrectorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CorrectorRep.h"


// Class XCorrectorRep
// ------------------------------------------------------------------------
/// Representation for an orbit corrector.
//  This derived class acts on the horizontal plane.

class XCorrectorRep: public CorrectorRep {

public:

    /// Constructor with given name.
    explicit XCorrectorRep(const string &name);

    XCorrectorRep();
    XCorrectorRep(const XCorrectorRep &);
    virtual ~XCorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey);

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get plane of action.
    //  Return the x-plane for this class.
    virtual Plane getPlane() const;

    /// Get field.
    //  Return horizontal component (always zero).
    virtual double getBx() const;

    /// Set field.
    //  Ignore the horizontal field value.
    virtual void setBx(double);

private:

    // Not implemented.
    void operator=(const XCorrectorRep &);
};

#endif // CLASSIC_XCorrectorRep_HH
