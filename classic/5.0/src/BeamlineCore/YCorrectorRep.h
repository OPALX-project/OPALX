#ifndef CLASSIC_YCorrectorRep_HH
#define CLASSIC_YCorrectorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: YCorrectorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: YCorrectorRep
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


// Class YCorrectorRep
// ------------------------------------------------------------------------
/// Representation for an orbit corrector.
//  Acts on the vertical plane.

class YCorrectorRep: public CorrectorRep {

public:

    /// Constructor with given name.
    explicit YCorrectorRep(const string &name);

    YCorrectorRep();
    YCorrectorRep(const YCorrectorRep &);
    virtual ~YCorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get plane.
    //  Return the y-plane for this class.
    virtual Plane getPlane() const;

    /// Get field.
    //  Return the vertical field component (always zero).
    virtual double getBy() const;

    /// Set field.
    //  Ignore the vertical field component.
    virtual void setBy(double);

private:

    // Not implemented.
    void operator=(const YCorrectorRep &);
};

#endif // CLASSIC_YCorrectorRep_HH
