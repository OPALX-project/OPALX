#ifndef CLASSIC_ParallelPlateRep_HH
#define CLASSIC_ParallelPlateRep_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelPlateRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelPlateRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ParallelPlate.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/AcceleratingField.h"


// Class ParallelPlateRep
// ------------------------------------------------------------------------
/// Representation for a RF cavity.

class ParallelPlateRep: public ParallelPlate {

public:

    /// Constructor with given name.
    explicit ParallelPlateRep(const string &name);

    ParallelPlateRep();
    ParallelPlateRep(const ParallelPlateRep &);
    virtual ~ParallelPlateRep();

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
    void operator=(const ParallelPlateRep &);

    /// The cavity's geometry.
    StraightGeometry geometry;

    /// The cavity's field.
    AcceleratingField field;

  
};

#endif // CLASSIC_ParallelPlateRep_HH
