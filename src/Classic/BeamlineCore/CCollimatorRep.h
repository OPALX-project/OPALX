#ifndef CLASSIC_CCollimatorRep_HH
#define CLASSIC_CCollimatorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: CCollimatorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CCollimatorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/CCollimator.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"


// Class CCollimatorRep
// ------------------------------------------------------------------------
/// Representation for a collimator.

class CCollimatorRep: public CCollimator {

public:

    /// Constructor with given name.
    explicit CCollimatorRep(const std::string &name);

    CCollimatorRep();
    CCollimatorRep(const CCollimatorRep &);
    virtual ~CCollimatorRep();

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

    /*
    /// Return the horizontal half-aperture.
    virtual double getXsize() = 0;//const;

    /// Return the vertical half-aperture.
    virtual double getYsize() = 0; const;

    /// Return the horizontal half-aperture.
    virtual void setXsize(double) = 0;

    /// Return the vertical half-aperture.
    virtual void setYsize(double) = 0;
    */
private:

    // Not implemented.
    void operator=(const CCollimatorRep &);

    // The zero magnetic field.
    NullField field;

    // The geometry.
    StraightGeometry geometry;

    // The horizontal half-size.
    double xSize;

    // The vertical half-size.
    double ySize;
};

#endif // CLASSIC_CCollimatorRep_HH