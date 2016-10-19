#ifndef CLASSIC_SolenoidRep_HH
#define CLASSIC_SolenoidRep_HH

// ------------------------------------------------------------------------
// $RCSfile: SolenoidRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SolenoidRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Solenoid.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/ConstBzField.h"


// Class SolenoidRep
// ------------------------------------------------------------------------
/// Tepresentation for a solenoid magnet.

class SolenoidRep: public Solenoid {

public:

    /// Constructor with given name.
    explicit SolenoidRep(const std::string &name);

    SolenoidRep();
    SolenoidRep(const SolenoidRep &);
    virtual ~SolenoidRep();

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
    virtual ConstBzField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const ConstBzField &getField() const;

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

    /// Get field.
    //  Return the solenoid field in Teslas.
    virtual double getBz() const;

    /// Set field.
    //  Assign the solenoid field in Teslas.
    virtual void setBz(double Bz);

private:

    // Not implemented.
    void operator=(const SolenoidRep &);

    /// The solenoid geometry.
    StraightGeometry geometry;

    /// The solenoid field.
    ConstBzField field;
};

#endif // CLASSIC_SolenoidRep_HH
