#ifndef CLASSIC_RFCavityRep_HH
#define CLASSIC_RFCavityRep_HH

// ------------------------------------------------------------------------
// $RCSfile: RFCavityRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavityRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFCavity.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/AcceleratingField.h"


// Class RFCavityRep
// ------------------------------------------------------------------------
/// Representation for a RF cavity.

class RFCavityRep: public RFCavity {

public:

    /// Constructor with given name.
    explicit RFCavityRep(const string &name);

    RFCavityRep();
    RFCavityRep(const RFCavityRep &);
    virtual ~RFCavityRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

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

    /// Get amplitude.
    //  Return the RF amplitude in Volts.
    virtual double getAmplitude() const;

    /// Get frequency.
    //  Return the RF frequency in Hertz.
    virtual double getFrequency() const;

    /// Get phase.
    //  Return the RF phase in radians.
    virtual double getPhase() const;

    /// Set amplitude.
    //  Assign the RF amplitude in Volts.
    virtual void setAmplitude(double V);

    /// Set frequency.
    //  Assign the RF frequency in Hertz.
    virtual void setFrequency(double f);

    /// Set phase.
    //  Assing the RF phase in radians.
    virtual void setPhase(double phi);

    /// Set ignore switch.
    //  If [b]ignore[/b] is true, cavities are ignored for static calculations.
    static void setIgnore(bool ignore = false);

private:

    // Not implemented.
    void operator=(const RFCavityRep &);

    /// The cavity's geometry.
    StraightGeometry geometry;

    /// The cavity's field.
    AcceleratingField field;

    /// Cavities are ignored (amplitude = 0) when this switch is set.
    static bool ignoreCavities;
};

#endif // CLASSIC_RFCavityRep_HH
