#ifndef CLASSIC_BeamBeam_HH
#define CLASSIC_BeamBeam_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamBeam.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
//   Defines the abstract interface for a beam-beam interaction.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


class Matrix3D;
class Vector3D;


// Class BeamBeam
// ------------------------------------------------------------------------
/// Abstract beam-beam interaction.
//  Class BeamBeam defines the abstract interface for beam-beam interactions.


class BeamBeam: public Component {

public:

    /// Constructor with given name.
    explicit BeamBeam(const string &name);

    BeamBeam();
    BeamBeam(const BeamBeam &right);
    virtual ~BeamBeam();

    /// Apply visitor to BeamBeam.
    virtual void accept(BeamlineVisitor &) const;

    /// Get bunch charge.
    //  Return the number of particles times the particle charge in the
    //  opposite bunch.  Units are proton charges.
    virtual double getBunchCharge() const = 0;

    /// Get moments.
    //  Return the moment matrix for the opposite bunch
    //  (the matrix of second momenta).  Units are square metres.
    virtual const Matrix3D &getBunchMoment() const = 0;

    /// Get displacement.
    //  Return the displacement vector for position of opposite bunch.
    //  Units are metres.
    virtual const Vector3D &getBunchDisplacement() const = 0;

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual const string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

    // Not implemented.
    void operator=(const BeamBeam &);
};

#endif // CLASSIC_BeamBeam_HH
