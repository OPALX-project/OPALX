#ifndef CLASSIC_Lambertson_HH
#define CLASSIC_Lambertson_HH

// ------------------------------------------------------------------------
// $RCSfile: Lambertson.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   *** MISSING *** Lambertson interface is still incomplete.
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


// Class Lambertson
// ------------------------------------------------------------------------
/// Interface for a Lambertson septum.
//  Class Lambertson defines the abstract interface for a Lambertson
//  septum magnet.

class Lambertson: public Component {

public:

    /// Constructor with given name.
    explicit Lambertson(const std::string &name);

    Lambertson();
    Lambertson(const Lambertson &);
    virtual ~Lambertson();

    /// Apply visitor to Lambertson.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

    // Not implemented.
    void operator=(const Lambertson &);
};

#endif // CLASSIC_Lambertson_HH
