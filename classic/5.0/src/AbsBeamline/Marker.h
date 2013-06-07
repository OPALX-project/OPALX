#ifndef CLASSIC_Marker_HH
#define CLASSIC_Marker_HH

// ------------------------------------------------------------------------
// $RCSfile: Marker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines the abstract interface for a marker element.
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


// Class Marker
// ------------------------------------------------------------------------
/// Interface for a marker.
//  Class Marker defines the abstract interface for a marker element.

class Marker: public Component {

public:

    /// Constructor with given name.
    explicit Marker(const string &name);

    Marker();
    Marker(const Marker &);
    virtual ~Marker();

    /// Apply visitor to Marker.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual const string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

    // Not implemented.
    void operator=(const Marker &);
};

#endif // CLASSIC_Marker_HH
