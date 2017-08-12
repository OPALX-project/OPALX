#ifndef CLASSIC_Degrader_HH
#define CLASSIC_Degrader_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//   *** MISSING *** Degrader interface is still incomplete.
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
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include <vector>
#include <tuple>

class LossDataSink;

// Class Degrader
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class Degrader defines the abstract interface for a collimator.

class Degrader: public Component {

public:
    /// Constructor with given name.
    explicit Degrader(const std::string &name);

    Degrader();
    Degrader(const Degrader &rhs);
    virtual ~Degrader();

    /// Apply visitor to Degrader.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    virtual bool isInMaterial(const Vector_t & R);

    void defineEllipticShape(double M, double m);
private:

    // Not implemented.
    void operator=(const Degrader &);

    virtual bool isInMaterialLong(const Vector_t &R);
    virtual bool isInMaterialTrans(const Vector_t &R);

    double position_m;
    double semiMinorAxis_m;
    double semiMajorAxis_m;
};

inline
bool Degrader::isInMaterial(const Vector_t & R)
{
 /**
     check if the particle is in the degarder material

  */
    return isInMaterialLong(R) && isInMaterialTrans(R);
}

inline
bool Degrader::isInMaterialLong(const Vector_t & R)
{
    return (R(2) > position_m &&
            R(2) <= position_m + getElementLength());
}

inline
bool Degrader::isInMaterialTrans(const Vector_t & R)
{
    double r = std::pow(R(0) / semiMajorAxis_m, 2) + std::pow(R(1) / semiMinorAxis_m, 2);
    return r < 1.0;
}

inline
void Degrader::defineEllipticShape(double M, double m)
{
    semiMinorAxis_m = m;
    semiMajorAxis_m = M;
}

#endif // CLASSIC_Degrader_HH