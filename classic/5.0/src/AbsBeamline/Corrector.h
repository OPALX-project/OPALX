#ifndef CLASSIC_Corrector_HH
#define CLASSIC_Corrector_HH

// ------------------------------------------------------------------------
// $RCSfile: Corrector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Corrector
//   Defines the abstract interface for a orbit corrector.
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
#include "Fields/BDipoleField.h"


// Class Corrector
// ------------------------------------------------------------------------
/// Interface for general corrector.
//  Class Corrector defines the abstract interface for closed orbit
//  correctors.

class Corrector: public Component {

public:

    /// Plane selection.
    enum Plane {
        /// Corrector is off (inactive).
        OFF,
        /// Corrector acts on x-plane.
        X,
        /// Corrector acts on y-plane.
        Y,
        /// Corrector acts on both planes.
        XY
    };

    /// Constructor with given name.
    explicit Corrector(const std::string &name);

    Corrector();
    Corrector(const Corrector &right);
    virtual ~Corrector();

    /// Apply a visitor to Corrector.
    virtual void accept(BeamlineVisitor &) const;


    /// Return the corrector field.
    //  Version for non-constant object.
    virtual BDipoleField &getField() = 0;

    /// Return the corrector field.
    //  Version for constant object.
    virtual const BDipoleField &getField() const = 0;

    /// Return the corrector geometry.
    virtual StraightGeometry &getGeometry() = 0;

    /// Return the corrector geometry. Version for const object.
    virtual const StraightGeometry &getGeometry() const = 0;

    /// Return the plane on which the corrector acts.
    virtual Plane getPlane() const = 0;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void goOnline(const double &kineticEnergy);

    virtual void finalise();

    virtual bool bends() const;

    virtual const std::string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void setKickX(double k);

    void setKickY(double k);

    double getKickX() const;

    double getKickY() const;


private:
    double startField_m;
    double endField_m;
    double kickX_m;
    double kickY_m;

    Vector_t kickField_m;

protected:

    // Not implemented.
    void operator=(const Corrector &);

    Plane plane_m;

    double position_m;
    bool   informed_m;


};

inline
void Corrector::setKickX(double k) {
    kickX_m = k;
}

inline
void Corrector::setKickY(double k) {
    kickY_m = k;
}

inline
double Corrector::getKickX() const {
    return kickX_m;
}

inline
double Corrector::getKickY() const {
    return kickY_m;
}

#endif // CLASSIC_Corrector_HH