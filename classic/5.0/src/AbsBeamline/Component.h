#ifndef CLASSIC_Component_HH
#define CLASSIC_Component_HH

// ------------------------------------------------------------------------
// $RCSfile: Component.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Component
//   An abstract base class which defines the common interface for all
//   CLASSIC components, i.e. beam line members which are not themselves
//   beam lines.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "Fields/EMField.h"
#include "Algorithms/Vektor.h"

#define EPS_MISALIGNMENT 1e-8

class PartData;
class PartBunch;
template <class T, int N> class FVps;

struct Point
{
  double x;
  double y;
};

// Class Component
// ------------------------------------------------------------------------
/// Interface for a single beam element.
//  Class Component defines the abstract interface for an arbitrary single
//  component in a beam line.  A component is the basic element in the
//  accelerator model, like a dipole, a quadrupole, etc.  It is normally
//  associated with an electro-magnetic field, which may be null.

class Component: public ElementBase {

public:

    /// Constructor with given name.
    explicit Component(const std::string &name);

    Component();
    Component(const Component &right);
    virtual ~Component();

    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for non-constant object).
    virtual EMField &getField() = 0;

    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for constant object).
    virtual const EMField &getField() const = 0;

    /// Return the field in a point.
    //  Return the value of the time-independent part of the electric
    //  field at point [b]P[/b].
    EVector Efield(const Point3D &P) const;

    /// Return the field in a point.
    //  Return the value of the time-independent part of the magnetic
    //  field at point [b]P[/b].
    BVector Bfield(const Point3D &P) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of the electric
    //  field at point [b]P[/b] for time [b]t[/b].
    EVector Efield(const Point3D &P, double t) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of the magnetic
    //  field at point [b]P[/b] for time [b]t[/b].
    BVector Bfield(const Point3D &P, double t) const;

    /// Return the field in a point.
    //  Return the value of the time-independent part of both electric
    //  and magnetic fields at point [b]P[/b].
    EBVectors EBfield(const Point3D &P) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of both electric
    //  and magnetic fields at point [b]P[/b] for time [b]t[/b].
    EBVectors EBfield(const Point3D &P, double t) const;

    virtual void addKR(int i, double t, Vector_t &K) {};

    virtual void addKT(int i, double t, Vector_t &K) {};

    virtual bool apply(const int &i, const double &t, double E[], double B[]) = 0;

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B) = 0;

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) = 0;

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) = 0;

    virtual void finalise() = 0;

    virtual bool bends() const = 0;

    virtual void goOnline();

    virtual void goOffline();

    virtual bool Online();

    void getOrientation(Vector_t &, double &) const;  // first component is alpha, second is beta. the third component is always neglected.
    // alpha is defined as the angle between projection of the normal of the face of the element onto
    // the s-u plane and the s vector. for beta the following is valid:
    // beta = arccos(|n_{parallel}|) where n_{parallel} is the projection of the normal onto the s-u
    // plane

    void setOrientation(const Vector_t &direction);

    virtual void getDimensions(double &zBegin, double &zEnd) const = 0;

    virtual const std::string &getType() const;

    virtual void setComponentType(std::string name) { };
    virtual std::string getComponentType() const { return ""; };

    /// Return design element.
    //  If a component is a wrapper, this method returns a pointer to
    //  its underlying design element, otherwise a pointer to this component.
    //  The default version returns ``this''.
    virtual const ElementBase &getDesign() const;

    /// Track particle bunch.
    //  This catch-all method implements a hook for tracking a particle
    //  bunch through a non-standard component.
    //  The default version throws a LogicalError.
    virtual void trackBunch(PartBunch &bunch, const PartData &,
                            bool revBeam, bool revTrack) const;

    /// Track a map.
    //  This catch-all method implements a hook for tracking a transfer
    //  map through a non-standard component.
    //  The default version throws a LogicalError.
    virtual void trackMap(FVps<double, 6> &map, const PartData &,
                          bool revBeam, bool revTrack) const;

    void setMisalignment(double x, double y, double s);

    void getMisalignment(double &x, double &y, double &s) const;

    void setExitFaceSlope(const double &);

protected:

    double dx_m;
    double dy_m;
    double ds_m;
    Vector_t Orientation_m;
    double exit_face_slope_m;

    PartBunch *RefPartBunch_m;
    bool online_m;
};


// Inline access functions to fields.
// ------------------------------------------------------------------------

inline EVector Component::Efield(const Point3D &P) const
{ return getField().Efield(P); }

inline BVector Component::Bfield(const Point3D &P) const
{ return getField().Bfield(P); }

inline EVector Component::Efield(const Point3D &P, double t) const
{ return getField().Efield(P, t); }

inline BVector Component::Bfield(const Point3D &P, double t) const
{ return getField().Bfield(P, t); }

inline EBVectors Component::EBfield(const Point3D &P) const
{ return getField().EBfield(P); }

inline EBVectors Component::EBfield(const Point3D &P, double t) const
{ return getField().EBfield(P, t); }

inline void Component::getOrientation(Vector_t &ori, double &m) const {
    ori = Orientation_m;
    m = exit_face_slope_m;
}

inline void Component::setOrientation(const Vector_t &direction)
{ Orientation_m = direction; }

inline void Component::setExitFaceSlope(const double &m)
{ exit_face_slope_m = m; }
#endif // CLASSIC_Component_HH
