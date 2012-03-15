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
#include "Algorithms/PartBunch.h"

class PartData;
//class PartBunch;
template <class T, int N> class FVps;


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
  explicit Component(const string &name);

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

  virtual bool getFieldstrength(double R[], double t, double E[], double B[]) const;

  virtual bool getFieldstrength(Vector_t R, double t, Vector_t &E, Vector_t &B) const;

  virtual bool readFieldMap(double &startField, double &endField, double scaleFactor);

  virtual void rescaleFieldMap(double scaleFactor);

  /**
     Methods for the cyclotron cmd.

  */

  virtual double getRinit() { } ;
  virtual double getPRinit() { } ;
  virtual double getPHIinit() { } ;

  virtual string getFieldMapFN() { } ;
  virtual double getRfFrequ() { } ;
  virtual double getSymmetry() { } ;

  virtual string getType() { } ;
  virtual double getCyclHarm() { } ;
  virtual void readFieldMap(double scaleFactor); 

  virtual double getRmax() { };
  virtual double getRmin() { };

  virtual void setComponentType(string name){};
  virtual string getComponentType()const{};

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
  virtual void trackMap(FVps<double,6> &map, const PartData &,
			bool revBeam, bool revTrack) const;
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

#endif // CLASSIC_Component_HH
