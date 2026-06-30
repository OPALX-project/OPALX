//
// Class TBeamline
//   Template class for beam lines.
//   Instantiation with different T types allows attachment of additional
//   data to each position in the line.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_TBeamline_HH
#define OPALX_TBeamline_HH

#include "AbsBeamline/BeamlineVisitor.h"
#include "Beamlines/Beamline.h"

#include <algorithm>
#include <list>
#include <string>
#include "Algorithms/Quaternion.hpp"

template <class T>
class TBeamline : public Beamline, public std::list<T> {
public:
    /// Default constructor.
    TBeamline();

    /// Constructor with given name.
    explicit TBeamline(const std::string& name);

    TBeamline(const TBeamline<T>& right);
    virtual ~TBeamline();

    /// Apply BeamlineVisitor to this line.
    virtual void accept(BeamlineVisitor&) const;

    /// Apply visitor to all elements of the line.
    //  If the parameter [b]r2l[/b] is true, the line is traversed from
    //  right (s=C) to left (s=0).
    //  If any error occurs, this method may throw an exception.
    virtual void iterate(BeamlineVisitor&, bool r2l) const;

    /// Make clone.
    virtual TBeamline<T>* clone() const;

    /// Make structure copy.
    virtual TBeamline<T>* copyStructure();

    /// Set sharable flag.
    //  The whole beamline and the elements depending on [b]this[/b] are
    //  marked as sharable.  After this call a [b]copyStructure()[/b] call
    //  reuses the element.
    virtual void makeSharable();

    /// Get geometry.
    //  A placeholder geometry that exists only to satisfy the ElementBase
    //  interface; the beamline's lengths are computed by summing its elements.
    virtual Geometry& getGeometry();
    virtual const Geometry& getGeometry() const;

    /// Get arc length.
    //  Return the length of the geometry, measured along the design orbit.
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the length of the geometry, measured along the design polygone.
    virtual double getElementLength() const;

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the
    //  geometry (non-negative).
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the
    //  geometry (non-positive).
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the
    //  geometry (non-negative).
    virtual double getExit() const;

    /// Get beamline type
    virtual ElementType getType() const;

    /// Append a T object.
    virtual void append(const T&);

    /// Prepend a T object.
    virtual void prepend(const T&);

    void setOrigin3D(const Vector_t<double, 3>& ori);
    Vector_t<double, 3> getOrigin3D() const;
    void setInitialDirection(const Quaternion& rot);
    Quaternion getInitialDirection() const;

    void setRelativeFlag(bool flag);
    bool getRelativeFlag() const;
    /// Get the number of elements in the TBeamline
    size_t size() const;

protected:
    /// Placeholder geometry to satisfy the ElementBase interface.
    Geometry itsGeometry_m;

    Vector_t<double, 3> itsOrigin_m;
    Quaternion itsCoordTrafoTo_m;
    bool relativePositions_m;
};

// Implementation of template class TBeamline
// ------------------------------------------------------------------------

template <class T>
TBeamline<T>::TBeamline()
    : Beamline(),
      std::list<T>(),
      itsOrigin_m(0),
      itsCoordTrafoTo_m(1.0, 0.0, 0.0, 0.0),
      relativePositions_m(false) {}

template <class T>
TBeamline<T>::TBeamline(const std::string& name)
    : Beamline(name),
      std::list<T>(),
      itsOrigin_m(0),
      itsCoordTrafoTo_m(1.0, 0.0, 0.0, 0.0),
      relativePositions_m(false) {}

template <class T>
TBeamline<T>::TBeamline(const TBeamline<T>& rhs)
    : Beamline(rhs),
      std::list<T>(rhs),
      itsOrigin_m(rhs.itsOrigin_m),
      itsCoordTrafoTo_m(rhs.itsCoordTrafoTo_m),
      relativePositions_m(rhs.relativePositions_m) {}

template <class T>
TBeamline<T>::~TBeamline() {}

template <class T>
void TBeamline<T>::accept(BeamlineVisitor& visitor) const {
    visitor.visitBeamline(*this);
}

template <class T>
void TBeamline<T>::iterate(BeamlineVisitor& visitor, bool r2l) const {
    if (r2l) {
        for (typename std::list<T>::const_reverse_iterator op = this->rbegin(); op != this->rend();
             ++op) {
            op->accept(visitor);
        }
    } else {
        for (typename std::list<T>::const_iterator op = this->begin(); op != this->end(); ++op) {
            op->accept(visitor);
        }
    }
}

template <class T>
TBeamline<T>* TBeamline<T>::clone() const {
    TBeamline<T>* line = new TBeamline(getName());

    for (typename std::list<T>::const_iterator op = this->begin(); op != this->end(); ++op) {
        // Make copy of the T object containing a deep copy of its child.
        T newObj(*op);
        newObj.setElement(op->getElement()->clone());
        line->append(newObj);
    }

    line->itsOrigin_m         = itsOrigin_m;
    line->itsCoordTrafoTo_m   = itsCoordTrafoTo_m;
    line->relativePositions_m = relativePositions_m;

    return line;
}

// 21-6-2000 ada change iterator to const_iterator
template <class T>
TBeamline<T>* TBeamline<T>::copyStructure() {
    if (isSharable()) {
        return this;
    } else {
        TBeamline<T>* line = new TBeamline(getName());

        for (typename std::list<T>::const_iterator iter = this->begin(); iter != this->end();
             ++iter) {
            // The copy constructor ensures proper transmission of data.
            T newObj(*iter);
            newObj.setElement(iter->getElement()->copyStructure());
            line->append(newObj);
        }

        line->itsOrigin_m         = itsOrigin_m;
        line->itsCoordTrafoTo_m   = itsCoordTrafoTo_m;
        line->relativePositions_m = relativePositions_m;

        return line;
    }
}

// 21-6-2000 ada change iterator to const_iterator
template <class T>
inline void TBeamline<T>::makeSharable() {
    shareFlag = true;

    for (typename std::list<T>::const_iterator iter = this->begin(); iter != this->end(); ++iter) {
        iter->getElement()->makeSharable();
    }
}

template <class T>
inline Geometry& TBeamline<T>::getGeometry() {
    return itsGeometry_m;
}

template <class T>
inline const Geometry& TBeamline<T>::getGeometry() const {
    return itsGeometry_m;
}

template <class T>
double TBeamline<T>::getArcLength() const {
    double length = 0.0;

    for (typename std::list<T>::const_iterator iter = this->begin(); iter != this->end(); ++iter) {
        length += iter->getElement()->getArcLength();
    }

    return length;
}

template <class T>
double TBeamline<T>::getElementLength() const {
    double length = 0.0;

    for (typename std::list<T>::const_iterator iter = this->begin(); iter != this->end(); ++iter) {
        length += iter->getElement()->getElementLength();
    }

    return length;
}

template <class T>
double TBeamline<T>::getOrigin() const {
    return (getArcLength() / 2.0);
}

template <class T>
double TBeamline<T>::getEntrance() const {
    return (-getOrigin());
}

template <class T>
double TBeamline<T>::getExit() const {
    return (getArcLength() / 2.0);
}


template <class T>
inline ElementType TBeamline<T>::getType() const {
    return ElementType::BEAMLINE;
}

template <class T>
inline void TBeamline<T>::append(const T& obj) {
    this->push_back(obj);
}

template <class T>
inline void TBeamline<T>::prepend(const T& obj) {
    this->push_front(obj);
}

template <class T>
inline void TBeamline<T>::setOrigin3D(const Vector_t<double, 3>& ori) {
    itsOrigin_m = ori;
}

template <class T>
inline Vector_t<double, 3> TBeamline<T>::getOrigin3D() const {
    return itsOrigin_m;
}

template <class T>
inline void TBeamline<T>::setInitialDirection(const Quaternion& trafoTo) {
    itsCoordTrafoTo_m = trafoTo;
}

template <class T>
inline Quaternion TBeamline<T>::getInitialDirection() const {
    return itsCoordTrafoTo_m;
}

template <class T>
inline void TBeamline<T>::setRelativeFlag(bool flag) {
    relativePositions_m = flag;
}

template <class T>
inline bool TBeamline<T>::getRelativeFlag() const {
    return relativePositions_m;
}

template <class T>
size_t TBeamline<T>::size() const {
    size_t blSize = 0;
    for (typename std::list<T>::const_iterator iter = this->begin(); iter != this->end(); ++iter) {
        blSize++;
    }

    return blSize;
}

#endif  // OPALX_TBeamline_HH
