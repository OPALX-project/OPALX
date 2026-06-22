//
// Class OpalBeamline
//   :FIXME: add class description
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
#ifndef OPAL_BEAMLINE_H
#define OPAL_BEAMLINE_H

#include <memory>
#include <set>
#include <string>

#include "AbsBeamline/Component.h"
#include "PartBunch/PartBunch.h"

#include "AbsBeamline/Marker.h"
#include "Beamlines/Beamline.h"
#include "Utilities/BeamlineFieldElement.h"

#include "Algorithms/CoordinateSystemTrafo.h"

#include "OPALTypes.h"

class ParticleMatterInteractionHandler;
class BoundaryGeometry;

class OpalBeamline {
public:
    OpalBeamline();
    OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation);
    ~OpalBeamline();

    void activateElements();
    std::set<std::shared_ptr<Component>> getElements(const Vector_t<double, 3>& x);

    /**
     * Get all elements in the beamline, regardless of their position.
     * @return Set of shared pointers to all elements in the beamline.
     */
    std::set<std::shared_ptr<Component>> getElements();

    Vector_t<double, 3> transformTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> transformFrom(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateFrom(const Vector_t<double, 3>& r) const;

    Vector_t<double, 3> transformToLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> transformFromLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateToLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateFromLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;

    /// Return the element's nominal body transform (global→local).
    CoordinateSystemTrafo getCSTrafoLab2Local(const std::shared_ptr<Component>& comp) const;
    CoordinateSystemTrafo getCSTrafoLab2Local() const;
    CoordinateSystemTrafo getMisalignment(const std::shared_ptr<Component>& comp) const;
    CoordinateSystemTrafo getNominalEntryTransform(const std::shared_ptr<Component>& comp) const;
    CoordinateSystemTrafo getNominalExitTransform(const std::shared_ptr<Component>& comp) const;

    double getStart(const Vector_t<double, 3>&) const;
    double getEnd(const Vector_t<double, 3>&) const;

    void switchElements(
            const double&, const double&, const double& kineticEnergy,
            const bool& nomonitors = false);
    void switchElementsOff();

    ParticleMatterInteractionHandler* getParticleMatterInteractionHandler(const unsigned int&);

    BoundaryGeometry* getBoundaryGeometry(const unsigned int&);

    unsigned long getFieldAt(
            const unsigned int&, const Vector_t<double, 3>&, const long&, const double&,
            Vector_t<double, 3>&, Vector_t<double, 3>&);
    unsigned long getFieldAt(
            const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
            Vector_t<double, 3>&, Vector_t<double, 3>&);

    template <class T>
    void visit(const T&, BeamlineVisitor&, PartBunch_t&);

    void prepareSections();
    void positionElementRelative(std::shared_ptr<ElementBase>);
    void compute3DLattice();
    void save3DLattice();
    void save3DInput();
    void print(Inform&) const;
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    FieldList getElementByType(ElementType);

    void swap(OpalBeamline& rhs);
    void merge(OpalBeamline& rhs);

private:
    /// Set one element's nominal body transform (its global→local
    /// CoordinateSystemTrafo) during beamline assembly.
    void setNominalPlacement(
            const std::shared_ptr<ElementBase>& element, const CoordinateSystemTrafo& parentToBody);

    /**
     * @brief Place ELEMEDGE-positioned elements along the reference path.
     *
     * An element placed with `ELEMEDGE` is stored only as a path-length position
     * `s` (plus an optional roll about the beam axis); its lab-frame pose is not
     * known until the full lattice is sorted by `s`. This setup-stage pass walks
     * the sorted elements once, accumulates the running reference-path transform
     * (including bends), and writes each element's nominal rigid placement
     * \f$T_i\f$. Elements already fixed in place by a 6D pose (X, Y, Z, THETA,
     * PHI, PSI) are skipped.
     */
    void placeElementsAlongReferencePath();

    FieldList elements_m;
    bool prepared_m;
    bool referencePathPlacementCompiled_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

template <class T>
inline void OpalBeamline::visit(const T& element, BeamlineVisitor&, PartBunch_t& bunch) {
    Inform msg("OPAL ");
    double startField = 0.0;
    double endField   = 0.0;
    std::shared_ptr<T> elptr(dynamic_cast<T*>(element.clone()));

    positionElementRelative(elptr);

    if (elptr->isElementPositionSet()) startField = elptr->getElementPosition();

    elptr->initialise(&bunch, startField, endField);
    elements_m.push_back(BeamlineFieldElement(elptr, startField, endField));
    prepared_m                       = false;
    referencePathPlacementCompiled_m = false;
}

template <>
inline void OpalBeamline::visit<Marker>(const Marker& /*element*/, BeamlineVisitor&, PartBunch_t&) {
}

inline Vector_t<double, 3> OpalBeamline::transformTo(const Vector_t<double, 3>& r) const {
    return coordTransformationTo_m.transformTo(r);
}

inline Vector_t<double, 3> OpalBeamline::transformFrom(const Vector_t<double, 3>& r) const {
    return coordTransformationTo_m.transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateTo(const Vector_t<double, 3>& r) const {
    return coordTransformationTo_m.rotateTo(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFrom(const Vector_t<double, 3>& r) const {
    return coordTransformationTo_m.rotateFrom(r);
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local(
        const std::shared_ptr<Component>& comp) const {
    return comp->getCSTrafoGlobal2Local();
}

inline Vector_t<double, 3> OpalBeamline::transformToLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).transformTo(r);
}

inline Vector_t<double, 3> OpalBeamline::transformFromLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).rotateTo(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).rotateFrom(r);
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline CoordinateSystemTrafo OpalBeamline::getMisalignment(
        const std::shared_ptr<Component>& comp) const {
    return comp->getMisalignment();
}

inline CoordinateSystemTrafo OpalBeamline::getNominalEntryTransform(
        const std::shared_ptr<Component>& comp) const {
    return comp->getEdgeToBegin() * getCSTrafoLab2Local(comp);
}

inline CoordinateSystemTrafo OpalBeamline::getNominalExitTransform(
        const std::shared_ptr<Component>& comp) const {
    return comp->getEdgeToEnd() * getCSTrafoLab2Local(comp);
}

#endif  // OPAL_BEAMLINE_H
