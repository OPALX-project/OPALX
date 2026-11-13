//
// Class OpalBeamline
//   Runtime tracking elements, placement and field selection.
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
#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "PartBunch/PartBunch.h"

#include "AbsBeamline/Marker.h"
#include "Beamlines/Beamline.h"
#include "Utilities/ElementList.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/GeneralElementMapping.h"

#include "OPALTypes.h"

class BoundaryGeometry;
class MeshGenerator;

/**
 * @brief Placed runtime lattice used for external-field queries.
 *
 * visit() clones concrete sequence occurrences into a flattened element list. Element lookups
 * by position go through one GeneralElementMapping built from the placed elements.
 */
class OpalBeamline {
public:
    OpalBeamline();
    OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation);
    ~OpalBeamline();

    void activateElements();
    /// Select runtime elements whose isInside() support contains lab position x [m].
    /// Includes aperture checks; selection does not imply a nonzero field. Candidates come
    /// from the element grid once placement is final, otherwise every element is tested.
    ElementSet getElements(const Vector_t<double, 3>& x);

    /**
     * @brief Build the element grid from the placed elements with the current OPTION values.
     *
     * Call once placement, misalignment and field-map extents are final. A later change of the
     * element list discards the grid.
     */
    void buildElementMapping();

    /// The element grid, built on first use.
    const GeneralElementMapping& getElementMapping();

    /**
     * Get all elements in the beamline, regardless of their position.
     * @return Set of shared pointers to all runtime elements, not a reference-path ordering.
     */
    ElementSet getElements();

    Vector_t<double, 3> transformTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> transformFrom(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateFrom(const Vector_t<double, 3>& r) const;

    /// Lab <-> element frame, misalignment included: the frame in which the tracker applies
    /// the element to the bunch.
    Vector_t<double, 3> transformToLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> transformFromLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateToLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateFromLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const;

    /// Return the element's nominal body transform (global→local).
    CoordinateSystemTrafo getCSTrafoLab2Local(const std::shared_ptr<ElementBase>& comp) const;
    CoordinateSystemTrafo getCSTrafoLab2Local() const;
    CoordinateSystemTrafo getMisalignment(const std::shared_ptr<ElementBase>& comp) const;
    CoordinateSystemTrafo getNominalEntryTransform(const std::shared_ptr<ElementBase>& comp) const;
    CoordinateSystemTrafo getNominalExitTransform(const std::shared_ptr<ElementBase>& comp) const;

    void switchElementsOff();

    BoundaryGeometry* getBoundaryGeometry(const unsigned int&);

    unsigned long getFieldAt(
            const unsigned int&, const Vector_t<double, 3>&, const long&, const double&,
            Vector_t<double, 3>&, Vector_t<double, 3>&);
    unsigned long getFieldAt(
            const Vector_t<double, 3>&, const Vector_t<double, 3>&, const double&,
            Vector_t<double, 3>&, Vector_t<double, 3>&);

    template <class T>
    void visit(const T&, BeamlineVisitor&, PartBunch_t&);

    /**
     * @brief Add a copy of @p element to the tracked elements and call its initialise().
     * @note The copy gets the next line index (see ElementBase::getLineIndex()).
     */
    void addElement(const ElementBase& element, PartBunch_t& bunch);

    void prepareSections();
    void compute3DLattice();
    void save3DLattice();
    /// Append the physical lattice to a display mesh containing global diagnostics.
    /// Diagnostic geometry remains outside elements_m and field/placement selection.
    void save3DLattice(MeshGenerator& mesh);
    void save3DInput();
    void print(Inform&) const;
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    ElementList getElementByType(ElementType);

    void swap(OpalBeamline& rhs);
    void merge(OpalBeamline& rhs);

private:
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

    ElementList elements_m;
    bool prepared_m;
    bool referencePathPlacementCompiled_m;
    /// Element grid for lookups by position; empty until built.
    std::shared_ptr<const GeneralElementMapping> elementMapping_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

inline void OpalBeamline::addElement(const ElementBase& element, PartBunch_t& bunch) {
    std::shared_ptr<ElementBase> elptr(element.clone());

    // Placement (both 6D-pose and ELEMEDGE) is resolved later, in one pass, by
    // placeElementsAlongReferencePath() during prepareSections().
    elptr->initialise(&bunch);
    elptr->setLineIndex(elements_m.size() + 1);
    elements_m.push_back(elptr);
    prepared_m                       = false;
    referencePathPlacementCompiled_m = false;
    elementMapping_m.reset();
}

template <class T>
inline void OpalBeamline::visit(const T& element, BeamlineVisitor&, PartBunch_t& bunch) {
    addElement(element, bunch);
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
        const std::shared_ptr<ElementBase>& comp) const {
    return comp->getCSTrafoGlobal2Local();
}

inline Vector_t<double, 3> OpalBeamline::transformToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    // Applied one after the other: composing the two frames costs far more per call.
    return comp->getMisalignment().transformTo(comp->getCSTrafoGlobal2Local().transformTo(r));
}

inline Vector_t<double, 3> OpalBeamline::transformFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return comp->getCSTrafoGlobal2Local().transformFrom(comp->getMisalignment().transformFrom(r));
}

inline Vector_t<double, 3> OpalBeamline::rotateToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return comp->getMisalignment().rotateTo(comp->getCSTrafoGlobal2Local().rotateTo(r));
}

inline Vector_t<double, 3> OpalBeamline::rotateFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return comp->getCSTrafoGlobal2Local().rotateFrom(comp->getMisalignment().rotateFrom(r));
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline CoordinateSystemTrafo OpalBeamline::getMisalignment(
        const std::shared_ptr<ElementBase>& comp) const {
    return comp->getMisalignment();
}

inline CoordinateSystemTrafo OpalBeamline::getNominalEntryTransform(
        const std::shared_ptr<ElementBase>& comp) const {
    return comp->getGeometry().getEdgeToBegin() * getCSTrafoLab2Local(comp);
}

inline CoordinateSystemTrafo OpalBeamline::getNominalExitTransform(
        const std::shared_ptr<ElementBase>& comp) const {
    return comp->getGeometry().getEdgeToEnd() * getCSTrafoLab2Local(comp);
}

#endif  // OPAL_BEAMLINE_H
