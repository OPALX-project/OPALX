//
// Class OpalBeamline
//   Runtime tracking elements, placement, field selection and transfer-map views.
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

#include "OPALTypes.h"

class BoundaryGeometry;
class MeshGenerator;

/**
 * @brief Pinned view of one representative nominal owner and its immutable map.
 *
 * Copies survive attachment-vector clearing/recalculation. Does not describe the
 * full owner set or ownerless intervals. Polynomial access ends at map-storage
 * release/Kokkos finalization; metadata and occurrence handles remain usable.
 */
struct ElementTransferMapRef {
    /// Lifetime handle to one representative occurrence.
    std::shared_ptr<const ElementBase> element;
    /// Shared immutable segment record, also retained by the complete pass.
    LinearTransferMapHandle map;
    std::size_t pass{0};
};

/**
 * @brief Placed runtime lattice used for external-field queries and map attachment.
 *
 * visit() clones concrete sequence occurrences into a flattened element list.
 * These clones retain logical BeamlineMembership, but start without calculated maps
 * or overlap flags. Nominal body ownership and field-support selection remain separate;
 * neither is inferred from that membership tag or from an element's name.
 */
class OpalBeamline {
public:
    OpalBeamline();
    OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation);
    /** Share runtime element occurrences, as merge() does; pass handles remain snapshots.
     * Copying a nonempty lattice marks both lists as potentially aliased. Setup
     * copies and mutations must be serialized, as for publication and placement.
     */
    OpalBeamline(const OpalBeamline&);
    /// Replace the list with shared occurrences and retain the conservative alias state.
    OpalBeamline& operator=(const OpalBeamline&);
    ~OpalBeamline();

    void activateElements();
    /// Select runtime elements whose isInside() support contains lab position x [m].
    /// Includes aperture checks; selection does not imply a nonzero field or body ownership.
    std::set<std::shared_ptr<ElementBase>> getElements(const Vector_t<double, 3>& x);

    /// Nominal body owners at lab position x [m], using ElementBase::isInsideBody().
    /// Used for map attachment; never use this set to restrict fringe-field evaluation.
    std::set<std::shared_ptr<ElementBase>> getBodyElements(const Vector_t<double, 3>& x) const;

    /**
     * Get all elements in the beamline, regardless of their position.
     * @return Set of shared pointers to all runtime elements, not a reference-path ordering.
     */
    std::set<std::shared_ptr<ElementBase>> getElements();

    Vector_t<double, 3> transformTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> transformFrom(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateTo(const Vector_t<double, 3>& r) const;
    Vector_t<double, 3> rotateFrom(const Vector_t<double, 3>& r) const;

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

    /**
     * @brief Ordered owner view; shared-owner segments appear once by scoped traversal ID.
     *
     * Sorted by entrance path length, then attachment ordinal, then owner name.
     * Ownerless intervals are absent. Use getTransferMapPass() for the complete
     * sequence or combined polynomial. Returned values pin maps across recalculation.
     */
    std::vector<ElementTransferMapRef> getLinearTransferMapsInReferenceOrder() const;

    /// Complete ordered snapshot, including ownerless segments. Pins TPS until teardown.
    std::shared_ptr<const TransferMapPass> getTransferMapPass() const { return transferMapPass_m; }
    /// Remove the local pass and attachments on its elements, including shared occurrences.
    /// Other aliases observe cleared attachments; retained immutable passes remain valid.
    void clearTransferMaps();
    /// Publish validated attachments atomically after staging allocations. Empty clears.
    /// Owners must belong to this lattice and segments to the supplied pass.
    void publishTransferMapPass(
            std::shared_ptr<const TransferMapPass> pass,
            const std::vector<std::pair<std::shared_ptr<ElementBase>, LinearTransferMapAttachment>>&
                    attachments);

    void swap(OpalBeamline& rhs);
    /// Append shared occurrences without consuming rhs; both lists may subsequently alias.
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

    std::shared_ptr<const TransferMapPass> transferMapPass_m;
    ElementList elements_m;
    // Copies/merges share attachment-bearing ElementBase objects. A local null pass
    // is then insufficient to skip invalidation. Retain this conservative flag
    // through clear/swap, even after other aliases die; no ownership scan is needed.
    mutable bool elementsMayBeShared_m = false;
    bool prepared_m;
    bool referencePathPlacementCompiled_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

template <class T>
inline void OpalBeamline::visit(const T& element, BeamlineVisitor&, PartBunch_t& bunch) {
    std::shared_ptr<T> elptr(dynamic_cast<T*>(element.clone()));

    // Placement (both 6D-pose and ELEMEDGE) is resolved later, in one pass, by
    // placeElementsAlongReferencePath() during prepareSections().
    elptr->initialise(&bunch);
    elements_m.push_back(elptr);
    // A different lattice can publish attachments onto shared occurrences.
    // Fresh unmapped lattices still avoid any all-element invalidation scan.
    if (transferMapPass_m || elementsMayBeShared_m) clearTransferMaps();
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
        const std::shared_ptr<ElementBase>& comp) const {
    return comp->getCSTrafoGlobal2Local();
}

inline Vector_t<double, 3> OpalBeamline::transformToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).transformTo(r);
}

inline Vector_t<double, 3> OpalBeamline::transformFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).rotateTo(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) const {
    return getCSTrafoLab2Local(comp).rotateFrom(r);
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
