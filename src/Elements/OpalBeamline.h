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

#include <array>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "PartBunch/PartBunch.h"

#include "Utilities/ElementList.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/GeneralElementMapping.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include "OPALTypes.h"

/**
 * @brief Placed runtime lattice used for external-field queries.
 *
 * addElement() takes the runtime elements into a flat element list. Element lookups by position
 * go through one GeneralElementMapping built from the placed elements.
 */
class OpalBeamline {
public:
    OpalBeamline();
    OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation);

    /// Select runtime elements whose isInside() support contains lab position x [m] in their
    /// nominal placement (transformToLocalCS()), as the reference particle sees them.
    /// Includes aperture checks; selection does not imply a nonzero field. Candidates come
    /// from the element grid once placement is final, otherwise every element is tested.
    ElementSet getElements(const Vector_t<double, 3>& x);

    /**
     * @brief Build the element grid from the placed elements.
     *
     * Call once placement, misalignment and field-map extents are final. A later change of the
     * element list discards the grid.
     * @param grid Cells along x, y and z (OPTION ELEMENTMAPPINGGRID).
     * @param reach Reach [m] of the elements without REACH (OPTION ELEMENTMAPPINGREACH).
     */
    void buildElementMapping(const std::array<unsigned, 3>& grid, double reach);

    /// The element grid. @throws LogicalError if buildElementMapping() has not built it.
    const GeneralElementMapping& getElementMapping() const;

    /**
     * Get all elements in the beamline, regardless of their position.
     * @return Set of shared pointers to all runtime elements, not a reference-path ordering.
     */
    ElementSet getElements();

    /// Lab <-> element frame of the nominal placement, without misalignment: the frame in which
    /// the reference particle sees the element. The bunch sees the misaligned element in
    /// elementFrame().
    static Vector_t<double, 3> transformToLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r);
    static Vector_t<double, 3> transformFromLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r);
    static Vector_t<double, 3> rotateToLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r);
    static Vector_t<double, 3> rotateFromLocalCS(
            const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r);

    /// Return the element's nominal body transform (global→local).
    static CoordinateSystemTrafo getCSTrafoLab2Local(const std::shared_ptr<ElementBase>& comp);
    CoordinateSystemTrafo getCSTrafoLab2Local() const;
    static CoordinateSystemTrafo getMisalignment(const std::shared_ptr<ElementBase>& comp);
    CoordinateSystemTrafo getNominalEntryTransform(const std::shared_ptr<ElementBase>& comp) const;
    CoordinateSystemTrafo getNominalExitTransform(const std::shared_ptr<ElementBase>& comp) const;

    void switchElementsOff();

    /// Turn reference recording of every MONITOR on or off (Monitor::setRecordReference()).
    void setMonitorsRecordReference(bool record);

    /**
     * @brief Apply the field of every element near each active container whose fieldFrame() is
     *        ElementFrame::Local, with the container in that element's local frame. Collective.
     *
     * Elements get the bunch time, time step and step number (ElementBase::apply()). Each
     * container goes from one element's frame directly to the next and returns to its reference
     * frame after its last element, also when an element throws.
     */
    void applyExternalFields(PartBunch_t& bunch);

    /**
     * @brief Mark the particles outside the transverse aperture of each element near their
     *        container. The particles stay in their frame; the marks wait for the next deletion.
     * @return Global number of newly marked particles. Collective.
     */
    size_t applyElementApertures(PartBunch_t& bunch);

    /// Beamline (emission source) frame for the particles of @p pc.
    ParticleFrame beamlineFrame(const ParticleContainer_t& pc) const;

    /**
     * @brief Mark the particles of the active containers that move backwards at or behind the
     *        source plane, read in the beamline frame; the particles stay where they are.
     *
     * The source plane is the Dirichlet plane of the space-charge solver at z = R0Z of the
     * emitting EMISSIONSOURCE; particles are lost 5 cm upstream of it, at the entrance of the
     * legacy OPAL SOURCE element.
     * @return Global number of newly marked particles; 0 without a Dirichlet plane. Collective.
     */
    size_t markBackwardParticlesAtSourcePlane(
            PartBunch_t& bunch, const opalx::spacecharge::DirichletPlaneConfig& plane) const;

    /**
     * @brief Add the fields of @p elements at one lab position to E and B, in lab coordinates.
     *
     * The one host loop over elements: each element gets the position in its frame
     * (transformToLocalCS()), adds its field with ElementBase::applyOnHost(), and the field is
     * rotated back. Elements are summed in line order.
     * @param elements The elements to sum, e.g. getElements(r).
     * @param r Lab position [m].
     * @param t Time [s].
     * @param E Electric field [V/m], added to.
     * @param B Magnetic field [T], added to.
     * @return Whether r is in the loss region of one of the elements, by the test of
     *         ElementBase::markOutsideAperture().
     */
    static bool applyOnHost(
            const ElementSet& elements, const Vector_t<double, 3>& r, double t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    /// The same over getElements(r).
    bool applyOnHost(
            const Vector_t<double, 3>& r, double t, Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    /**
     * @brief Let every MONITOR in @p elements record the reference particle at lab position @p r
     *        with lab momentum @p p and time @p t (Monitor::recordReference()).
     */
    void recordReference(
            const ElementSet& elements, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p,
            double t);

    /**
     * @brief The external fields {B (T), E (V/m)} at the reference particle of each container of
     *        @p bunch, at the middle of the last time step; zero for an empty container.
     */
    std::vector<std::array<Vector_t<double, 3>, 2>> getReferenceFields(const PartBunch_t& bunch);

    /**
     * @brief Add @p element itself (no copy) to the tracked elements and call its initialise().
     * @note For elements built for one RUN only (opalx::run::buildElements()).
     */
    void addElement(std::shared_ptr<ElementBase> element, PartBunch_t& bunch);

    void prepareSections();
    void compute3DLattice();
    /**
     * @brief Write the element positions to data/<inputBasename>_ElementPositions.txt and the
     *        element meshes to data/<inputBasename>_ElementPositions.py (rank 0 only).
     * @param append Append to an existing positions file instead of replacing it.
     */
    void save3DLattice(const std::string& inputBasename, bool append);
    /// Write data/<inputBasename>_3D.opal: the input file @p inputFileName with every element
    /// placed by absolute coordinates (rank 0 only).
    void save3DInput(const std::string& inputFileName, const std::string& inputBasename);

    ElementList getElementByType(ElementType);

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

    /**
     * @brief Call @p func for each selected element intersecting each active container, with
     * the container in that element's local frame.
     *
     * The container goes from one selected element's frame directly to the next and returns to
     * its reference frame after its last element, so no frame outlives the call. Unselected
     * elements cost no transform.
     *
     * @param selects Whether @p func runs for an element.
     */
    void forEachElementInBunchFrame(
            PartBunch_t& bunch, const std::function<bool(const ElementBase&)>& selects,
            const std::function<
                    void(const std::shared_ptr<ElementBase>&,
                         const std::shared_ptr<ParticleContainer_t>&)>& func);

    /// Elements near container @p ci in the order of the element grid; empty for an inactive
    /// container or one without its own or its root's reference particle. Collective.
    std::vector<std::shared_ptr<ElementBase>> candidateElements(
            const PartBunch_t& bunch, size_t ci) const;

    /// Local frame of @p element, misalignment included, for the particles of @p pc.
    ParticleFrame elementFrame(
            const std::shared_ptr<ElementBase>& element, const ParticleContainer_t& pc) const;

    ElementList elements_m;
    bool prepared_m;
    bool referencePathPlacementCompiled_m;
    /// Element grid for lookups by position; empty until built.
    std::shared_ptr<const GeneralElementMapping> elementMapping_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

inline void OpalBeamline::addElement(std::shared_ptr<ElementBase> elptr, PartBunch_t& bunch) {
    // Placement (both 6D-pose and ELEMEDGE) is resolved later, in one pass, by
    // placeElementsAlongReferencePath() during prepareSections().
    elptr->initialise(&bunch);
    elptr->setLineIndex(elements_m.size() + 1);
    elements_m.push_back(std::move(elptr));
    prepared_m                       = false;
    referencePathPlacementCompiled_m = false;
    elementMapping_m.reset();
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local(
        const std::shared_ptr<ElementBase>& comp) {
    return comp->getCSTrafoGlobal2Local();
}

inline Vector_t<double, 3> OpalBeamline::transformToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) {
    return getCSTrafoLab2Local(comp).transformTo(r);
}

inline Vector_t<double, 3> OpalBeamline::transformFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) {
    return getCSTrafoLab2Local(comp).transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) {
    return getCSTrafoLab2Local(comp).rotateTo(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromLocalCS(
        const std::shared_ptr<ElementBase>& comp, const Vector_t<double, 3>& r) {
    return getCSTrafoLab2Local(comp).rotateFrom(r);
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline CoordinateSystemTrafo OpalBeamline::getMisalignment(
        const std::shared_ptr<ElementBase>& comp) {
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
