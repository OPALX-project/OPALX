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

#include <cmath>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <string>

#include "AbsBeamline/BendBase.h"
#include "AbsBeamline/Component.h"
#include "PartBunch/PartBunch.h"

#include "AbsBeamline/Marker.h"
#include "BeamlineGeometry/PlacedElement.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
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

    /**
     * @brief Return the rigid runtime field-frame transform used for tracking queries.
     *
     * The placement bridge distinguishes the nominal body frame used for
     * geometry/export from the local field chart used by runtime field
     * evaluation. Most straight elements use the body frame directly. Analytic
     * bends use this rigid transform only as an auxiliary frame; the actual
     * longitudinal coordinate used for runtime field queries may be
     * curvilinear, as implemented in `transformToFieldLocalCS()`.
     */
    CoordinateSystemTrafo getFieldCSTrafoLab2Local(const std::shared_ptr<Component>& comp) const;

    /**
     * @brief Transform a lab-space position into the component's runtime field chart.
     *
     * This first applies the rigid field-frame transform returned by
     * `getFieldCSTrafoLab2Local()` and then lets the component convert that rigid-frame point into
     * its final field chart.
     */
    Vector_t<double, 3> transformToFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;

    /**
     * @brief Transform a rigid field-frame point back to lab coordinates.
     *
     * This helper only inverts the rigid field-frame transform. For curvilinear components the
     * caller must provide a point already expressed in the rigid frame, not in the final field
     * chart.
     */
    Vector_t<double, 3> transformFromFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;

    /**
     * @brief Rotate a lab-space vector into the component field basis using the default chart
     * origin.
     *
     * This convenience overload preserves the historical rigid behavior for callers that do not
     * have a field-local position available.
     */
    Vector_t<double, 3> rotateToFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;

    /**
     * @brief Rotate a lab-space vector into the component field basis at a given field-local point.
     *
     * The supplied `fieldLocalPosition` is used by curvilinear elements to evaluate the local
     * tangent basis.
     */
    Vector_t<double, 3> rotateToFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& fieldLocalPosition,
            const Vector_t<double, 3>& r) const;

    /**
     * @brief Rotate a field-basis vector back to lab space using the default chart origin.
     *
     * This convenience overload preserves the historical rigid behavior for callers that do not
     * have a field-local position available.
     */
    Vector_t<double, 3> rotateFromFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const;

    /**
     * @brief Rotate a field-basis vector back to lab space at a given field-local point.
     *
     * The supplied `fieldLocalPosition` is used by curvilinear elements to evaluate the inverse
     * tangent-basis rotation.
     */
    Vector_t<double, 3> rotateFromFieldLocalCS(
            const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& fieldLocalPosition,
            const Vector_t<double, 3>& r) const;

    /**
     * @brief Return the placed-element view used by the bridge stage.
     *
     * This is the beamline-facing access point to the geometric placement
     * record. It keeps nominal placement, local correction, and port geometry
     * together while the legacy runtime still stores them on ElementBase.
     */
    PlacedElement getPlacedElement(const std::shared_ptr<Component>& comp) const;

    /**
     * @brief Return the nominal rigid placement transform \f$T_i\f$.
     *
     * In the language of the placement note, this is the nominal rigid
     * placement transform from the parent/reference base into the element body
     * frame. During the bridge stage this remains identical to the legacy
     * lab-to-local transform storage.
     */
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
    bool reportPortContinuityDiagnostics(std::ostream& out) const;
    void print(Inform&) const;
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    FieldList getElementByType(ElementType);

    void swap(OpalBeamline& rhs);
    void merge(OpalBeamline& rhs);

private:
    using PlacementAssembly = std::map<const ElementBase*, PlacedElement>;

    /**
     * @brief Update the nominal rigid placement transform of one element.
     *
     * This is the write-side bridge from the legacy beamline assembly code to
     * the new placement vocabulary. The stored quantity is the nominal rigid
     * placement transform \f$T_i\f$; local survey or misalignment corrections
     * remain on the element and are not composed here.
     */
    void setNominalPlacement(
            const std::shared_ptr<ElementBase>& element, const CoordinateSystemTrafo& parentToBody);

    /**
     * @brief Refresh the beamline-owned placed-element assembly record.
     *
     * Beamline assembly owns a snapshot of the geometric placement records used
     * by placement/export queries. This keeps the assembled placement model
     * distinct from the legacy storage that still lives on ElementBase during
     * the bridge stage.
     */
    void storePlacedElement(const std::shared_ptr<ElementBase>& element);

    /**
     * @brief Compile legacy reference-order placement into explicit nominal poses.
     *
     * This setup-stage bridge consumes the compatibility placement encoded via
     * `ELEMEDGE` and converts it once into nominal rigid placement transforms
     * \f$T_i\f$ stored on the elements. Existing callers may still invoke
     * `compute3DLattice()`, but new code should treat this as a setup
     * conversion, not as a runtime geometry query.
     */
    void compileCompatibilityPlacement();

    /**
     * @brief Return true when an unplaced element still requires `ELEMEDGE`.
     *
     * Analytic `SBEND` and `RBEND` elements need either an explicit reference
     * path coordinate from `ELEMEDGE` or an explicit body pose from `Z`. Without
     * either one, the compatibility bridge has no longitudinal start coordinate
     * for the bend field support.
     */
    static bool requiresElementEdgeWithoutExplicitPosition(const ElementBase& element);

    /**
     * @brief Resolve the longitudinal field interval start used during visit().
     *
     * If `ELEMEDGE` is set, this returns that value. Otherwise, explicitly
     * positioned elements, including analytic bends, use the nominal body-origin
     * coordinate \f$Z_\mathrm{body}\f$ as the start of their active interval.
     * Unplaced `SBEND` and `RBEND` elements remain strict and throw without
     * `ELEMEDGE`.
     */
    double resolveVisitStartField(const std::shared_ptr<ElementBase>& element) const;

    FieldList elements_m;
    std::vector<std::shared_ptr<Component>> declaredOrder_m;
    PlacementAssembly placementAssembly_m;
    bool prepared_m;
    bool compatibilityPlacementCompiled_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

template <class T>
inline void OpalBeamline::visit(const T& element, BeamlineVisitor&, PartBunch_t& bunch) {
    Inform msg("OPAL ");
    double startField = 0.0;
    double endField   = 0.0;
    std::shared_ptr<T> elptr(dynamic_cast<T*>(element.clone()));

    positionElementRelative(elptr);

    startField = resolveVisitStartField(elptr);

    elptr->initialise(&bunch, startField, endField);
    elements_m.push_back(BeamlineFieldElement(elptr, startField, endField));
    declaredOrder_m.push_back(elptr);
    placementAssembly_m.insert_or_assign(elptr.get(), elptr->getPlacedElement());
    prepared_m                       = false;
    compatibilityPlacementCompiled_m = false;
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

inline Vector_t<double, 3> OpalBeamline::transformToLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getPlacedElement(comp).getNominalBodyTransform().transformTo(r);
}

inline Vector_t<double, 3> OpalBeamline::transformFromLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getPlacedElement(comp).getNominalBodyTransform().transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getPlacedElement(comp).getNominalBodyTransform().rotateTo(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getPlacedElement(comp).getNominalBodyTransform().rotateFrom(r);
}

inline CoordinateSystemTrafo OpalBeamline::getFieldCSTrafoLab2Local(
        const std::shared_ptr<Component>& comp) const {
    return comp->getFieldCSTrafoLab2Local(getPlacedElement(comp));
}

inline Vector_t<double, 3> OpalBeamline::transformToFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return comp->transformFieldFrameToLocal(getFieldCSTrafoLab2Local(comp).transformTo(r));
}

inline Vector_t<double, 3> OpalBeamline::transformFromFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return getFieldCSTrafoLab2Local(comp).transformFrom(r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return rotateToFieldLocalCS(comp, Vector_t<double, 3>(0.0), r);
}

inline Vector_t<double, 3> OpalBeamline::rotateToFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& fieldLocalPosition,
        const Vector_t<double, 3>& r) const {
    return comp->rotateFieldFrameToLocal(
            getFieldCSTrafoLab2Local(comp).rotateTo(r), fieldLocalPosition);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& r) const {
    return rotateFromFieldLocalCS(comp, Vector_t<double, 3>(0.0), r);
}

inline Vector_t<double, 3> OpalBeamline::rotateFromFieldLocalCS(
        const std::shared_ptr<Component>& comp, const Vector_t<double, 3>& fieldLocalPosition,
        const Vector_t<double, 3>& r) const {
    return getFieldCSTrafoLab2Local(comp).rotateFrom(
            comp->rotateFieldLocalToFieldFrame(r, fieldLocalPosition));
}

inline PlacedElement OpalBeamline::getPlacedElement(const std::shared_ptr<Component>& comp) const {
    const auto found = placementAssembly_m.find(comp.get());
    if (found != placementAssembly_m.end()) {
        return found->second;
    }
    return comp->getPlacedElement();
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local(
        const std::shared_ptr<Component>& comp) const {
    return getPlacedElement(comp).getNominalBodyTransform();
}

inline CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline CoordinateSystemTrafo OpalBeamline::getMisalignment(
        const std::shared_ptr<Component>& comp) const {
    return getPlacedElement(comp).getMisalignment().getNominalToActual();
}

inline CoordinateSystemTrafo OpalBeamline::getNominalEntryTransform(
        const std::shared_ptr<Component>& comp) const {
    return getPlacedElement(comp).getNominalEntryTransform();
}

inline CoordinateSystemTrafo OpalBeamline::getNominalExitTransform(
        const std::shared_ptr<Component>& comp) const {
    return getPlacedElement(comp).getNominalExitTransform();
}

#endif  // OPAL_BEAMLINE_H
