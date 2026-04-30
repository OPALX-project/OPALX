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

#include <map>
#include <memory>
#include <set>
#include <string>

#include "AbsBeamline/Component.h"
#include "PartBunch/PartBunch.h"

#include "AbsBeamline/Marker.h"
#include "BeamlineGeometry/PlacedElement.h"
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

    FieldList elements_m;
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

    if (elptr->isElementPositionSet()) startField = elptr->getElementPosition();

    elptr->initialise(&bunch, startField, endField);
    elements_m.push_back(BeamlineFieldElement(elptr, startField, endField));
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
