//
// Class FieldmapElement
//   A beamline element whose only field source is a tabulated field map.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_FieldmapElement_HH
#define OPALX_FieldmapElement_HH

#include "AbsBeamline/ElementBase.h"

#include <string>

class Fieldmap;

/**
 * @class FieldmapElement
 * @brief A straight element whose field comes entirely from a magnetostatic field map.
 *
 * Unlike SBEND, RBEND or QUADRUPOLE, this element has no analytic field of its own and
 * derives nothing from the input file: the map supplies the field, and the map's own extent
 * supplies the element's box. Nothing about the element's orientation is computed from a
 * bend angle, which is what makes it usable for a measured magnet whose real deflection is
 * whatever the tabulated field produces.
 *
 * @note The element's local frame IS the map's coordinate frame. The map is queried at the
 *       particle position without any shift, so the 6D lab pose (X, Y, Z, THETA, PHI, PSI)
 *       positions the map's own origin.
 *       This is forced by Fieldmap::applyField(), which reads the particle positions
 *       directly and takes no offset.
 *
 * @note Because of that, a map whose z range does not start at zero has its field window
 *       offset from the body window [0, L] that ElementBase assumes. markOutsideAperture()
 *       and getBoundingBoxInLabCoords() are overridden to use the field window instead.
 *
 * @note Placement is restricted to the 6D lab pose. ELEMEDGE positions an element by path
 *       length along the reference orbit, which needs the element's length along that orbit;
 *       this element's extent comes from its map and is not a path length. The rejection
 *       lives in OpalFieldmapElement::update().
 */
class FieldmapElement : public ElementBase {
public:
    /* ============================== Constructors ============================== */
    explicit FieldmapElement(const std::string& name);
    FieldmapElement();
    FieldmapElement(const FieldmapElement&);
    virtual ~FieldmapElement();

    /* ============================== Apply functions =========================== */
    /**
     * @brief Add the map's field to every particle in the bunch.
     *
     * Delegates straight to Fieldmap::applyField(). Particles outside the map get no
     * contribution; the reader kernels already skip them.
     */
    virtual void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /**
     * @brief Add the map's field at a single position.
     *
     * @param R Position in the element's local frame, which is the map's frame.
     * @param P Momentum (unused).
     * @param t Time (unused, the map is static).
     * @param E Electric field (untouched, the map is magnetostatic).
     * @param B Magnetic field, accumulated into.
     */
    virtual void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /**
     * @brief Add the map's field for the reference particle.
     *
     * @return true only for a genuine aperture hit. Being outside the map's transverse
     *         extent is NOT a loss: the particle simply gets no field. Returning true there
     *         would make OrbitThreader treat a benign miss as hitting material.
     */
    virtual bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /**
     * @brief Mark particles inside the field window but outside the transverse aperture.
     *
     * Overridden because ElementBase gates on the body window [0, L], which is offset from
     * the field window whenever the map's z range does not start at zero.
     */
    virtual size_t markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /* ============================== Functions ================================= */
    /// @brief Apply visitor to FieldmapElement.
    virtual void accept(BeamlineVisitor&) const override;

    /**
     * @brief Load the map, check its type, and take the element's box from its extent.
     *
     * The map header is parsed by the reader's constructor inside Fieldmap::getFieldmap(),
     * so both extents are available here -- before goOnline() reads the data. This runs on
     * the tracked clone in OpalBeamline::visit(), which is before the element list is sorted
     * and before PlacementResolver runs, so the map-derived length reaches every consumer.
     */
    virtual void initialise(PartBunch_t* bunch) override;

    virtual void finalise() override;

    /// @brief Read the map data and go online.
    virtual void goOnline(const double& kineticEnergy) override;

    /// @brief Release the map data and go offline.
    virtual void goOffline() override;

    virtual ElementType getType() const override;

    /// @brief The longitudinal field-support interval, taken from the map.
    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    /**
     * @brief Lab-frame bounding box, built from the field window rather than [0, L].
     */
    virtual BoundingBox getBoundingBoxInLabCoords() const override;

    /**
     * @brief A finite transverse support envelope for placement and visualisation.
     *
     * Prefers a finite configured aperture, otherwise falls back to the map's own transverse
     * extent. A one-dimensional map has no transverse extent, so this returns false for one
     * unless an APERTURE was given.
     *
     * @param horizontalRadius Output half-width in x [m].
     * @param verticalRadius Output half-width in y [m].
     * @return true if a finite envelope is available.
     */
    bool getSupportEnvelope(double& horizontalRadius, double& verticalRadius) const;

    /* ============================== Accessors ================================= */
    /// @brief Set the field map filename.
    void setFieldMapFN(const std::string& fn);

    const std::string& getFieldMapFN() const;

    /**
     * @brief Set the plain multiplier applied to the tabulated magnetic field.
     *
     * This is not a normalised strength: the G4beamline readers store absolute Tesla and do
     * not normalise, so a scale of 1 reproduces the map as written. It is the same quantity
     * as the `current=` given where the map is placed in a G4beamline input.
     */
    void setScale(double scale);

    double getScale() const;

    /**
     * @brief Set the plain multiplier applied to the tabulated electric field.
     *
     * Separate from setScale() because G4beamline scales the two fields independently:
     * `current` and `normB` for the magnetic field, `gradient` and `normE` for the electric
     * one. This is the same quantity as the `gradient=` given where the map is placed. It
     * has no effect on a map that carries no electric field.
     */
    void setEScale(double escale);

    double getEScale() const;

    /**
     * @brief Read the map back to front: mirror it in z and negate the longitudinal
     *        component. Only G4beamline cylinder maps support this.
     */
    void setIsZReversed(bool zReverse);

    bool getIsZReversed() const;

    /// @brief True if the loaded map declares a transverse extent.
    bool hasTransverseExtent() const;

private:
    void operator=(const FieldmapElement&) = delete;

    /* ============================== Variables ================================= */

    /// Name of the field map file.
    std::string filename_m;

    /// The loaded map. Null until initialise(), and again after goOffline().
    Fieldmap* fieldmap_m;

    /// Plain multiplier on the tabulated magnetic field.
    double scale_m;

    /// Plain multiplier on the tabulated electric field.
    double escale_m;

    /// Load the map mirrored in z.
    bool isZReversed_m;

    /// Start of the field support in the local chart [m].
    double startField_m;

    /// End of the field support in the local chart [m].
    double endField_m;

    /// True if the map answered the transverse-extent query (2D and 3D maps do, 1D do not).
    bool hasTransverseExtent_m;

    /// Half-width of the map's box in x [m], valid only when hasTransverseExtent_m.
    double halfWidthX_m;

    /// Half-width of the map's box in y [m], valid only when hasTransverseExtent_m.
    double halfWidthY_m;
};

#endif  // OPALX_FieldmapElement_HH
