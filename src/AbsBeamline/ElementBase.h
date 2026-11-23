//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPALX_ElementBase_HH
#define OPALX_ElementBase_HH

#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/ElementRegions.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"
#include "BeamlineGeometry/Geometry.h"
#include "OPALTypes.h"
#include "Structure/BoundingBox.h"
#include "Utilities/GeneralOpalException.h"
#include "VectorMath.h"

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

class BeamlineVisitor;
class Channel;
class ConstChannel;

namespace opalx::run {
    struct OptionValues;
    struct RunContext;
}  // namespace opalx::run

enum class ElementType : unsigned short {
    ANY,
    BEAMLINE,
    BOX,
    COLLIMATOR,
    DRIFT,
    FIELDMAP,
    LASER,
    MARKER,
    MONITOR,
    MULTIPOLE,
    MULTIPOLET,
    RFCAVITY,
    TRAVELINGWAVE,
    SBEND,
    RBEND,
    RBEND3D,
    VACUUM,
    SOLENOID,
    SOURCE,
    CONSTANTEFIELDCAVITY,
    CONSTANTFOCUSING
};

/// Oriented box: centre, orthonormal axes and half-extents along those axes [m].
struct SupportBox {
    Vector_t<double, 3> centre;
    std::array<Vector_t<double, 3>, 3> axes;
    Vector_t<double, 3> half;

    /// Whether @p point lies in the box, enlarged by @p tolerance [m] on every side.
    bool contains(const Vector_t<double, 3>& point, double tolerance = 0.0) const {
        const Vector_t<double, 3> d = point - centre;
        for (unsigned i = 0; i < 3; ++i) {
            if (std::abs(dot(d, axes[i])) > half(i) + tolerance) return false;
        }
        return true;
    }
};

/// Frame in which an element callback needs the particle data.
enum class ElementFrame : unsigned char {
    None,  ///< The callback leaves the particles untouched; the tracker neither calls nor moves.
    Local  ///< Element-local coordinates, misalignment included.
};

/**
 * @class ElementBase
 * @brief Base class for all beam line elements.
 *
 * An element owns its identity (name, type), its placement (global-to-local
 * coordinate transform, ELEMEDGE position, misalignment), its transverse
 * aperture, and its field model (the apply() family, called by the tracker
 * every time step, and getFieldExtent(), the field-support interval used for
 * element selection). All geometric state — lengths, bend angle, pole-face
 * rotations, edge transforms — lives in the element's Geometry, accessed
 * through getGeometry().
 *
 * The parser builds lattices by cloning registered prototypes: clone() makes
 * a deep copy, copyStructure() copies a structure while reusing elements
 * marked sharable (makeSharable()). Algorithms dispatch on the concrete
 * element type through accept(BeamlineVisitor&).
 */
/**
 * @brief Time of the current tracking step, passed to ElementBase::apply().
 *
 * Elements used to read these values from the bunch (RefPartBunch_m).
 */
struct StepTime {
    double t       = 0.0;  ///< Bunch time at the start of the step (s).
    double dt      = 0.0;  ///< Time step (s).
    long long step = 0;    ///< Global step number.
};

/// How ParallelTracker uses an element of the tracked line.
enum class TrackedAs {
    FieldElement,  ///< Copied into OpalBeamline: applies its field or aperture, or records.
    Ignored,       ///< Left out of tracking (MARKER).
    Unsupported    ///< Tracking it is an error (LASER).
};

class ElementBase : public std::enable_shared_from_this<ElementBase> {
public:
    /* ========================= Construction & lifecycle ====================== */

    /// @brief Constructor with given name.
    /// @param name The element name.
    explicit ElementBase(const std::string& name);

    /// @brief Default constructor.
    ElementBase();

    /// @brief Copy constructor.
    ElementBase(const ElementBase&);

    /// @brief Destructor.
    virtual ~ElementBase();

    /// @brief Return an identical deep copy of the element.
    /// @return The cloned element.
    virtual ElementBase* clone() const = 0;

    /// @brief Make a structural copy.
    /// @note A fresh copy of the beam line structure is made, but sharable
    ///       elements remain shared.
    /// @return The copied structure.
    virtual ElementBase* copyStructure();

    /// @brief Test if the element can be shared.
    /// @return True if the element is sharable.
    bool isSharable() const;

    /// @brief Set the sharable flag.
    /// @note The whole structure depending on this is marked as sharable.
    ///       After this call a copyStructure() call reuses the element.
    virtual void makeSharable();

    /* =============================== Identity ================================ */

    /// @brief Get the element name.
    /// @return The element name.
    virtual const std::string& getName() const;

    /// @brief Set the element name.
    /// @param name The element name.
    virtual void setName(const std::string& name);

    /// @brief Get the element type.
    /// @note Default returns ElementType::ANY; concrete elements override.
    /// @return The element type.
    virtual ElementType getType() const;

    /// @brief Get the element type as a string.
    /// @return The element type string.
    std::string getTypeString() const;

    /// @brief Get the string for a given element type.
    /// @param type The element type.
    /// @return The element type string.
    static std::string getTypeString(ElementType type);

    /// @brief Apply a visitor.
    /// @note This method must be overridden by derived classes. It should call
    ///       the method of the visitor corresponding to the element class. If
    ///       any error occurs, this method throws an exception.
    /// @param visitor The visitor to apply.
    virtual void accept(BeamlineVisitor& visitor) const = 0;

    /* =============================== Geometry ================================ */

    /// @brief Get the geometry. Version for non-constant object.
    /// @note Supplied by the representation layer.
    /// @return The element geometry.
    virtual Geometry& getGeometry() = 0;

    /// @brief Get the geometry. Version for constant object.
    /// @return The element geometry.
    virtual const Geometry& getGeometry() const = 0;

    /// @brief Whether the element's field applies at r: getRegions().inFieldRegion(r).
    /// @param r Position [m] in the element's local frame.
    virtual bool isInside(const Vector_t<double, 3>& r) const;

    /**
     * @brief Transverse range of the element's own field gate, if limited.
     *
     * MULTIPOLET returns its HAPERT/VAPERT gate and VARIABLE_RF_CAVITY its WIDTH/HEIGHT. Without
     * an APERTURE this range is the aperture, and particles outside it within the reach are lost.
     * @return Range in the local frame (arc coordinates for a sector bend); empty if unlimited.
     */
    virtual std::optional<ApertureShape> getFieldRange() const;

    /**
     * @brief Transverse range of the element's field map, if it has one with a transverse extent.
     *
     * Without an APERTURE this range is the element's body across and its aperture, so the field
     * acts inside the map and nothing is deleted. The body along z is the field window.
     * @return Range in the local frame; empty without a field map or for an on-axis (1D) map.
     */
    virtual std::optional<ApertureShape> getFieldMapRange() const;

    /**
     * @brief Where the element applies its field and where it deletes particles.
     *
     * Field window: getFieldExtent(). Body window: [startZ, startZ + L), from 0 in the arc
     * coordinates of a sector bend. Reach: getReach(), else defaultElementReach, across both
     * axes. buildTracking() sets OPTION ELEMENTMAPPINGREACH on every tracked element without
     * REACH.
     * Aperture, in this order:
     * - the APERTURE (values of 1e5 m or more count as none); particles outside it are lost;
     * - getFieldMapRange(); nothing is lost;
     * - getFieldRange(); particles outside it are lost;
     * - the reach box; nothing is lost.
     */
    ElementRegions getRegions() const;

    /**
     * @brief Check that an aperture that deletes particles lies inside the reach.
     * @throws OpalException naming the element if the APERTURE, or the field range used in its
     *         place, is wider than the reach.
     */
    void checkRegions() const;

    /// @brief Get the bounding box.
    /// @return The bounding box in lab coordinates.
    virtual BoundingBox getBoundingBoxInLabCoords() const;

    /**
     * @brief Local bounding box of the region where this element is taken to act.
     *
     * Each transverse half-width is max(aperture, @p reach), where an aperture of 1e5 m or
     * more counts as none, and z joins getFieldExtent() with the body [startZ, startZ + L].
     * A field map without APERTURE is its map instead: getFieldMapRange() across and
     * getFieldExtent() along z. Overrides add their own support where it differs (bend orbits,
     * cavity widths, annular sectors, the monitor plane). The box contains getRegions(): the
     * field region and the loss region. A sector-bend geometry uses the arc chart of
     * GeometryHelper::toBendArcCoords(). Valid after placement and initialise().
     * @param reach Minimum transverse half-width [m]: the element's REACH or the global default.
     * @return Box [m] in the local frame used by isInside() and apply().
     */
    virtual BoundingBox getLocalSupportBox(double reach) const;

    /**
     * @brief The same region as getLocalSupportBox() as one or more oriented boxes.
     *
     * GeneralElementMapping activates the element for a bunch whose box intersects one of
     * them. A sector-bend arc is split into chord-aligned segments, so a straight line along
     * the entrance tangent stays outside; every other element returns its bounding box.
     * @param reach Minimum transverse half-width [m], as for getLocalSupportBox().
     * @return Boxes [m] in the local frame used by isInside() and apply().
     */
    virtual std::vector<SupportBox> getLocalSupportBoxes(double reach) const;

    /// @brief Set the reach [m]: the element's REACH, or OPTION ELEMENTMAPPINGREACH.
    void setReach(double reach);

    /// @brief The reach [m]; empty if setReach() was not called.
    std::optional<double> getReach() const;

    /* ===================== Coordinate system & placement ===================== */

    /// @brief Set the lab -> element entrance frame trafo.
    /// @param ori The trafo.
    void setCSTrafoGlobal2Local(const CoordinateSystemTrafo& ori);

    /// @brief Get the lab->element entrance frame trafo.
    /// @return The trafo.
    CoordinateSystemTrafo getCSTrafoGlobal2Local() const;

    /// @brief Set the misaligment.
    /// @param cst The misalignment.
    void setMisalignment(const CoordinateSystemTrafo& cst);

    /// Not implemented.
    void getMisalignment(double& x, double& y, double& s) const;

    /// @brief Get the misalignment.
    /// @return The misalignment.
    CoordinateSystemTrafo getMisalignment() const;

    /// @brief Unlock the position so the global -> local transform can change.
    void releasePosition();

    /// @brief Lock the position so the global -> local transform cannot change.
    void fixPosition();

    /// @brief Test if the position is locked.
    /// @return True if the position is fixed.
    bool isPositioned() const;

    /// @brief Set the rotation about the z axis in the bend frame.
    /// @param rotation The rotation angle.
    void setRotationAboutZ(double rotation);

    /// @brief Get the rotation about the z axis in the bend frame.
    /// @return The rotation angle.
    double getRotationAboutZ() const;

    /// @brief Set the ELEMEDGE position of the element.
    /// @param elemedge The element edge position.
    void setElementPosition(double elemedge);

    /// @brief Get the ELEMEDGE position of the element.
    /// @return The element edge position.
    double getElementPosition() const;

    /// @brief Test if the ELEMEDGE position has been set.
    /// @return True if ELEMEDGE has been set.
    bool isElementPositionSet() const;

    /**
     * @brief Path length this element adds to the reference path when placing the elements
     *        after it by ELEMEDGE. Default: the arc length of its geometry.
     */
    virtual double getPlacementPathLength() const { return getGeometry().getArcLength(); }

    /// How ParallelTracker uses this element; most elements are field elements.
    virtual TrackedAs trackedAs() const { return TrackedAs::FieldElement; }

    /* =============================== Aperture ================================ */

    void setAperture(const ApertureType& type, const std::vector<double>& args);
    std::pair<ApertureType, std::vector<double>> getAperture() const;

    /* ===================== Field application & physics ======================= */
    /**
     * Apply functions apply the element's electromagnetic field to the
     * particles. They are called inside OpalBeamline::applyExternalFields().
     */

    /**
     * @brief Apply to all particles. Kernel launch moved inside the function.
     * @param pc   The particles, in the frame fieldFrame() requests.
     * @param time Bunch time, time step and step number of the current step.
     */
    virtual void apply(const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& time);

    /**
     * @brief Frame that apply(pc) needs; the tracker moves the particles there before the call.
     *
     * Only FramelessElement returns ElementFrame::None; the tracker then skips both the call and
     * the transform. The default keeps a new element correct.
     */
    virtual ElementFrame fieldFrame() const { return ElementFrame::Local; }

    /**
     * @brief Apply to particle with position R and momentum P
     *
     * Fields are only applied to particles inside the element's field bounds
     * and transverse aperture.
     *
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     */
    virtual void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    /**
     * @brief Apply to reference particle with position R and momemtum P
     *
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    /**
     * @brief Mark local particles inside the element's longitudinal window but
     * outside its transverse aperture in pc->InvalidMask.
     *
     * Called by OpalBeamline::applyElementApertures. The kernel reads the element-local
     * position through @p toLocal; the stored particles do not change frame. ORs into
     * InvalidMask and never clears existing marks; the tracker's end-of-step deletion flushes
     * all marks. A particle is marked where getRegions().inLossRegion() holds: inside the body
     * window and the reach, outside the aperture.
     * No-op when getFlagDeleteOnTransverseExit() is false.
     *
     * @param toLocal Transform from the particles' current frame to the element-local frame.
     * @return local number of newly marked particles
     */
    virtual size_t markOutsideAperture(
            const std::shared_ptr<ParticleContainer_t>& pc, const RigidTransform& toLocal);

    // Design energy for elements such as RF-cavities
    virtual double getDesignEnergy() const;
    virtual void setDesignEnergy(const double& energy, bool changeable = true);

    /**
     * @brief Take the OPTION values and the run context of the RUN that tracks this element.
     * @note buildBeamline() calls it for every element before initialise(). Elements that read
     *       run values in initialise() or later override it; the default ignores them.
     */
    virtual void setRunValues(
            const opalx::run::OptionValues& /*options*/,
            const opalx::run::RunContext& /*context*/) {}

    // Setup
    virtual void initialise(PartBunch_t* bunch) = 0;

    // Clean-up
    virtual void finalise() = 0;

    /// Release runtime resources / flush element output (e.g. loss data); clears online_m.
    /// initialise() prepares them (field maps, output files) and sets online_m where it is read.
    virtual void goOffline();

    /**
     * @brief Return the field-support extent of the element.
     *
     * This is the longitudinal interval
     * \f$[z_\mathrm{field}^{\mathrm{begin}}, z_\mathrm{field}^{\mathrm{end}}]\f$
     * on which the external field model is evaluated in the element-local
     * chart. It may extend past the body interval [0, L] of the geometry, for
     * example when fringe fields extend beyond the hardware body or when a
     * field map occupies only part of the body.
     */
    virtual void getFieldExtent(double& zBegin, double& zEnd) const = 0;

    virtual int getRequiredNumberOfTimeSteps() const;

    /* ========================= User-defined attributes ======================= */

    /// Get attribute value.
    //  If the attribute does not exist, return zero.
    virtual double getAttribute(const std::string& aKey) const;

    /// Test for existence of an attribute.
    //  If the attribute exists, return true, otherwise false.
    virtual bool hasAttribute(const std::string& aKey) const;

    /// Remove an existing attribute.
    virtual void removeAttribute(const std::string& aKey);

    /// Set value of an attribute.
    virtual void setAttribute(const std::string& aKey, double val);

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel* getChannel(const std::string& aKey, bool create = false);

    /// Construct a read-only channel.
    //  This method constructs a Channel permitting read-only access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual const ConstChannel* getConstChannel(const std::string& aKey) const;

    /// Update element.
    //  This method stores all attributes contained in the AttributeSet to
    //  "*this".  The return value [b]true[/b] indicates success.
    bool update(const AttributeSet&);

    /* ============================= Miscellaneous ============================= */

    /// Set output filename
    void setOutputFN(std::string fn);
    /// Get output filename
    std::string getOutputFN() const;
    /// The output file name, without extension, of element @p name with OUTFN @p outputFN.
    static std::string outputFileName(const std::string& name, const std::string& outputFN);

    void setFlagDeleteOnTransverseExit(bool = true);
    bool getFlagDeleteOnTransverseExit() const;

protected:
    /// Transverse half-widths max(aperture, reach) [m]; an aperture of 1e5 m or more counts as
    /// none.
    std::pair<double, double> getSupportHalfWidths(double reach) const;

    /// Whether an APERTURE was given: a half-width below 1e5 m.
    bool hasApertureSet() const;

    /**
     * @brief Box enclosing a band around a sector arc in the chart of toBendArcCoords().
     *
     * The band is |x| <= halfX, |y| <= halfY in arc coordinates for arc length s in
     * [sBegin, sEnd]: straight along the entrance tangent for s < 0, curved with
     * @p curvature on [0, bodyLength], and straight along the exit tangent beyond.
     */
    static BoundingBox getArcSupportBox(
            double curvature, double bodyLength, double sBegin, double sEnd, double halfX,
            double halfY);

    /// The band of getArcSupportBox() as oriented boxes: the straight pieces before and after
    /// the body, and the body in segments of at most 0.1 rad, each enclosing its annular piece.
    static std::vector<SupportBox> getArcSupportBoxes(
            double curvature, double bodyLength, double sBegin, double sEnd, double halfX,
            double halfY);

    /// The axis-aligned box @p box as a SupportBox.
    static SupportBox toSupportBox(const BoundingBox& box);

    // Sharable flag.
    // If this flag is true, the element is always shared.
    mutable bool shareFlag;

    // --- Coordinate system & placement ---
    CoordinateSystemTrafo csTrafoGlobal2Local_m;
    CoordinateSystemTrafo misalignment_m;
    double rotationZAxis_m;

    // --- Aperture ---
    std::pair<ApertureType, std::vector<double>> aperture_m;
    // Default aperture - Needs to be changed to Kokkos::View
    static const std::vector<double> defaultAperture_m;

    // --- Field / physics ---
    // The reference bunch (not owned). Only Monitor still reads it, to save the bunch when
    // the reference particle crosses a TEMPORAL monitor.
    PartBunch_t* RefPartBunch_m;
    bool online_m;

    /// Rest mass (eV) and charge (elementary charges) of the first beam, stored by
    /// takeDesignParticle(); the charge stays 1 for an element without a bunch.
    double designMass_m   = 0.0;
    double designCharge_m = 1.0;
    /// Store the first beam's mass and charge (a charge of zero keeps 1).
    void takeDesignParticle(const PartBunch_t* bunch);

private:
    // Not implemented.
    void operator=(const ElementBase&);

    // --- Identity ---
    // The element's name
    std::string elementID;
    static const std::map<ElementType, std::string> elementTypeToString_s;

    // --- User-defined attributes ---
    AttributeSet userAttribs;

    // --- Placement ---
    bool positionIsFixed;
    double elementPosition_m;  // S position of the element entrance
    bool elemedgeSet_m;

    // --- Miscellaneous ---
    std::string outputfn_m; /**< The name of the outputfile*/
    bool deleteOnTransverseExit_m = true;
    std::optional<double> reach_m;

    /// Position in the tracked line, counted from 1 by OpalBeamline; 0 if not in a line.
    std::size_t lineIndex_m = 0;

public:
    /// Position in the tracked line, counted from 1; 0 for an element outside a tracked line.
    std::size_t getLineIndex() const { return lineIndex_m; }
    /// Set by OpalBeamline when it takes its copy of the element.
    void setLineIndex(std::size_t index) { lineIndex_m = index; }
};

/// Element without a bunch field: apply(pc) stays the no-op, so the tracker skips the call
/// and the transform.
class FramelessElement : public ElementBase {
public:
    using ElementBase::apply;
    using ElementBase::ElementBase;
    ElementFrame fieldFrame() const final { return ElementFrame::None; }
    void apply(const std::shared_ptr<ParticleContainer_t>&, const StepTime&) final {}
};

/**
 * @brief Orders elements by their position in the tracked line.
 *
 * Loops over the elements at one point sum their fields; with this order the sum is the same in
 * every run. Ordering by address made the last digits of the results change from run to run.
 * Elements outside a tracked line (index 0) fall back to their addresses.
 */
struct ByLinePosition {
    bool operator()(
            const std::shared_ptr<ElementBase>& a, const std::shared_ptr<ElementBase>& b) const {
        // A null pointer sorts first, by address like before.
        if (a && b && a->getLineIndex() != b->getLineIndex())
            return a->getLineIndex() < b->getLineIndex();
        return std::less<std::shared_ptr<ElementBase>>()(a, b);
    }
};

/// A set of elements in line order; see ByLinePosition.
using ElementSet = std::set<std::shared_ptr<ElementBase>, ByLinePosition>;

// Inline functions.
// ------------------------------------------------------------------------

/* ============================ Lifecycle & identity ========================= */

inline bool ElementBase::isSharable() const { return shareFlag; }

inline std::string ElementBase::getTypeString() const { return getTypeString(getType()); }

/* ================================ Geometry ================================= */

inline bool ElementBase::isInside(const Vector_t<double, 3>& r) const {
    return getRegions().inFieldRegion(r);
}

/* ===================== Coordinate system & placement ====================== */

inline void ElementBase::setCSTrafoGlobal2Local(const CoordinateSystemTrafo& trafo) {
    if (positionIsFixed) return;

    csTrafoGlobal2Local_m = trafo;
}

inline CoordinateSystemTrafo ElementBase::getCSTrafoGlobal2Local() const {
    return csTrafoGlobal2Local_m;
}

inline void ElementBase::setMisalignment(const CoordinateSystemTrafo& cst) { misalignment_m = cst; }

inline CoordinateSystemTrafo ElementBase::getMisalignment() const { return misalignment_m; }

inline void ElementBase::releasePosition() { positionIsFixed = false; }

inline void ElementBase::fixPosition() { positionIsFixed = true; }

inline bool ElementBase::isPositioned() const { return positionIsFixed; }

inline void ElementBase::setRotationAboutZ(double rotation) { rotationZAxis_m = rotation; }

inline double ElementBase::getRotationAboutZ() const { return rotationZAxis_m; }

inline void ElementBase::setElementPosition(double elemedge) {
    elementPosition_m = elemedge;
    elemedgeSet_m     = true;
}

inline double ElementBase::getElementPosition() const {
    if (elemedgeSet_m) return elementPosition_m;

    throw GeneralOpalException(
            "ElementBase::getElementPosition()",
            std::string("ELEMEDGE for \"") + getName() + "\" not set");
}

inline bool ElementBase::isElementPositionSet() const { return elemedgeSet_m; }

/* ================================ Aperture ================================ */

inline void ElementBase::setAperture(const ApertureType& type, const std::vector<double>& args) {
    aperture_m.first  = type;
    aperture_m.second = args;
}

inline std::pair<ApertureType, std::vector<double>> ElementBase::getAperture() const {
    return aperture_m;
}

/* ========================= Field application & physics ==================== */

inline void ElementBase::setDesignEnergy(const double& /*energy*/, bool /*changeable*/) { return; }

inline double ElementBase::getDesignEnergy() const { return -1.0; }

inline int ElementBase::getRequiredNumberOfTimeSteps() const { return 10; }

inline void ElementBase::setFlagDeleteOnTransverseExit(bool flag) {
    deleteOnTransverseExit_m = flag;
}

inline bool ElementBase::getFlagDeleteOnTransverseExit() const { return deleteOnTransverseExit_m; }

inline void ElementBase::setReach(double reach) { reach_m = reach; }

inline std::optional<double> ElementBase::getReach() const { return reach_m; }

#endif  // OPALX_ElementBase_HH
