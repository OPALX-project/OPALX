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
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"
#include "BeamlineGeometry/Geometry.h"
#include "OPALTypes.h"
#include "Structure/BoundingBox.h"
#include "Utilities/GeneralOpalException.h"
#include "VectorMath.h"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class BeamlineVisitor;
class Channel;
class ConstChannel;

enum class ElementType : unsigned short {
    ANY,
    BEAMLINE,
    BEAMBEAM,
    DRIFT,
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
    RING,
    PROBE,
    VACUUM,
    SOLENOID,
    SOURCE,
    CONSTANTEFIELDCAVITY
};

enum class ApertureType : unsigned short {
    RECTANGULAR,
    ELLIPTICAL,
    CONIC_RECTANGULAR,
    CONIC_ELLIPTICAL
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

    /// @brief Check if the point r is inside the S interval with a field.
    /// @param r The point to test.
    /// @return True if r is inside the field interval.
    virtual bool isInside(const Vector_t<double, 3>& r) const;

    /// @brief Get the bounding box.
    /// @return The bounding box in lab coordinates.
    virtual BoundingBox getBoundingBoxInLabCoords() const;

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

    /* =============================== Aperture ================================ */

    void setAperture(const ApertureType& type, const std::vector<double>& args);
    std::pair<ApertureType, std::vector<double>> getAperture() const;

    /* ===================== Field application & physics ======================= */
    /**
     * Apply functions apply the element's electromagnetic field to the
     * particles. They are called inside ParallelTracker::computeExternalFields().
     */

    /**
     * @brief Apply to all particles. Kernel launch moved inside the function.
     */
    virtual void apply(const std::shared_ptr<ParticleContainer_t>& pc);

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

    // Design energy for elements such as RF-cavities
    virtual double getDesignEnergy() const;
    virtual void setDesignEnergy(const double& energy, bool changeable = true);

    // Setup
    virtual void initialise(PartBunch_t* bunch) = 0;

    // Clean-up
    virtual void finalise() = 0;

    /// Prepare runtime resources for tracking (e.g. load a field map); sets online_m.
    virtual void goOnline(const double& kineticEnergy);
    /// Release runtime resources / flush element output (e.g. loss data); clears online_m.
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

    void setFlagDeleteOnTransverseExit(bool = true);
    bool getFlagDeleteOnTransverseExit() const;

protected:
    bool isInsideTransverse(const Vector_t<double, 3>& r) const;

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
    // The reference bunch (not owned)
    PartBunch_t* RefPartBunch_m;
    bool online_m;

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
};

// Inline functions.
// ------------------------------------------------------------------------

/* ============================ Lifecycle & identity ========================= */

inline bool ElementBase::isSharable() const { return shareFlag; }

inline std::string ElementBase::getTypeString() const { return getTypeString(getType()); }

/* ================================ Geometry ================================= */

inline bool ElementBase::isInside(const Vector_t<double, 3>& r) const {
    // Selection uses the field-support extent (the longitudinal interval where the element's
    // field model is non-zero), not the body length, so a particle in a fringe region beyond
    // the body is still attributed to the element. getFieldExtent() is the single source for
    // that interval, in the same local chart as r.
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);
    return r(2) >= zBegin && r(2) < zEnd && isInsideTransverse(r);
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

#endif  // OPALX_ElementBase_HH
