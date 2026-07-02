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

    /// Constructor with given name.
    explicit ElementBase(const std::string& name);

    /// Default Constructor
    ElementBase();

    /// Copy Constructor
    ElementBase(const ElementBase&);

    /// Destructor
    virtual ~ElementBase();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const = 0;

    /// Make a structural copy.
    //  Return a fresh copy of any beam line structure is made,
    //  but sharable elements remain shared.
    virtual ElementBase* copyStructure();

    /// Test if the element can be shared.
    bool isSharable() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /* =============================== Identity ================================ */

    /// Get element name.
    virtual const std::string& getName() const;

    /// Set element name.
    virtual void setName(const std::string& name);

    /// Get element type std::string.
    //  Default returns ElementType::ANY; concrete elements override.
    virtual ElementType getType() const;

    std::string getTypeString() const;
    static std::string getTypeString(ElementType type);

    /// Apply visitor.
    //  This method must be overridden by derived classes. It should call the
    //  method of the visitor corresponding to the element class.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor& visitor) const = 0;

    /* =============================== Geometry ================================ */

    /// Get geometry.
    //  Return the element geometry, supplied by the representation layer.
    //  Version for non-constant object.
    virtual Geometry& getGeometry() = 0;

    /// Get geometry. Version for constant object.
    virtual const Geometry& getGeometry() const = 0;

    ///
    virtual bool isInside(const Vector_t<double, 3>& r) const;

    virtual BoundingBox getBoundingBoxInLabCoords() const;

    /* ===================== Coordinate system & placement ===================== */

    void setCSTrafoGlobal2Local(const CoordinateSystemTrafo& ori);
    CoordinateSystemTrafo getCSTrafoGlobal2Local() const;

    void setMisalignment(const CoordinateSystemTrafo& cst);
    void getMisalignment(double& x, double& y, double& s) const;
    CoordinateSystemTrafo getMisalignment() const;

    void releasePosition();
    void fixPosition();
    bool isPositioned() const;

    /// Set rotation about z axis in bend frame.
    void setRotationAboutZ(double rotation);
    double getRotationAboutZ() const;

    ///@{ Access to ELEMEDGE attribute
    void setElementPosition(double elemedge);
    double getElementPosition() const;
    bool isElementPositionSet() const;
    ///@}

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
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply(const std::shared_ptr<ParticleContainer_t>& pc);

    /**
     * @brief Apply to particle with position R and momentum P
     *
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply(
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
    ///@{ ELEMEDGE attribute
    double elementPosition_m;
    bool elemedgeSet_m;
    ///@}

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
