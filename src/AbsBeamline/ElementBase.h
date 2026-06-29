//
// Class ElementBase
//   The very base class for beam line representation objects. A beam line
//   is modelled as a composite structure having a single root object
//   (the top level beam line), which contains both ``single'' leaf-type
//   elements (Components), as well as sub-lines (composites).
//
//   Interface for basic beam line object.
//   This class defines the abstract interface for all objects which can be
//   contained in a beam line. ElementBase forms the base class for two distinct
//   but related hierarchies of objects:
//   [OL]
//   [LI]
//   A set of concrete accelerator element classes, which compose the standard
//   accelerator component library (SACL).
//   [LI]
//   [/OL]
//   Instances of the concrete classes for single elements are by default
//   sharable. Instances of beam lines and integrators are by
//   default non-sharable, but they may be made sharable by a call to
//   [b]makeSharable()[/b].
//   [p]
//   An ElementBase object can return two lengths, which may be different:
//   [OL]
//   [LI]
//   The arc length along the geometry.
//   [LI]
//   The design length, often measured along a straight line.
//   [/OL]
//   Class ElementBase contains a map of name versus value for user-defined
//   attributes (see file AbsBeamline/AttributeSet.hh). The map is primarily
//   intended for processes that require algorithm-specific data in the
//   accelerator model.
//   [P]
//   The class ElementBase is a base class.
//   Virtual derivation is used to make multiple inheritance possible for
//   derived concrete classes. ElementBase implements three copy modes:
//   [OL]
//   [LI]
//   Copy by reference: Use std::shared_ptr to share ownership.
//   [LI]
//   Copy structure: use ElementBase::copyStructure().
//   During copying of the structure, all sharable items are re-used, while
//   all non-sharable ones are cloned.
//   [LI]
//   Copy by cloning: use ElementBase::clone().
//   This returns a full deep copy.
//   [/OL]
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
#include <queue>
#include <string>
#include <utility>
#include <vector>

class BeamlineVisitor;
class BoundaryGeometry;
class Channel;
class ConstChannel;
class ParticleMatterInteractionHandler;
class WakeFunction;
class PartData;

template <class T, int N>
class FVps;

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

class ElementBase : public std::enable_shared_from_this<ElementBase> {
public:
    /// Constructor with given name.
    explicit ElementBase(const std::string& name);

    ElementBase();
    ElementBase(const ElementBase&);
    virtual ~ElementBase();

    /// Get element name.
    virtual const std::string& getName() const;

    /// Set element name.
    virtual void setName(const std::string& name);

    /// Get element type std::string.
    //  Default returns ElementType::ANY; concrete elements override.
    virtual ElementType getType() const;

    std::string getTypeString() const;
    static std::string getTypeString(ElementType type);

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual BGeometryBase& getGeometry() = 0;

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const BGeometryBase& getGeometry() const = 0;

    /// Get arc length.
    //  Return the entire arc length measured along the design orbit
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the design length defined by the geometry.
    //  This may be the arc length or the straight length.
    virtual double getElementLength() const;

    /// Set design length.
    //  Set the design length defined by the geometry.
    //  This may be the arc length or the straight length.
    virtual void setElementLength(double length);

    /**
     * @brief Return the nominal body extent of the element.
     *
     * The first placement redesign stage distinguishes between the nominal body
     * extent and the field-support extent. The body extent is the canonical
     * longitudinal interval of the placed hardware,
     * \f$[z_\mathrm{body}^{\mathrm{begin}}, z_\mathrm{body}^{\mathrm{end}}]\f$,
     * and therefore drives ports, placement, and visualization. By default it
     * coincides with the geometry length \f$[0, L]\f$ in the local chart.
     */
    virtual void getElementDimensions(double& begin, double& end) const {
        begin = 0.0;
        end   = getElementLength();
    }

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the element
    //  (origin >= 0)
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the element
    //  (entrance <= 0)
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the element
    //  (exit >= 0)
    virtual double getExit() const;

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

    /// Apply visitor.
    //  This method must be overridden by derived classes. It should call the
    //  method of the visitor corresponding to the element class.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor& visitor) const = 0;

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const = 0;

    /// Make a structural copy.
    //  Return a fresh copy of any beam line structure is made,
    //  but sharable elements remain shared.
    virtual ElementBase* copyStructure();

    /* ============================== Apply Functions =========================== */
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
     * @brief Apply to particle i
     *
     * @param i Particle index
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply(
            const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B);

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

    /**
     * @brief Calculate the four-potential at some position relative to the
     * element
     *
     * @param R position in the local coordinate system of the element
     * @param t time
     * @param A filled with the calculated magnetic vector potential
     * @param phi filled with the calculated electric potential
     * Note that any existing values in A and phi may be overwritten by this
     * method.
     *
     * @returns true if particle is outside the field map, else false
     * Default is to return false and make no change to A and phi
     */
    virtual bool getPotential(
            const Vector_t<double, 3>& /*R*/, const double& /*t*/, Vector_t<double, 3>& /*A*/,
            double& /*phi*/) {
        return false;
    }

    // Design energy for elements such as RF-cavities
    virtual double getDesignEnergy() const;
    virtual void setDesignEnergy(const double& energy, bool changeable = true);

    // Setup
    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) = 0;

    // Clean-up
    virtual void finalise() = 0;

    // Does the element bend?
    virtual bool bends() const = 0;

    /// @name Bend queries
    /// Overridden by SBend/RBend; straight-element defaults otherwise. Let callers
    /// query bend geometry through an ElementBase pointer without downcasting.
    ///@{
    virtual double getBendAngle() const;
    virtual double getEntranceAngle() const;
    virtual double getChordLength() const;
    virtual std::vector<Vector_t<double, 3>> getDesignPath(std::size_t minSamples = 32) const;
    ///@}

    // Read & free fieldmaps
    virtual void goOnline(const double& kineticEnergy);
    virtual void goOffline();

    // Is the element online (been initialised)?
    virtual bool Online();

    /**
     * @brief Return the field-support extent of the element.
     *
     * This is the longitudinal interval
     * \f$[z_\mathrm{field}^{\mathrm{begin}}, z_\mathrm{field}^{\mathrm{end}}]\f$
     * on which the external field model is evaluated in the element-local
     * chart. In the first extent model this may differ from the nominal body
     * extent returned by `getElementDimensions()`, for example when fringe
     * fields extend beyond the hardware body or when a field map occupies only
     * part of the body.
     */
    virtual void getFieldExtend(double& zBegin, double& zEnd) const = 0;

    /**
     * @brief Track a borrowed particle bunch through a non-standard element.
     *
     * The default implementation throws a LogicalError.
     *
     * @param bunch Particle bunch to track. The element does not take ownership.
     */
    virtual void trackBunch(PartBunch_t& bunch, const PartData&, bool revBeam, bool revTrack) const;

    /// Track a map.
    //  This catch-all method implements a hook for tracking a transfer
    //  map through a non-standard element.
    //  The default version throws a LogicalError.
    virtual void trackMap(FVps<double, 6>& map, const PartData&, bool revBeam, bool revTrack) const;

    void setExitFaceSlope(const double&);
    /* ========================================================================== */

    /// Test if the element can be shared.
    bool isSharable() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Update element.
    //  This method stores all attributes contained in the AttributeSet to
    //  "*this".  The return value [b]true[/b] indicates success.
    bool update(const AttributeSet&);

    ///@{ Access to ELEMEDGE attribute
    void setElementPosition(double elemedge);
    double getElementPosition() const;
    bool isElementPositionSet() const;
    ///@}
    /// attach a boundary geometry field to the element
    virtual void setBoundaryGeometry(BoundaryGeometry* geo);

    /// return the attached boundary geometrt object if there is any
    virtual BoundaryGeometry* getBoundaryGeometry() const;

    virtual bool hasBoundaryGeometry() const;

    /// attach a wake field to the element
    virtual void setWake(WakeFunction* wf);

    /// return the attached wake object if there is any
    virtual WakeFunction* getWake() const;

    virtual bool hasWake() const;

    virtual void setParticleMatterInteraction(ParticleMatterInteractionHandler* spys);

    virtual ParticleMatterInteractionHandler* getParticleMatterInteraction() const;

    virtual bool hasParticleMatterInteraction() const;

    void setCSTrafoGlobal2Local(const CoordinateSystemTrafo& ori);
    CoordinateSystemTrafo getCSTrafoGlobal2Local() const;
    void releasePosition();
    void fixPosition();
    bool isPositioned() const;

    /// Body→entrance-edge transform of the element's local chart (identity for
    /// straight elements; overridden by bends).
    virtual CoordinateSystemTrafo getEdgeToBegin() const;
    /// Body→exit-edge transform of the element's local chart (a +length shift in
    /// local z for straight elements; overridden by bends).
    virtual CoordinateSystemTrafo getEdgeToEnd() const;

    void setAperture(const ApertureType& type, const std::vector<double>& args);
    std::pair<ApertureType, std::vector<double> > getAperture() const;

    virtual bool isInside(const Vector_t<double, 3>& r) const;

    void setMisalignment(const CoordinateSystemTrafo& cst);

    void getMisalignment(double& x, double& y, double& s) const;
    CoordinateSystemTrafo getMisalignment() const;

    void setActionRange(const std::queue<std::pair<double, double> >& range);
    void setCurrentSCoordinate(double s);

    /// Set rotation about z axis in bend frame.
    void setRotationAboutZ(double rotation);
    double getRotationAboutZ() const;

    virtual BoundingBox getBoundingBoxInLabCoords() const;

    virtual int getRequiredNumberOfTimeSteps() const;

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

    CoordinateSystemTrafo csTrafoGlobal2Local_m;
    CoordinateSystemTrafo misalignment_m;

    std::pair<ApertureType, std::vector<double> > aperture_m;

    double elementEdge_m;

    double rotationZAxis_m;

    // Default aperture - Needs to be changed to Kokkos::View
    static const std::vector<double> defaultAperture_m;
    double exit_face_slope_m;

    // The reference bunch (not owned)
    PartBunch_t* RefPartBunch_m;
    bool online_m;

private:
    // Not implemented.
    void operator=(const ElementBase&);

    // The element's name
    std::string elementID;

    static const std::map<ElementType, std::string> elementTypeToString_s;

    // The user-defined set of attributes.
    AttributeSet userAttribs;

    WakeFunction* wake_m;

    BoundaryGeometry* bgeometry_m;

    ParticleMatterInteractionHandler* parmatint_m;

    bool positionIsFixed;
    ///@{ ELEMEDGE attribute
    double elementPosition_m;
    bool elemedgeSet_m;
    ///@}
    std::queue<std::pair<double, double> > actionRange_m;

    std::string outputfn_m; /**< The name of the outputfile*/

    bool deleteOnTransverseExit_m = true;
};

// Inline functions.
// ------------------------------------------------------------------------

inline double ElementBase::getArcLength() const { return getGeometry().getArcLength(); }

inline double ElementBase::getElementLength() const { return getGeometry().getElementLength(); }

inline void ElementBase::setElementLength(double length) { getGeometry().setElementLength(length); }

inline double ElementBase::getOrigin() const { return getGeometry().getOrigin(); }

inline double ElementBase::getEntrance() const { return getGeometry().getEntrance(); }

inline double ElementBase::getExit() const { return getGeometry().getExit(); }

inline bool ElementBase::isSharable() const { return shareFlag; }

inline WakeFunction* ElementBase::getWake() const { return wake_m; }

inline bool ElementBase::hasWake() const { return wake_m != nullptr; }

inline BoundaryGeometry* ElementBase::getBoundaryGeometry() const { return bgeometry_m; }

inline bool ElementBase::hasBoundaryGeometry() const { return bgeometry_m != nullptr; }

inline ParticleMatterInteractionHandler* ElementBase::getParticleMatterInteraction() const {
    return parmatint_m;
}

inline bool ElementBase::hasParticleMatterInteraction() const { return parmatint_m != nullptr; }

inline void ElementBase::setCSTrafoGlobal2Local(const CoordinateSystemTrafo& trafo) {
    if (positionIsFixed) return;

    csTrafoGlobal2Local_m = trafo;
}

inline CoordinateSystemTrafo ElementBase::getCSTrafoGlobal2Local() const {
    return csTrafoGlobal2Local_m;
}

inline CoordinateSystemTrafo ElementBase::getEdgeToBegin() const {
    return getGeometry().getEdgeToBegin();
}

inline CoordinateSystemTrafo ElementBase::getEdgeToEnd() const {
    return getGeometry().getEdgeToEnd();
}

inline void ElementBase::setAperture(const ApertureType& type, const std::vector<double>& args) {
    aperture_m.first  = type;
    aperture_m.second = args;
}

inline std::pair<ApertureType, std::vector<double> > ElementBase::getAperture() const {
    return aperture_m;
}

inline bool ElementBase::isInside(const Vector_t<double, 3>& r) const {
    const double length = getElementLength();
    return r(2) >= 0.0 && r(2) < length && isInsideTransverse(r);
}

inline void ElementBase::setMisalignment(const CoordinateSystemTrafo& cst) { misalignment_m = cst; }

inline CoordinateSystemTrafo ElementBase::getMisalignment() const { return misalignment_m; }

inline void ElementBase::releasePosition() { positionIsFixed = false; }

inline void ElementBase::fixPosition() { positionIsFixed = true; }

inline bool ElementBase::isPositioned() const { return positionIsFixed; }

inline void ElementBase::setActionRange(const std::queue<std::pair<double, double> >& range) {
    actionRange_m = range;

    if (!actionRange_m.empty()) elementEdge_m = actionRange_m.front().first;
}

inline void ElementBase::setRotationAboutZ(double rotation) { rotationZAxis_m = rotation; }

inline double ElementBase::getRotationAboutZ() const { return rotationZAxis_m; }

inline std::string ElementBase::getTypeString() const { return getTypeString(getType()); }

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

inline int ElementBase::getRequiredNumberOfTimeSteps() const { return 10; }

inline void ElementBase::setFlagDeleteOnTransverseExit(bool flag) {
    deleteOnTransverseExit_m = flag;
}

inline bool ElementBase::getFlagDeleteOnTransverseExit() const { return deleteOnTransverseExit_m; }

inline double ElementBase::getBendAngle() const { return 0.0; }

inline double ElementBase::getEntranceAngle() const { return 0.0; }

inline double ElementBase::getChordLength() const { return getElementLength(); }

inline std::vector<Vector_t<double, 3>> ElementBase::getDesignPath(std::size_t) const {
    return {Vector_t<double, 3>({0.0, 0.0, 0.0}),
            Vector_t<double, 3>({0.0, 0.0, getElementLength()})};
}

inline void ElementBase::setExitFaceSlope(const double& m) { exit_face_slope_m = m; }

inline void ElementBase::setDesignEnergy(const double& /*energy*/, bool /*changeable*/) { return; }

inline double ElementBase::getDesignEnergy() const { return -1.0; }

#endif  // OPALX_ElementBase_HH
