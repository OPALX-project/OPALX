/**
 * @file Component.h
 * @brief Defines the abstract interface for a single beamline component in the
 * accelerator model.
 *
 * @details
 * The Component class provides an abstract interface for an arbitrary single
 * component in a beam line.
 * A component is the basic element in the accelerator model, such as a dipole,
 * quadrupole, etc.
 */
#ifndef CLASSIC_Component_HH
#define CLASSIC_Component_HH

#include "AbsBeamline/ElementBase.h"
#include "Fields/EMField.h"
#include "OPALTypes.h"

using ParticleContainer_t = ParticleContainer<double, 3>;

class PartData;

template <class T, int N>
class FVps;

struct Point {
    double x;
    double y;
};

class Component : public ElementBase {
public:
    /* ============================== Constructors ============================== */
    explicit Component(const std::string& name);
    Component();
    Component(const Component& right);
    virtual ~Component();
    /* ========================================================================== */
    /* ============================== Field Functions =========================== */
    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for non-constant object).
    virtual EMField& getField() = 0;

    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for constant object).
    virtual const EMField& getField() const = 0;

    /// Return the field in a point.
    //  Return the value of the time-independent part of the electric
    //  field at point [b]P[/b].
    EVector Efield(const Point3D& P) const;

    /// Return the field in a point.
    //  Return the value of the time-independent part of the magnetic
    //  field at point [b]P[/b].
    BVector Bfield(const Point3D& P) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of the electric
    //  field at point [b]P[/b] for time [b]t[/b].
    EVector Efield(const Point3D& P, double t) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of the magnetic
    //  field at point [b]P[/b] for time [b]t[/b].
    BVector Bfield(const Point3D& P, double t) const;

    /// Return the field in a point.
    //  Return the value of the time-independent part of both electric
    //  and magnetic fields at point [b]P[/b].
    EBVectors EBfield(const Point3D& P) const;

    /// Return the field in a point.
    //  Return the value of the time-dependent part of both electric
    //  and magnetic fields at point [b]P[/b] for time [b]t[/b].
    EBVectors EBfield(const Point3D& P, double t) const;
    /* ========================================================================== */
    /* ============================== Apply Functions =========================== */
    /**
     * Apply functions apply the components electromagnetic field to the
     * particles.
     * They are called inside ParalleTracker::computeExternalFields()
     */

    /**
     * @brief Apply to all particles. Kernel launch moved inside the function.
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply();

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
        const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B
    );

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
        Vector_t<double, 3>& E, Vector_t<double, 3>& B
    );

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
        Vector_t<double, 3>& E, Vector_t<double, 3>& B
    );

    /* ========================================================================== */
    /* ============================== Functions ================================= */
    /**
     * @brief Calculate the four-potential at some position relative to the
     * component
     *
     * @param R position in the local coordinate system of the component
     * @param t time
     * @param A filled with the calculated magnetic vector potential
     * @param phi filled with the calculated electric potential
     * Note that any existing values in A and phi may be overwritten by this
     * method.
     *
     * @returns true if particle is outside the field map, else false
     * Default for component is to return false and make no change to A and phi
     */
    virtual bool getPotential(
        const Vector_t<double, 3>& /*R*/, const double& /*t*/, Vector_t<double, 3>& /*A*/,
        double& /*phi*/
    ) {
        return false;
    }

    // Design energy for components such as RF-cavities
    virtual double getDesignEnergy() const;
    virtual void setDesignEnergy(const double& energy, bool changeable = true);

    // Setup
    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) = 0;

    // Clean-up
    virtual void finalise() = 0;

    // Does the component bend?
    virtual bool bends() const = 0;

    // Read & free fieldmaps
    virtual void goOnline(const double& kineticEnergy);
    virtual void goOffline();

    // Is the component online (been initialised)?
    virtual bool Online();

    // Gets the size (along Z) of the component
    virtual void getDimensions(double& zBegin, double& zEnd) const = 0;

    // Gets the element type defined in ElementBase
    virtual ElementType getType() const;

    /// Return design element.
    //  If this method returns a pointer to this component.
    //  The default version returns ``this''.
    virtual const ElementBase& getDesign() const;

    /// Track particle bunch.
    //  This catch-all method implements a hook for tracking a particle
    //  bunch through a non-standard component.
    //  The default version throws a LogicalError.
    virtual void trackBunch(PartBunch_t* bunch, const PartData&, bool revBeam, bool revTrack) const;

    /// Track a map.
    //  This catch-all method implements a hook for tracking a transfer
    //  map through a non-standard component.
    //  The default version throws a LogicalError.
    virtual void trackMap(FVps<double, 6>& map, const PartData&, bool revBeam, bool revTrack) const;

    void setExitFaceSlope(const double&);
    /* ========================================================================== */
protected:
    // Aperture - Needs to be changed to Kokkos::View
    static const std::vector<double> defaultAperture_m;
    double exit_face_slope_m;

    // The reference bunch
    PartBunch_t* RefPartBunch_m;
    bool online_m;
};

// Inline access functions to fields.
// ------------------------------------------------------------------------

inline EVector Component::Efield(const Point3D& P) const { return getField().Efield(P); }

inline BVector Component::Bfield(const Point3D& P) const { return getField().Bfield(P); }

inline EVector Component::Efield(const Point3D& P, double t) const {
    return getField().Efield(P, t);
}

inline BVector Component::Bfield(const Point3D& P, double t) const {
    return getField().Bfield(P, t);
}

inline EBVectors Component::EBfield(const Point3D& P) const { return getField().EBfield(P); }

inline EBVectors Component::EBfield(const Point3D& P, double t) const {
    return getField().EBfield(P, t);
}

inline void Component::setExitFaceSlope(const double& m) { exit_face_slope_m = m; }

inline void Component::setDesignEnergy(const double& /*energy*/, bool /*changeable*/) { return; }

inline double Component::getDesignEnergy() const { return -1.0; }

#endif  // CLASSIC_Component_HH
