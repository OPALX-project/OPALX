#ifndef CLASSIC_RBend_HH
#define CLASSIC_RBend_HH

// ------------------------------------------------------------------------
// $RCSfile: RBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Algorithms/PartPusher.h"
#include "Fields/Fieldmap.hh"
#include <vector>

class Fieldmap;

// Class RBend
// ------------------------------------------------------------------------
/// Interface for rectangular bend.
//  Class RBend defines the abstract interface for rectangular bend magnets.
//  A rectangular bend magnet has a rectilinear geometry about which its
//  multipole components are specified.

class RBend: public Component {

public:

    /// Constructor with given name.
    explicit RBend(const std::string &name);

    RBend();
    RBend(const RBend &);
    virtual ~RBend();

    /// Apply visitor to RBend.
    virtual void accept(BeamlineVisitor &) const;

    /// Get dipole field of RBend.
    virtual double getB() const = 0;

    /// Get RBend geometry.
    //  Version for non-constant object.
    virtual RBendGeometry &getGeometry() = 0;

    /// Get RBend geometry
    //  Version for constant object.
    virtual const RBendGeometry &getGeometry() const = 0;

    /// Get multipole expansion of field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField() = 0;

    /// Get multipole expansion of field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const = 0;

    /// Get normal component.
    //  Return the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getNormalComponent(int) const;

    /// Get skew component.
    //  Return the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getSkewComponent(int) const;

    /// Set normal component.
    //  Set the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setNormalComponent(int, double);

    /// Set skew component.
    //  Set the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setSkewComponent(int, double);

    /// Get pole entry face rotation.
    //  Return the rotation of the entry pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getEntryFaceRotation() const = 0;

    /// Get exit pole face rotation.
    //  Return the rotation of the exit pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getExitFaceRotation() const = 0;

    /// Get entry pole face curvature.
    //  Return the curvature of the entry pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getEntryFaceCurvature() const = 0;

    /// Get exit pole face curvature.
    //  Return the curvature of the exit pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getExitFaceCurvature() const = 0;

    /// Get number of slices.
    //  Slices and stepsize used to determine integration step.
    virtual double getSlices() const = 0;

    /// Get stepsize.
    //  Slices and stepsize used to determine integration step.
    virtual double getStepsize() const = 0;

    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    void setBendAngle(const double &angle);
    void setAmplitudem(double vPeak);

    void setFullGap(const double &gap);
    void setLength(const double &length);

    void setLongitudinalRotation(const double &rotation);
    void setLongitudinalRotation(const double &k0, const double &k0s);

    /// Set the name of the "field map";
    /// For a bend this file contains the coefficients for the Enge function
    void setFieldMapFN(std::string fmapfn);

    std::string getFieldMapFN() const;

    void setEngeCoefs(const std::vector<double> EngeCoefs);

    void setAlpha(const double &alpha);

    void setBeta(const double &beta);

    void setDesignEnergy(const double &energy);

    virtual const std::string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    double getEffectiveLength() const;
    double getEffectiveCenter() const;
    double getBendAngle() const;

    double getStartElement() const;
    double getR() const;

private:

    // Magnet length.
    double length_m;

    // Magnet gap.
    double gap_m;

    // Flag to reinitialize the field the first time the
    // magnet is applied.
    bool reinitialize_m;

    // Not implemented.
    void operator=(const RBend &);
    std::string filename_m;             /**< The name of the inputfile*/
    Fieldmap *myFieldmap;
    double scale_m;                /**< scale multiplier*/
    Fieldmap *fieldmap_m;
    double amplitude_m;
    Vektor<double, 2> field_orientation_m;   // (cos(alpha), sin(alpha))
    double startField_m;
    double endField_m;

    bool fast_m;

    double sin_face_alpha_m; // alpha is the angle between the projection of the normal of the face onto the
    double cos_face_alpha_m; // s-u plane
    double tan_face_alpha_m;

    double sin_face_beta_m; // beta is defined as arctan(|n_parallel|/|n_perpendicular|) where n_parallel is the component
    double cos_face_beta_m; // of the normal of the face which is parallel to the s-u plane and n_perpendicular is
    double tan_face_beta_m; // perpendicular to it

    double design_energy_m;
    double angle_m; // Bend angle of central trajectory particle (radians).
    double *map_m;
    int map_size_m;
    double map_step_size_m;
    BorisPusher pusher_m;
    double startElement_m;
    double R_m;

    // Effective length of hard edge dipole approximation.
    double effectiveLength_m;

    // Effective center of hard edge dipole approximation.
    //
    // This is defined as the z point at which the field integral
    // is half of the full field integral starting at the upstream
    // end of the dipole.
    double effectiveCenter_m;

    static int RBend_counter_m;

    void calculateEffectiveLength();
    void calculateEffectiveCenter();

    // Set the bend strength based on the desired bend angle
    // and the reference energy.
    void setBendStrength();

    // Calculate bend angle given reference energy and field strength.
    double calculateBendAngle(double bendLength);

    // Calculate the reference particle trajectory map.
    double calculateRefTrajectory(const double zBegin);
};

inline void RBend::setAlpha(const double &alpha) {
    Orientation_m(0) = -alpha * Physics::pi / 180.0;
    sin_face_alpha_m = sin(Orientation_m(0));
    cos_face_alpha_m = cos(Orientation_m(0));
    tan_face_alpha_m = tan(Orientation_m(0));
}

inline void RBend::setBeta(const double &beta) {
    Orientation_m(1) = beta * Physics::pi / 180.0;
    sin_face_beta_m = sin(Orientation_m(1));
    cos_face_beta_m = cos(Orientation_m(1));
    tan_face_beta_m = tan(Orientation_m(1));
}

inline void RBend::setDesignEnergy(const double &energy)
{ design_energy_m = energy; }

inline void RBend::setLongitudinalRotation(const double &rotation) {
    Orientation_m(2) = rotation * Physics::pi / 180.0;
}

inline void RBend::setLongitudinalRotation(const double &k0, const double &k0s) {

    Orientation_m(2) = atan2(k0s, k0);
    amplitude_m = sqrt(pow(k0, 2.0) + pow(k0s, 2.0));

    if(fabs(k0s) > 1.e-8) {
        if(amplitude_m > 1.e-8) {
            field_orientation_m(0) = k0 / amplitude_m;
            field_orientation_m(1) = k0s / amplitude_m;
        }
    } else {
        field_orientation_m(0) = 1.0;
        field_orientation_m(1) = 0.0;
    }
}
#endif // CLASSIC_RBend_HH
