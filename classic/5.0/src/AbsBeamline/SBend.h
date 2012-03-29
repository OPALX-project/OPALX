#ifndef CLASSIC_SBend_HH
#define CLASSIC_SBend_HH

// ------------------------------------------------------------------------
// $RCSfile: SBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
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
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Algorithms/PartPusher.h"
#include "Physics/Physics.h"
#include <vector>


class Fieldmap;

// Class SBend
// ------------------------------------------------------------------------
/// Interface for sector bend.
//  This class defines the abstract interface for sector bend magnets.
//  A sector bend magnet has a curved geometry about which its multipole
//  components are specified.

class SBend: public Component {

public:

    /// Constructor with given name.
    explicit SBend(const std::string &name);

    SBend();
    SBend(const SBend &);
    virtual ~SBend();

    /// Apply visitor to SBend.
    virtual void accept(BeamlineVisitor &) const;

    /// Get dipole field of SBend.
    virtual double getB() const = 0;

    /// Get SBend geometry.
    //  Version for non-constant object.
    virtual PlanarArcGeometry &getGeometry() = 0;

    /// Get SBend geometry
    //  Version for constant object.
    virtual const PlanarArcGeometry &getGeometry() const = 0;

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

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    void setBendAngle(const double &);
    void setAmplitudem(double vPeak);

    void setFullGap(const double &);
    void setLength(const double &);

    void setLongitudinalRotation(const double &);
    void setLongitudinalRotation(const double &, const double &);

    void setFaceAngleEntry(double angle);

    void setFaceAngleExit(double angle);
    /// Set the name of the "field map";
    /// For a bend this file contains the coefficients for the Enge function
    void setFieldMapFN(std::string fmapfn);

    std::string getFieldMapFN() const;

    void setAlpha(const double &alpha);

    void setBeta(const double &beta);

    void setDesignEnergy(const double &energy);

    void setK1(const double &k1);

    void setExitAngle(const double &exitAngle);

    virtual const std::string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    double getEffectiveLength() const;
    double getEffectiveCenter() const;
    double getBendAngle() const;

    double getStartElement() const;

    double getR() const;

private:

    // Not implemented.
    void operator=(const SBend &);

    // Name of field map that defines magnet.
    std::string filename_m;
    // Magnet field map.
    Fieldmap *fieldmap_m;
    // Flag to turn on
    bool fast_m;

    double ElementEdge_m;
    double startElement_m;
    // Start of magnet field map.
    double startField_m;
    // End of magnet field map.
    double endField_m;
    // Magnet length.
    double length_m;
    // Magnet gap.
    double gap_m;

    // Flag to reinitialize the field the first time the
    // magnet is applied.
    bool reinitialize_m;

    double alpha_m; // Angle between incoming beam and the entrance face of the magnet.
    double exitAngle_m; // Angle between the outgoing, reference trajectory and exit face of the magnet.

    double sin_face_alpha_m; // alpha is the angle between the projection of the normal of the face onto the
    double cos_face_alpha_m; // s-u plane
    double tan_face_alpha_m;

    double sin_face_beta_m; // beta is defined as arctan(|n_parallel|/|n_perpendicular|) where n_parallel is the component
    double cos_face_beta_m; // of the normal of the face which is parallel to the s-u plane and n_perpendicular is
    double tan_face_beta_m; // perpendicular to it

    double gradient_m;
    // Bend angle of central trajectory.

    // Amplitude of magnet field (Tesla).
    double amplitude_m;
    double angle_m;
    double R_m;
    double design_energy_m;

    Vektor<double, 2> field_orientation_m;

    double *map_m;
    int map_size_m;
    double map_step_size_m;
    // Effective length of hard edge dipole approximation.
    // This is in s-coordinates, so is the effective length
    // along the design energy central trajectory.
    double effectiveLength_m;
    // Effective center of hard edge dipole approximation.
    // This is defined as the s point at which the bend angle is
    // one half the full bend angle.
    double effectiveCenter_m;

    BorisPusher pusher_m;

    void calculateEffectiveLength();
    void calculateEffectiveCenter();

    // Set the bend strength based on the desired bend angle
    // and the reference energy.
    void setBendStrength();

    // Calculate bend angle.
    double calculateBendAngle(double bendLength, bool modifyField);

    // Calculate the reference particle trajectory map.
    double calculateRefTrajectory(const double zBegin);
};

inline void SBend::setAlpha(const double &alpha) {
    alpha_m = alpha * Physics::pi / 180.0;
    Orientation_m(0) = -alpha * Physics::pi / 180.0;
    //    Orientation_m(0) = 0.0;
    sin_face_alpha_m = sin(Orientation_m(0));
    cos_face_alpha_m = cos(Orientation_m(0));
    tan_face_alpha_m = tan(Orientation_m(0));
}

inline void SBend::setBeta(const double &beta) {
    Orientation_m(1) = beta * Physics::pi / 180.0;
    sin_face_beta_m = sin(Orientation_m(1));
    cos_face_beta_m = cos(Orientation_m(1));
    tan_face_beta_m = tan(Orientation_m(1));
}

inline void SBend::setK1(const double &k1) {
    gradient_m = k1;
}

inline void SBend::setExitAngle(const double &exitAngle) {
    exitAngle_m = exitAngle * Physics::pi / 180.0;
}

inline void SBend::setDesignEnergy(const double &energy)
{ design_energy_m = energy; }

inline void SBend::setLongitudinalRotation(const double &rotation) {
    Orientation_m(2) = rotation * Physics::pi / 180.0;
}

inline void SBend::setLongitudinalRotation(const double &k0, const double &k0s) {

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

#endif // CLASSIC_SBend_HH
