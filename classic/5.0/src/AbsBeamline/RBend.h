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
#include "Fields/Fieldmap.hh"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Algorithms/PartPusher.h"
#include <vector>


// Class RBend
// ------------------------------------------------------------------------
/// Interface for rectangular bend.
//  Class RBend defines the abstract interface for rectangular bend magnets.
//  A rectangular bend magnet has a rectilinear geometry about which its
//  multipole components are specified.

class RBend: public Component {

public:

    /// Constructor with given name.
    explicit RBend(const string &name);

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

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);
  
    virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual void rescaleFieldMap(const double &scaleFactor);

    virtual bool bends() const;

    void setAmplitudem(double vPeak);

    void setAngle(const double&, const double&);

    /// Set the name of the "field map"; 
    /// For a bend this file contains the coefficients for the Enge function 
    void setFieldMapFN(string fmapfn);

    string getFieldMapFN() const;

    void setEngeCoefs(const vector<double> EngeCoefs);

    void setAlpha(const double &alpha);

    void setBeta(const double &beta);

    void setDesignEnergy(const double &energy);

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    double getR() const;

private:

    // Not implemented.
    void operator=(const RBend &);
    string filename_m;             /**< The name of the inputfile*/
    Fieldmap *fieldmap_m;
    double amplitude_m;
    Vektor<double, 2> field_orientation_m;   // (cos(alpha), sin(alpha))
    double startField_m;
    double endField_m;
    double length_m;  

    bool fast_m;

    double sin_face_alpha_m; // alpha is the angle between the projection of the normal of the face onto the
    double cos_face_alpha_m; // s-u plane
    double tan_face_alpha_m;

    double sin_face_beta_m; // beta is defined as arctan(|n_parallel|/|n_perpendicular|) where n_parallel is the component
    double cos_face_beta_m; // of the normal of the face which is parallel to the s-u plane and n_perpendicular is
    double tan_face_beta_m; // perpendicular to it

    double design_energy_m;
    double *map_m;
    int map_size_m;
    double map_step_size_m;
    BorisPusher pusher_m;

    double R_m;
};

inline void RBend::setAlpha(const double &alpha)
{
    Orientation_m(0) = alpha * Physics::pi / 180.0;
    sin_face_alpha_m = sin(Orientation_m(0));
    cos_face_alpha_m = cos(Orientation_m(0));
    tan_face_alpha_m = tan(Orientation_m(0));
}

inline void RBend::setBeta(const double &beta)
{
    Orientation_m(1) = beta * Physics::pi / 180.0;
    sin_face_beta_m = sin(Orientation_m(1));
    cos_face_beta_m = cos(Orientation_m(1));
    tan_face_beta_m = tan(Orientation_m(1));
}

inline void RBend::setDesignEnergy(const double &energy)
{ design_energy_m = energy; }

inline void RBend::setAngle(const double& k0, const double& k0s) 
{
    if (fabs(k0s) > 1.e-8) {
        amplitude_m = sqrt(k0*k0 + k0s*k0s);
        if (amplitude_m > 1.e-8) {
            field_orientation_m(0) = k0 / amplitude_m;
            field_orientation_m(1) = k0s / amplitude_m;
        }
    } else {
        amplitude_m = k0;
        field_orientation_m(0) = 1.0;
        field_orientation_m(1) = 0.0;
    }
}
#endif // CLASSIC_RBend_HH
