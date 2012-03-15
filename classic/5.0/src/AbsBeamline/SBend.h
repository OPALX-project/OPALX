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


// Class SBend
// ------------------------------------------------------------------------
/// Interface for sector bend.
//  This class defines the abstract interface for sector bend magnets.
//  A sector bend magnet has a curved geometry about which its multipole
//  components are specified.

class SBend: public Component {

public:

    /// Constructor with given name.
    explicit SBend(const string &name);

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
  
    virtual bool apply( const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);
  
    virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual void rescaleFieldMap(const double &scaleFactor);

    virtual bool bends() const;

    void setAmplitudem(double vPeak);

    void setGapWidth(double gapwidth);

    /// Set the name of the "field map"; 
    /// For a bend this file contains the coefficients for the Enge function 
    void setFieldMapFN(string fmapfn);

    string getFieldMapFN() const;

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    double getR() const;

private:

    // Not implemented.
    void operator=(const SBend &);
    string filename_m;             /**< The name of the inputfile*/
    double amplitude_m;
    double B_xy_angle_m;
    double z_entry_m;
    double z_exit_m;
    double fringe_length_entry_m;
    double fringe_length_exit_m;
    double gap_width_m;
    double* enge_coefficients_m;
    unsigned int enge_polynomial_order_m;
    double enge_shift_entry_m;
    double enge_shift_exit_m;

    double R_m;
};

#endif // CLASSIC_SBend_HH
