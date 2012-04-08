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

/*
 * Class SBend
 *
 * Interface for sector bend magnet.
 *
 * A sector bend magnet has a curved geometry. A sector magnet with zero degree edge angles is
 * simply a section of a circle when projected onto the y axis.
 *
 * The standard sector magnet, for purposes of definitions, has a field in the y direction.
 * This produces a bend in the horizontal (x) plane. Bends in other planes can be accomplished
 * by rotating the magnet about the axes.
 *
 * A positive bend angle is defined as one that bends a beam to the right when looking down
 * (in the negative y direction) so that the beam is bent in the negative x direction.
 *
 * A zero degree entrance edge angle is parallel to the x direction in an x/y/s coordinate system.
 * A positive entrance edge angle is defined as one that rotates the positive edge (in x) of the angle
 * toward the positive s axis.
 *
 * A zero degree exit edge angle is parallel to the x direction in an x/y/s coordinate system. A
 * positive exit edge angle is defined as one that rotates the positive edge (in x) of the angle toward
 * the negative s axis.
 *
 * ------------------------------------------------------------------------
 *
 * This class defines two interfaces:
 *
 * 1) Interface for sector magnets for OPAL-MAP.
 *
 *  Here we specify multipole components about the curved magnet trajectory.
 *
 *
 * 2) Interface for sector magnets for OPAL-T.
 *
 * Here we defined the magnet as a field map.
 */

class SBend: public Component {

public:

    /// Constructor with given name.
    explicit SBend(const std::string &name);

    SBend();
    SBend(const SBend &);
    virtual ~SBend();

    /// Apply visitor to SBend.
    virtual void accept(BeamlineVisitor &) const;


    /*
     * Methods for OPAL-MAP
     * ====================
     */

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

    void setFaceAngleEntry(double angle);
    void setFaceAngleExit(double angle);


    /*
     * Methods for OPAL-T.
     * ===================
     */

    /// Apply field to particles with coordinates in magnet frame.
    virtual bool apply(const int &i, const double &t, double E[], double B[]);
    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    /// Apply field to particles in beam frame.
    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    /// Setup bend.
    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    /// Indicates that element bends the beam.
    virtual bool bends() const;

    /// Get field map file name.
    std::string getFieldMapFN() const;

    /// Get element type.
    virtual const std::string &getType() const;

    /// Get bend dimensions. Includes fringe fields.
    virtual void getDimensions(double &zBegin, double &zEnd) const;

    double getBendAngle() const;
    double getEffectiveLength() const;
    double getEffectiveCenter() const;
    double getR() const;
    double getStartElement() const;

    void setAmplitudem(double vPeak);
    void setBendAngle(const double &);

    void setAlpha(const double &alpha);
    void setExitAngle(const double &exitAngle);
    void setBeta(const double &beta);
    void setDesignEnergy(const double &energy);

    /// Set quadrupole field component.
    void setK1(const double &k1);

    /*
     * Set the name of the field map.
     *
     * For now this means a file that contains Enge function coefficients
     * that describe the fringe fields at the entrance and exit.
     */
    void setFieldMapFN(std::string fmapfn);

    void setFullGap(const double &);
    void setLength(const double &);

    /// Set rotation about z axis in bend frame.
    void setLongitudinalRotation(const double &);
    void setLongitudinalRotation(const double &, const double &);

private:

    // Not implemented.
    void operator=(const SBend &);

    BorisPusher pusher_m;
    Vektor<double, 2> field_orientation_m;

    /// Name of field map that defines magnet.
    std::string filename_m;

    /// Magnet field map.
    Fieldmap *fieldmap_m;

    /// Flag to turn on fast field calculation. (Not currently used.)
    bool fast_m;

    /*
     * Flag to reinitialize the bend the first time the magnet
     * is called. This redefines the design energy of the bend
     * to the current average beam energy, keeping the bend angle
     * constant.
     */
    bool reinitialize_m;

    /// Start of magnet field map.
    double startField_m;

    /// End of magnet field map.
    double endField_m;

    /*
     * The magnet length and gap are used to define the magnet field
     * when the default, internal field map is used. Otherwise these
     * are effectively defined by the field map file.
     */
    double length_m;
    double gap_m;

    double ElementEdge_m;
    double startElement_m;

    /// Amplitude of magnet field (Tesla).
    double amplitude_m;

    /// Quadrupole component of field.
    double gradient_m;

    /// Bend angle for reference particle with bend design energy.
    double angle_m;
    double design_energy_m;

    /*
     * Edge angle parameters. We also store trig functions of these parameters
     * as they are used repeatedly.
     */
    double alpha_m; // Angle between incoming beam and the entrance face of the magnet.
    double exitAngle_m; // Angle between the outgoing, reference trajectory and exit face of the magnet.

    double sin_face_alpha_m; // alpha is the angle between the projection of the normal of the face onto the
    double cos_face_alpha_m; // s-u plane
    double tan_face_alpha_m;

    double sin_face_beta_m; // beta is defined as arctan(|n_parallel|/|n_perpendicular|) where n_parallel is the component
    double cos_face_beta_m; // of the normal of the face which is parallel to the s-u plane and n_perpendicular is
    double tan_face_beta_m; // perpendicular to it

    /// Map of reference particle trajectory.
    double *map_m;
    int map_size_m;
    double map_step_size_m;

    /*
     * Parameters that define an effective, hard edge bend that
     * is the approximate equivalent to the actual bend. Radius.
     * Effective length and effective center in s coordinates.
     * Effective start in floor coordinates.
     */
    double R_m;
    double effectiveLength_m;
    double effectiveCenter_m;
    double effectiveStart_m;

    double calculateBendAngle(double bendLength, bool modifyField);
    void calculateDistFromRef(Vector_t X, double &deltaX, double &angle);
    void calculateEffectiveLength();
    void calculateEffectiveCenter();
    void calculateMapField(Vector_t X, double &bX, double &bY, double &bZ);
    double calculateRefTrajectory(const double zBegin);
    bool reinitialize();
    void setBendStrength();

};

#endif // CLASSIC_SBend_HH
