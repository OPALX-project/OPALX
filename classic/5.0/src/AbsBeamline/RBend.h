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

/*
 * Class RBend
 *
 * Interface for a rectangular bend magnet.
 *
 * A rectangular bend magnet physically has a rectangular shape.
 *
 * The standard rectangular magnet, for purposes of definitions, has a field in
 * the y direction. This produces a bend in the horizontal (x) plane. Bends in
 * other planes can be accomplished by rotating the magnet about the axes.
 *
 * A positive bend angle is defined as one that bends a beam to the right when
 * looking down (in the negative y direction) so that the beam is bent in the
 * negative x direction. (This definition of a positive bend is the same whether
 * the charge is positive or negative.)
 *
 * A zero degree entrance edge angle is parallel to the x direction in an x/y/s
 * coordinate system. A positive entrance edge angle is defined as one that
 * rotates the positive edge (in x) of the angle toward the positive s axis.
 *
 * Since the magnet geometry is a fixed rectangle, the exit edge angle is
 * defined by the bend angle of the magnet and the entrance edge angle. In
 * general, the exit edge angle is equal to the bend angle minus the entrance
 * edge angle.
 *
 * ------------------------------------------------------------------------
 *
 * This class defines three interfaces:
 *
 * 1) Interface for rectangular magnets for OPAL-MAP.
 *
 *  Here we specify multipole components about the curved magnet trajectory.
 *
 *
 * 2) Interface for rectangular magnets in OPAL-SLICE.
 *
 * Here we define a rectangular magnet for use in the OPAL
 * slice model.
 *
 *
 * 3) Interface for rectangular magnets for OPAL-T.
 *
 * Here we defined the magnet as a field map.
 */

class RBend: public Component {

public:

    /// Constructor with given name.
    explicit RBend(const std::string &name);

    RBend();
    RBend(const RBend &);
    virtual ~RBend();

    /// Apply visitor to RBend.
    virtual void accept(BeamlineVisitor &) const;

    /*
     * Methods for OPAL-MAP
     * ====================
     */

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


    /*
     * Methods for OPAL-SLICE.
     */
    virtual void addKR(int i, double t, Vector_t &K);
    virtual void addKT(int i, double t, Vector_t &K);


    /*
     * Methods for OPAL-T.
     */

    /// Apply field to particles with coordinates in magnet frame.
    virtual bool apply(const size_t &i,
                       const double &t,
                       double E[],
                       double B[]);
    virtual bool apply(const size_t &i,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B);

    /// Apply field to particles in beam frame.
    virtual bool apply(const Vector_t &R,
                       const Vector_t &centroid,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B);

    /// Indicates that element bends the beam.
    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void finalise();
    virtual void getDimensions(double &sBegin, double &sEnd) const;
    virtual const std::string &getType() const;
    virtual void initialise(PartBunch *bunch,
                            double &startField,
                            double &endField,
                            const double &scaleFactor);

    double GetBendAngle() const;
    double GetBendRadius() const;
    double GetEffectiveCenter() const;
    double GetEffectiveLength() const;
    std::string GetFieldMapFN() const;
    double GetStartElement() const;
    void SetAperture(double aperture);
    void SetBendAngle(double angle);
    void SetBeta(double beta);
    void SetDesignEnergy(double energy);
    void SetEntranceAngle(double entranceAngle);
    void SetFieldAmplitude(double k0, double k0s);

    /*
     * Set the name of the field map.
     *
     * For now this means a file that contains Enge function coefficients
     * that describe the fringe fields at the entrance and exit.
     */
    void SetFieldMapFN(std::string fileName);
    void SetFullGap(double gap);

    /// Set quadrupole field component.
    void SetK1(double k1);
    void SetLength(double length);

    /// Set rotation about z axis in bend frame.
    void SetRotationAboutZ(double rotation);

    void resetReinitializeFlag();
    void resetRecalcRefTrajFlag();

private:

    // Not implemented.
    void operator=(const RBend &);

    void AdjustFringeFields(double ratio);
    double CalculateBendAngle();
    void CalcCentralField(Vector_t R,
                          double deltaX,
                          double angle,
                          Vector_t &B);
    void CalcDistFromRefOrbitCentralField(Vector_t R,
                                          double &deltaX,
                                          double &angle);
    void CalcEngeFunction(double zNormalized,
                          std::vector<double> engeCoeff,
                          int polyOrder,
                          double &engeFunc,
                          double &engeFuncDeriv,
                          double &engeFuncSecDerivNorm);
    void CalcEntranceFringeField(Vector_t REntrance,
                                 double deltaX,
                                 Vector_t &B);
    void CalcExitFringeField(Vector_t RExit, double deltaX, Vector_t &B);
    void CalculateMapField(Vector_t R, Vector_t &E, Vector_t &B);
    void CalculateRefTrajectory(double &angleX, double &angleY);
    double EstimateFieldAdjustmentStep(double actualBendAngle,
                                       double mass,
                                       double betaGamma);
    void FindBendEffectiveLength(double startField, double endField);
    bool FindBendLength(Inform &msg,
                        double &bendLength,
                        bool &bendLengthFromMap);
    void FindBendStrength(double mass,
                          double gamma,
                          double betaGamma,
                          double charge);
    bool FindIdealBendParameters(double bendLength);
    void FindReferenceExitOrigin(double &x, double &z);
    bool InitializeFieldMap(Inform &msg);
    bool InMagnetCentralRegion(Vector_t R, double &deltaX, double &angle);
    bool InMagnetEntranceRegion(Vector_t R, double &deltaX);
    bool InMagnetExitRegion(Vector_t R, double &deltaX);
    bool IsPositionInEntranceField(Vector_t R, Vector_t &REntrance);
    bool IsPositionInExitField(Vector_t R, Vector_t &RExit);
    void Print(Inform &msg, double bendAngleX, double bendAngle);
    void ReadFieldMap(Inform &msg);
    bool Reinitialize();
    Vector_t RotateOutOfBendFrame(Vector_t X);
    Vector_t RotateToBendFrame(Vector_t X);
    void SetBendEffectiveLength(double startField, double endField);
    void SetBendStrength();
    void SetEngeOriginDelta(double delta);
    void SetFieldCalcParam(bool lengthFromMap);
    void SetGapFromFieldMap();
    bool SetupBendGeometry(Inform &msg, double &startField, double &endField);
    bool SetupDefaultFieldMap(Inform &msg);
    void SetFieldBoundaries(double startField, double endField);
    void SetupPusher(PartBunch *bunch);
    bool TreatAsDrift(Inform &msg);

    BorisPusher pusher_m;       /// Pusher used to integrate reference particle
    /// through the bend.
    std::string fileName_m;     /// Name of field map that defines magnet.
    Fieldmap *fieldmap_m;       /// Magnet field map.
    bool fast_m;                /// Flag to turn on fast field calculation.
    /// (Not currently used.)
    double angle_m;             /// Bend angle for reference particle with bend
    /// design energy (radians).
    double aperture_m;          /// Aperture of magnet in non-bend plane.
    double designEnergy_m;      /// Bend design energy (eV).
    double designRadius_m;      /// Bend design radius (m).
    double fieldAmplitude_m;    /// Amplitude of magnet field (T).
    double bX_m;                /// Amplitude of Bx field (T).
    double bY_m;                /// Amplitude of By field (T).
    double entranceAngle_m;     /// Angle between incoming reference trajectory
    /// and the entrance face of the magnet (radians).
    double exitAngle_m;         /// Angle between outgoing reference trajectory
    /// and the exit face of the magnet (radians).
    double fieldIndex_m;        /// Dipole field index.
    double elementEdge_m;       /// Physical start of magnet in s coordinates (m).
    double startField_m;        /// Start of magnet field map in s coordinates (m).
    double endField_m;          /// End of magnet field map in s coordinates (m).

    /*
     * Flag to reinitialize the bend the first time the magnet
     * is called. This redefines the design energy of the bend
     * to the current average beam energy, keeping the bend angle
     * constant.
     */
    bool reinitialize_m;

    bool recalcRefTraj_m;       /// Re-calculate reference trajectory.

    /*
     * These two parameters are used to set up the bend geometry.
     *
     * When using the default, internal field map, the gap of the
     * magnet must be set in the input file. Otherwise, this parameter
     * is defined by the magnet field map.
     *
     * The magnet length is used to define the length of the magnet.
     * This parameter must be set in the input file when using the default,
     * internal field map. Otherwise defining it is optional.
     */
    double length_m;
    double gap_m;

    /// Map of reference particle trajectory.
    std::vector<double> refTrajMapX_m;
    std::vector<double> refTrajMapY_m;
    std::vector<double> refTrajMapZ_m;
    int refTrajMapSize_m;
    double refTrajMapStepSize_m;

    /*
     * Enge function field map members.
     */

    /*
     * Entrance and exit position parameters. Ultimately they are used to
     * determine the origins of the entrance and exit edge Enge functions and
     * the extent of the field map. However, how they are used to do this
     * depends on how the bend using the map is setup in the OPAL input file.
     * So, we use generic terms to start.
     */
    double entranceParameter1_m;
    double entranceParameter2_m;
    double entranceParameter3_m;
    double exitParameter1_m;
    double exitParameter2_m;
    double exitParameter3_m;

    /// Enge coefficients for map entry and exit regions.
    std::vector<double> engeCoeffsEntry_m;
    std::vector<double> engeCoeffsExit_m;

    /*
     * All coordinates are with respect to (x, z) = (0, 0). It is
     * assumed the ideal reference trajectory passes through this point.
     */
    double xOriginEngeEntry_m;      /// x coordinate of entry Enge function origin.
    double zOriginEngeEntry_m;      /// z coordinate of entry Enge function origin.
    double deltaBeginEntry_m;       /// Perpendicular distance from entrance Enge
    /// function origin where Enge function starts.
    double deltaEndEntry_m;         /// Perpendicular distance from entrance Enge
    /// function origin that Enge function ends.
    int polyOrderEntry_m;           /// Enge function order for entry region.

    /*
     * The ideal reference trajectory passes through (xExit_m, zExit_m).
     */
    double xExit_m;
    double zExit_m;

    double xOriginEngeExit_m;       /// x coordinate of exit Enge function origin.
    double zOriginEngeExit_m;       /// z coordinate of exit Enge function origin.
    double deltaBeginExit_m;        /// Perpendicular distance from exit Enge
    /// function origin that Enge function starts.
    double deltaEndExit_m;          /// Perpendicular distance from exit Enge
    /// function origin that Enge function ends.
    int polyOrderExit_m;            /// Enge function order for entry region.

    double cosEntranceAngle_m;
    double sinEntranceAngle_m;
    double exitEdgeAngle_m;         /// Exit edge angle (radians.
    double cosExitAngle_m;
    double sinExitAngle_m;

    // For OPAL-SLICE.
    static int RBend_counter_m;
};

inline
void RBend::resetReinitializeFlag() {
    reinitialize_m = true;
}

inline
void RBend::resetRecalcRefTrajFlag() {
    recalcRefTraj_m = true;
}

#endif // CLASSIC_RBend_HH