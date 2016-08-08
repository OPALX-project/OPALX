#ifndef CLASSIC_BEND_H
#define CLASSIC_BEND_H

// ------------------------------------------------------------------------
// $RCSfile: SBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Bend
//   Defines the abstract interface for a general bend magnet.
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
#include "Fields/BMultipoleField.h"
#include "Algorithms/PartPusher.h"
#include "Physics/Physics.h"
#include <vector>


class Fieldmap;

/*
 * Class Bend
 *
 * Interface for general bend magnet.
 *
 * ------------------------------------------------------------------------
 *
 */

class Bend: public Component {

public:

    /// Constructor with given name.
    explicit Bend(const std::string &name);

    Bend();
    Bend(const Bend &);
    virtual ~Bend();

    /// Apply visitor to Bend.
    virtual void accept(BeamlineVisitor &) const = 0;

    /*
     * Methods for OPAL-SLICE.
     */
    virtual void addKR(int i, double t, Vector_t &K) { };
    virtual void addKT(int i, double t, Vector_t &K) { };


    /*
     * Methods for OPAL-T.
     * ===================
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
    virtual ElementBase::ElementType getType() const = 0;
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
    void SetAngleGreaterThanPiFlag(bool angleGreaterThanPi);
    void SetAperture(double aperture);
    virtual void SetBendAngle(double angle);
    void SetBeta(double beta);
    void SetDesignEnergy(double energy);
    virtual void SetEntranceAngle(double entranceAngle);
    void setExitAngle(double exitAngle);
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

protected:
    void setMessageHeader(const std::string & header);
    double getStartField() const;
    Fieldmap* getFieldmap();
    double getFieldAmplitude() const;

    double & getAngle();
    double & getEntranceAngle();
    double & getExitAngle();
    double getLength();
    double getMapLength();

private:

    // Not implemented.
    void operator=(const Bend &);

    void AdjustFringeFields(double ratio);
    double CalculateBendAngle();
    void CalcCentralField(const Vector_t &R,
                          double deltaX,
                          double angle,
                          Vector_t &B);
    void CalcDistFromRefOrbitCentralField(const Vector_t &R,
                                          double &deltaX,
                                          double &angle);
    void CalcEngeFunction(double zNormalized,
                          const std::vector<double> &engeCoeff,
                          int polyOrder,
                          double &engeFunc,
                          double &engeFuncDeriv,
                          double &engeFuncSecDerivNorm);
    void CalcEntranceFringeField(const Vector_t &REntrance,
                                 double deltaX,
                                 Vector_t &B);
    void CalcExitFringeField(const Vector_t &RExit, double deltaX, Vector_t &B);
    void CalculateMapField(const Vector_t &R, Vector_t &E, Vector_t &B);
    void CalculateRefTrajectory(double &angleX, double &angleY);
    double EstimateFieldAdjustmentStep(double actualBendAngle,
                                       double mass,
                                       double betaGamma);
    void FindBendEffectiveLength(double startField, double endField);
    void FindBendStrength(double mass,
                          double gamma,
                          double betaGamma,
                          double charge);
    virtual bool FindChordLength(Inform &msg,
                                 double &chordLength,
                                 bool &chordLengthFromMap) = 0;
    bool FindIdealBendParameters(double chordLength);
    void FindReferenceExitOrigin(double &x, double &z);
    bool InitializeFieldMap(Inform &msg);
    bool InMagnetCentralRegion(const Vector_t &R, double &deltaX, double &angle);
    bool InMagnetEntranceRegion(const Vector_t &R, double &deltaX);
    bool InMagnetExitRegion(const Vector_t &R, double &deltaX);
    bool IsPositionInEntranceField(const Vector_t &R, Vector_t &REntrance);
    bool IsPositionInExitField(const Vector_t &R, Vector_t &RExit);
    void Print(Inform &msg, double bendAngleX, double bendAngle);
    void ReadFieldMap(Inform &msg);
    bool Reinitialize();
    Vector_t RotateOutOfBendFrame(const Vector_t &X);
    Vector_t RotateToBendFrame(const Vector_t &X);
    void SetBendEffectiveLength(double startField, double endField);
    void SetBendStrength();
    void SetEngeOriginDelta(double delta);
    void SetFieldCalcParam(bool chordLengthFromMap);
    void SetGapFromFieldMap();
    bool SetupBendGeometry(Inform &msg, double &startField, double &endField);
    bool SetupDefaultFieldMap(Inform &msg);
    void SetFieldBoundaries(double startField, double endField);
    void SetupPusher(PartBunch *bunch);
    bool TreatAsDrift(Inform &msg, double chordlength);
    void retrieveDesignEnergy(double startField);

    std::string messageHeader_m;

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
    bool angleGreaterThanPi_m;  /// Set to true if bend angle is greater than
    /// 180 degrees.
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
     * magnet must be set in the input file.
     *
     * The magnet length is used initially to define the chord length of
     * the reference particle arc length. It is used at the start to
     * define the magnet's geometry. This parameter must be set in the
     * input file when using the default, internal field map. Otherwise
     * defining it is optional.
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
};

inline
void Bend::resetReinitializeFlag() {
    reinitialize_m = true;
}

inline
void Bend::resetRecalcRefTrajFlag() {
    recalcRefTraj_m = true;
}

inline
void Bend::setMessageHeader(const std::string & header)
{
    messageHeader_m = header;
}

inline
double Bend::getStartField() const
{
    return startField_m;
}

inline
Fieldmap* Bend::getFieldmap()
{
    return fieldmap_m;
}

inline
double Bend::getFieldAmplitude() const
{
    return fieldAmplitude_m;
}

inline
double & Bend::getAngle()
{
    return angle_m;
}

inline
double & Bend::getEntranceAngle()
{
    return entranceAngle_m;
}

inline
double & Bend::getExitAngle()
{
    return exitAngle_m;
}

inline
double Bend::getLength()
{
    return length_m;
}

inline
double Bend::getMapLength()
{
    return exitParameter2_m - entranceParameter2_m;
}

#endif // CLASSIC_SBend_HH