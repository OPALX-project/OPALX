#ifndef OPALX_BendBase_HH
#define OPALX_BendBase_HH

#include "AbsBeamline/Component.h"
#include "Fields/BMultipoleField.h"

#include <cstddef>
#include <string>
#include <vector>

/**
 * @class BendBase
 * @brief Common OPALX interface for analytic horizontal bending magnets.
 *
 * This class provides the shared bookkeeping for sector and rectangular bends
 * in the first OPALX-native bend port. The placed hardware body is represented
 * by the element geometry, while the magnetic field is described by a local
 * multipole expansion evaluated in the element chart.
 *
 * For the analytic field model the dominant dipole contribution satisfies the
 * Lorentz-force design relation
 * \f[
 * \rho = \frac{\beta \gamma m}{|q| c |B|}
 * \f]
 * with design radius \f$\rho\f$, rest mass \f$m\f$, charge \f$q\f$, and
 * reference relativistic factor \f$\beta \gamma\f$.
 *
 * The first implementation keeps the nominal body extent and the field-support
 * extent identical. Fringe-field support will extend this model later without
 * changing the body-placement semantics.
 */
class BendBase : public Component {
public:
    BendBase();
    explicit BendBase(const std::string& name);
    BendBase(const BendBase&);
    ~BendBase() override;

    bool bends() const override;

    void initialise(PartBunch_t* bunch, double& startField, double& endField) override;
    void finalise() override;

    bool apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    bool apply(const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B)
            override;
    bool apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    void getFieldExtend(double& zBegin, double& zEnd) const override;
    bool isInside(const Vector_t<double, 3>& r) const override;

    CoordinateSystemTrafo getEdgeToBegin() const override;
    CoordinateSystemTrafo getEdgeToEnd() const override;

    /**
     * @brief Set the nominal body length.
     *
     * The bend body length is the placed hardware length and therefore drives
     * ports and visualization. The geometry itself remains the source of truth
     * for the placed body extent.
     */
    void setLength(double length);
    double getLength() const;

    /**
     * @brief Return the geometric chord between entrance and exit frames.
     *
     * For sector bends this is shorter than the arc length, while for
     * rectangular bends it equals the straight body length.
     */
    double getChordLength() const;

    virtual void setBendAngle(double angle);
    double getBendAngle() const;

    virtual void setEntranceAngle(double entranceAngle);
    double getEntranceAngle() const;

    virtual void setExitAngle(double exitAngle);
    virtual double getExitAngle() const = 0;

    void setFullGap(double gap);
    double getFullGap() const;

    void setDesignEnergy(const double& energy, bool changeable = true) override;
    double getDesignEnergy() const override;

    /**
     * @brief Sample the local curved reference path of the bend body.
     *
     * The returned points lie in the canonical local chart and can be mapped to
     * lab space via the placed body transform. The sampling is uniform in the
     * body parameter \f$s\f$ between the entrance and exit frames.
     */
    std::vector<Vector_t<double, 3>> getDesignPath(std::size_t minSamples = 32) const;

    /**
     * @brief Store the dipole design amplitudes used by the analytic field.
     *
     * The normal and skew dipole strengths are the leading coefficients of the
     * local multipole expansion
     * \f[
     * B_y + i B_x = \sum_{n=0}^{N} (B_n + i A_n) (x + i y)^n .
     * \f]
     */
    void setFieldAmplitude(double k0, double k0s);
    double getFieldAmplitude() const;

    void setFieldMapFN(std::string fileName);
    std::string getFieldMapFN() const;

    double getB() const;
    void setB(double B);

    void setEntryFaceRotation(double rotation);
    double getEntryFaceRotation() const;

    void setExitFaceRotation(double rotation);
    double getExitFaceRotation() const;

    void setEntryFaceCurvature(double curvature);
    double getEntryFaceCurvature() const;

    void setExitFaceCurvature(double curvature);
    double getExitFaceCurvature() const;

    void setSlices(double slices);
    double getSlices() const;

    void setStepsize(double stepSize);
    double getStepsize() const;

    void setNSlices(const std::size_t& nSlices);
    std::size_t getNSlices() const;

    void setK1(double k1);
    double getK1() const;

    int getRequiredNumberOfTimeSteps() const override;

    virtual BMultipoleField& getField() override             = 0;
    virtual const BMultipoleField& getField() const override = 0;

protected:
    double calcDesignRadius(double fieldAmplitude) const;
    double calcFieldAmplitude(double radius) const;
    double calcBendAngle(double chordLength, double radius) const;
    double calcDesignRadius(double chordLength, double angle) const;
    double calcGamma() const;
    double calcBetaGamma() const;
    double getStoredExitAngle() const;

private:
    static void computeFieldHost(
            const Vector_t<double, 3>& R, const BMultipoleField& field, Vector_t<double, 3>& B);

    double startField_m;
    double endField_m;
    double angle_m;
    double entranceAngle_m;
    double exitAngle_m;
    double gap_m;
    double designEnergy_m;
    bool designEnergyChangeable_m;
    double fieldAmplitudeX_m;
    double fieldAmplitudeY_m;
    double fieldAmplitude_m;
    std::string fileName_m;
    double entryFaceRotation_m;
    double exitFaceRotation_m;
    double entryFaceCurvature_m;
    double exitFaceCurvature_m;
    double slices_m;
    double stepSize_m;
    std::size_t nSlices_m;
    double k1_m;
};

inline bool BendBase::bends() const { return true; }

inline void BendBase::setLength(double length) { setElementLength(std::abs(length)); }

inline double BendBase::getLength() const { return getElementLength(); }

inline void BendBase::setBendAngle(double angle) { angle_m = angle; }

inline double BendBase::getBendAngle() const { return angle_m; }

inline void BendBase::setEntranceAngle(double entranceAngle) { entranceAngle_m = entranceAngle; }

inline double BendBase::getEntranceAngle() const { return entranceAngle_m; }

inline void BendBase::setExitAngle(double exitAngle) { exitAngle_m = exitAngle; }

inline void BendBase::setFullGap(double gap) { gap_m = std::abs(gap); }

inline double BendBase::getFullGap() const { return gap_m; }

inline void BendBase::setDesignEnergy(const double& energy, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m           = std::abs(energy) * 1e6;
        designEnergyChangeable_m = changeable;
    }
}

inline double BendBase::getDesignEnergy() const { return designEnergy_m; }

inline void BendBase::setFieldAmplitude(double k0, double k0s) {
    fieldAmplitudeY_m = k0;
    fieldAmplitudeX_m = k0s;
    fieldAmplitude_m  = std::hypot(k0, k0s);
}

inline double BendBase::getFieldAmplitude() const { return fieldAmplitude_m; }

inline void BendBase::setFieldMapFN(std::string fileName) { fileName_m = std::move(fileName); }

inline std::string BendBase::getFieldMapFN() const { return fileName_m; }

inline double BendBase::getB() const { return getField().getNormalComponent(0); }

inline void BendBase::setB(double B) { getField().setNormalComponent(0, B); }

inline void BendBase::setEntryFaceRotation(double rotation) { entryFaceRotation_m = rotation; }

inline double BendBase::getEntryFaceRotation() const { return entryFaceRotation_m; }

inline void BendBase::setExitFaceRotation(double rotation) { exitFaceRotation_m = rotation; }

inline double BendBase::getExitFaceRotation() const { return exitFaceRotation_m; }

inline void BendBase::setEntryFaceCurvature(double curvature) { entryFaceCurvature_m = curvature; }

inline double BendBase::getEntryFaceCurvature() const { return entryFaceCurvature_m; }

inline void BendBase::setExitFaceCurvature(double curvature) { exitFaceCurvature_m = curvature; }

inline double BendBase::getExitFaceCurvature() const { return exitFaceCurvature_m; }

inline void BendBase::setSlices(double slices) { slices_m = slices; }

inline double BendBase::getSlices() const { return slices_m; }

inline void BendBase::setStepsize(double stepSize) { stepSize_m = stepSize; }

inline double BendBase::getStepsize() const { return stepSize_m; }

inline void BendBase::setNSlices(const std::size_t& nSlices) { nSlices_m = nSlices; }

inline std::size_t BendBase::getNSlices() const { return nSlices_m; }

inline void BendBase::setK1(double k1) { k1_m = k1; }

inline double BendBase::getK1() const { return k1_m; }

inline int BendBase::getRequiredNumberOfTimeSteps() const { return 10; }

inline double BendBase::getStoredExitAngle() const { return exitAngle_m; }

#endif  // OPALX_BendBase_HH
