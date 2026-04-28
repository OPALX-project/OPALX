#ifndef OPALX_BendBase_HH
#define OPALX_BendBase_HH

#include "AbsBeamline/Component.h"
#include "Fields/BMultipoleField.h"

#include <cmath>
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
 * The bend body extent remains the nominal placed hardware extent \f$[0, L]\f$
 * in the local chart, while the analytic field support may extend beyond the
 * body through hard-edge fringe regions
 * \f[
 * z \in [-\ell_\mathrm{entry},\, L + \ell_\mathrm{exit}] .
 * \f]
 * The fringe support lengths are modeled from the half gap \f$h\f$, fringe
 * integral \f$F_\mathrm{int}\f$, and pole-face angles \f$E_1, E_2\f$ as
 * \f[
 * \ell_\mathrm{entry} = \frac{h F_\mathrm{int}}{|\cos E_1|}, \qquad
 * \ell_\mathrm{exit} = \frac{h F_\mathrm{int}}{|\cos E_2|}.
 * \f]
 * The corresponding effective integrated field length of the linear-ramp
 * hard-edge model is
 * \f[
 * L_\mathrm{eff} = L_\mathrm{ref}
 *                + \frac{\ell_\mathrm{entry} + \ell_\mathrm{exit}}{2},
 * \f]
 * where \f$L_\mathrm{ref}\f$ is the body length measured along the design
 * reference path. For sector bends this equals the placed arc length, while
 * for rectangular bends it is the reference-path length between the rotated
 * entrance and exit faces rather than the straight body chord.
 *
 * The first optics-level fringe model starts from the standard hard-edge
 * edge-focusing relations
 * \f[
 * x'_\mathrm{out} = x'_\mathrm{in} + h \tan(E)\,x, \qquad
 * y'_\mathrm{out} = y'_\mathrm{in} - h \tan(E)\,y,
 * \f]
 * with signed curvature \f$h = \theta / L_\mathrm{eff}\f$. The fringe
 * integral modifies the vertical edge focusing through the standard effective
 * edge-angle correction
 * \f[
 * \psi = h\,h_\mathrm{gap}\,F_\mathrm{int}\,
 *        \frac{1 + \sin^2(E)}{\cos(E)},
 * \qquad
 * y'_\mathrm{out} = y'_\mathrm{in} - h \tan(E - \psi)\,y.
 * \f]
 * OPALX realizes these first-order optics effects through an equivalent
 * distributed fringe-field correction across the entry and exit support
 * regions while keeping the bend body-placement semantics unchanged.
 *
 * The analytic multipole strengths are stored internally in the normalized
 * lattice convention \f$k_n\f$ and converted to physical Tesla coefficients
 * only when a reference rigidity is available. If `DESIGNENERGY` is provided,
 * the runtime scaling uses that kinetic energy together with the bunch species
 * mass and charge. Otherwise the scaling uses the bunch reference momentum
 * directly. This avoids coupling bend tracking strength to the parser-global
 * variable \f$P_0\f$.
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

    /**
     * @brief Set the half gap used for analytic fringe support.
     *
     * The half gap \f$h\f$ determines the extent of the hard-edge fringe
     * support together with the fringe integral and pole-face angles.
     */
    void setFringeHalfGap(double halfGap);
    double getFringeHalfGap() const;

    /**
     * @brief Set the hard-edge fringe integral used for support extension.
     *
     * The analytic fringe model uses the thin-lens hard-edge convention
     * \f[
     * \ell = \frac{h F_\mathrm{int}}{|\cos E|}.
     * \f]
     */
    void setFringeIntegral(double fringeIntegral);
    double getFringeIntegral() const;

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

    /**
     * @brief Store the normalized analytic multipole strengths.
     *
     * The coefficients are stored in lattice convention and later converted to
     * physical Tesla strengths from either the explicit design energy or the
     * runtime bunch reference momentum.
     */
    void setNormalizedField(const BMultipoleField& field);

    /**
     * @brief Resolve the reference momentum used for runtime field scaling.
     *
     * If `DESIGNENERGY` is set, the bend uses the corresponding momentum for
     * the given particle mass. Otherwise it uses the bunch reference momentum.
     *
     * @param bunchReferenceMomentumEV Bunch reference momentum in eV/c.
     * @param massEV Particle rest mass in eV.
     * @return Selected reference momentum in eV/c.
     */
    double resolveReferenceMomentumEV(double bunchReferenceMomentumEV, double massEV) const;

    /**
     * @brief Convert the stored normalized multipole strengths to Tesla units.
     *
     * The supplied momentum is interpreted as the reference momentum in eV/c.
     * The physical multipole coefficients are obtained from
     * \f[
     * B_n = \frac{p_0}{q c} k_n.
     * \f]
     */
    void updatePhysicalFieldFromMomentumEV(double referenceMomentumEV, double charge);

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

    /**
     * @brief Rigid tracking slice used by many-particle bend tracking.
     *
     * The slice coordinate system is attached to the design reference path at
     * longitudinal coordinate \f$s_\mathrm{center}\f$ measured from the bend
     * entry face. The stored transform maps entry-frame coordinates into the
     * slice-local tangent frame. The slice interval itself is
     * \f$[s_\mathrm{begin}, s_\mathrm{end}]\f$ on the same entry-based
     * reference-path coordinate.
     */
    struct TrackingSlice {
        CoordinateSystemTrafo entryToSliceLocal;
        double sBegin;
        double sCenter;
        double sEnd;
    };

    /**
     * @brief Build the rigid tangent-frame slices used for many-particle tracking.
     *
     * The slice decomposition is defined on the entry-based design-path
     * coordinate and covers the full field-support extent
     * \f$[-\ell_\mathrm{entry}, L_\mathrm{ref} + \ell_\mathrm{exit}]\f$.
     * The slices are uniform in \f$s\f$ and are used only for the many-particle
     * tracking path. The reference-particle path remains field-based and
     * unsliced.
     */
    std::vector<TrackingSlice> buildTrackingSlices() const;

    /**
     * @brief Apply the analytic bend field to particles already expressed in one slice frame.
     *
     * The particles are assumed to be in the slice-local tangent frame
     * associated with `slice`. The kernel evaluates the existing analytic bend
     * multipole field at the approximated global reference-path coordinate
     * \f$s = s_\mathrm{center} + z_\mathrm{local}\f$ and accumulates the field
     * contribution for particles inside the slice support.
     *
     * @note This path is used only for many-particle tracking. It does not
     * alter the reference-particle field evaluation.
     */
    bool applySlice(
            const std::shared_ptr<ParticleContainer_t>& pc, const TrackingSlice& slice) const;

    /**
     * @brief Return the local hard-edge fringe support in front of the body.
     */
    double getEntryFringeSupportLength() const;

    /**
     * @brief Return the local hard-edge fringe support behind the body.
     */
    double getExitFringeSupportLength() const;

    /**
     * @brief Return the effective integrated field length of the fringe model.
     *
     * The effective length is defined on the design reference path used to
     * normalize the analytic dipole field. For sector bends this is the arc
     * length of the body. For rectangular bends it is the design-path length
     * between the rotated entrance and exit faces, which differs from the
     * straight body length.
     */
    double getEffectiveFieldLength() const;

    /**
     * @brief Return the local field scale for the hard-edge fringe model.
     *
     * The scale is zero outside the field support, ramps linearly from zero to
     * one across the entry fringe, remains one on the body interval, and ramps
     * linearly back to zero across the exit fringe.
     */
    double getFieldScale(double z) const;

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
    double getSignedCurvature() const;
    double getEntryEdgeHorizontalStrength() const;
    double getExitEdgeHorizontalStrength() const;
    double getEntryEdgeVerticalStrength() const;
    double getExitEdgeVerticalStrength() const;
    void updateFieldSupportExtent();
    void updatePhysicalFieldFromReference();
    /**
     * @brief Return the hard-edge design-path body length used for curvature.
     *
     * This is the longitudinal extent of the bend body measured along the
     * reference path used to define the nominal deflection angle. The default
     * implementation uses the nominal body length itself; `RBEND` overrides
     * this with its arc-equivalent design-path length through the rotated
     * faces.
     */
    virtual double getReferencePathLength() const;

private:
    static void computeFieldHost(
            const Vector_t<double, 3>& R, const BMultipoleField& field, Vector_t<double, 3>& B);
    double getReferenceMomentumEV() const;

    double startField_m;
    double endField_m;
    double fieldBegin_m;
    double fieldEnd_m;
    double angle_m;
    double entranceAngle_m;
    double exitAngle_m;
    double gap_m;
    double fringeHalfGap_m;
    double fringeIntegral_m;
    double designEnergy_m;
    bool designEnergyChangeable_m;
    double fieldAmplitudeX_m;
    double fieldAmplitudeY_m;
    double fieldAmplitude_m;
    BMultipoleField normalizedField_m;
    bool hasNormalizedField_m;
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

inline void BendBase::setLength(double length) {
    setElementLength(std::abs(length));
    updateFieldSupportExtent();
}

inline double BendBase::getLength() const { return getElementLength(); }

inline void BendBase::setBendAngle(double angle) {
    angle_m = angle;
    updateFieldSupportExtent();
}

inline double BendBase::getBendAngle() const { return angle_m; }

inline void BendBase::setEntranceAngle(double entranceAngle) {
    entranceAngle_m = entranceAngle;
    updateFieldSupportExtent();
}

inline double BendBase::getEntranceAngle() const { return entranceAngle_m; }

inline void BendBase::setExitAngle(double exitAngle) {
    exitAngle_m = exitAngle;
    updateFieldSupportExtent();
}

inline void BendBase::setFullGap(double gap) { gap_m = std::abs(gap); }

inline double BendBase::getFullGap() const { return gap_m; }

inline void BendBase::setFringeHalfGap(double halfGap) {
    fringeHalfGap_m = std::abs(halfGap);
    updateFieldSupportExtent();
}

inline double BendBase::getFringeHalfGap() const { return fringeHalfGap_m; }

inline void BendBase::setFringeIntegral(double fringeIntegral) {
    fringeIntegral_m = std::abs(fringeIntegral);
    updateFieldSupportExtent();
}

inline double BendBase::getFringeIntegral() const { return fringeIntegral_m; }

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

inline void BendBase::setNormalizedField(const BMultipoleField& field) {
    normalizedField_m    = field;
    hasNormalizedField_m = true;
}

inline double BendBase::resolveReferenceMomentumEV(
        const double bunchReferenceMomentumEV, const double massEV) const {
    if (designEnergy_m > 0.0) {
        const double gamma = designEnergy_m / massEV + 1.0;
        return std::sqrt(gamma * gamma - 1.0) * massEV;
    }

    return bunchReferenceMomentumEV;
}

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

inline double BendBase::getEntryFringeSupportLength() const { return -fieldBegin_m; }

inline double BendBase::getExitFringeSupportLength() const {
    return std::max(0.0, fieldEnd_m - getElementLength());
}

inline double BendBase::getReferencePathLength() const { return getElementLength(); }

inline double BendBase::getEffectiveFieldLength() const {
    return getReferencePathLength()
           + 0.5 * (getEntryFringeSupportLength() + getExitFringeSupportLength());
}

inline double BendBase::getStoredExitAngle() const { return exitAngle_m; }

#endif  // OPALX_BendBase_HH
