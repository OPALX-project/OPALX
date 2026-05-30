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
 * body through the same default Enge fringe profile used by OPAL's
 * `Bend2D`/`FM1DProfile1` model for `1DPROFILE1-DEFAULT`:
 * \f[
 * z \in [-\ell_\mathrm{entry},\, L + \ell_\mathrm{exit}] .
 * \f]
 * OPAL's default map defines field-profile points
 * \f$p_1=-0.1\,\mathrm{m}\f$, \f$p_2=0\f$, and
 * \f$p_3=0.1\,\mathrm{m}\f$ for a
 * \f$g_0=0.02\,\mathrm{m}\f$ full gap. OPALX scales these points by the
 * configured full gap \f$g\f$ (`GAP`, or \f$2\,\mathrm{HGAP}\f$ if `GAP` is
 * absent), so the physical half width of one fringe is
 * \f[
 * a = 0.1\,\mathrm{m}\,\frac{g}{g_0} = 5g .
 * \f]
 * The corresponding local support lengths are projected through the pole-face
 * angles,
 * \f[
 * \ell_\mathrm{entry} = \frac{a}{|\cos E_1|}, \qquad
 * \ell_\mathrm{exit} = \frac{a}{|\cos E_2|}.
 * \f]
 *
 * In the fringe regions the longitudinal field envelope is the OPAL default
 * Enge function
 * \f[
 * f(u) = \frac{1}{1 + \exp\!\left(\sum_{i=0}^{5} c_i u^i\right)},
 * \f]
 * with
 * \f[
 * (c_0,\ldots,c_5) =
 * (0.478959,\;1.911289,\;-1.185953,\;1.630554,\;-1.082657,\;0.318111).
 * \f]
 * If \f$\chi_1=|\cos E_1|\f$ and \f$\chi_2=|\cos E_2|\f$, the entry and exit
 * envelopes are
 * \f[
 * F_\mathrm{entry}(z)=f\!\left(-\frac{\chi_1 z}{g}\right), \qquad
 * F_\mathrm{exit}(z)=f\!\left(\frac{\chi_2(z-L_\mathrm{ref})}{g}\right).
 * \f]
 * OPALX uses the smaller of the applicable entry/exit envelopes for short
 * bends where the fringe regions overlap. The on-axis dipole field is
 * multiplied by \f$F(z)\f$ and the leading source-free fringe correction
 * \f[
 * B_y(x,y,z) = B_0\!\left(F(z) - \frac{1}{2}F''(z)y^2\right), \qquad
 * B_z(x,y,z) = B_0 F'(z)y
 * \f]
 * is applied to the normal dipole term.
 *
 * The corresponding effective integrated field length is now computed from the
 * Enge envelope itself,
 * \f[
 * L_\mathrm{eff} = \int_{-\ell_\mathrm{entry}}^{L_\mathrm{ref}
 *                   + \ell_\mathrm{exit}} F(z)\,dz ,
 * \f]
 * where \f$L_\mathrm{ref}\f$ is the body length measured along the design
 * reference path. For sector bends this equals the placed arc length, while
 * for rectangular bends it is the reference-path length between the rotated
 * entrance and exit faces rather than the straight body chord.
 *
 * The distributed fringe model follows OPAL's `Bend2D` field shape: no
 * additional hard-edge \f$B_y \propto x\f$ term is superposed on the Enge
 * field. The optics-level correction retained here is the vertical edge
 * focusing term. With signed curvature \f$h=\theta/L_\mathrm{eff}\f$, the
 * fringe integral modifies the vertical edge focusing through the standard
 * effective edge-angle correction
 * \f[
 * \psi = h\,h_\mathrm{gap}\,F_\mathrm{int}\,
 *        \frac{1 + \sin^2(E)}{\cos(E)},
 * \qquad
 * y'_\mathrm{out} = y'_\mathrm{in} - h \tan(E - \psi)\,y.
 * \f]
 * The `FINT` parameter remains part of this first-order pole-face focusing
 * correction. It does not set the Enge field-profile width in OPAL's
 * `FM1DProfile1` model, and it does not introduce a separate horizontal field
 * gradient in the OPAL-compatible map.
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
     * If no full `GAP` is provided, the OPAL default Enge fringe model uses
     * \f$g=2h\f$ as the full gap when scaling the `FM1DProfile1` profile
     * points.
     */
    void setFringeHalfGap(double halfGap);
    double getFringeHalfGap() const;

    /**
     * @brief Set the pole-face fringe integral used by the edge-focusing model.
     *
     * `FINT` enters the first-order vertical edge-angle correction
     * \f[
     * \psi = h_\mathrm{eff}\,h_\mathrm{gap}\,F_\mathrm{int}
     *        \frac{1 + \sin^2 E}{\cos E}.
     * \f]
     * It does not define the longitudinal Enge field-support width; that width
     * follows OPAL's `FM1DProfile1` default profile and the configured gap.
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

    /**
     * @brief Set the legacy bend slice-count hint.
     *
     * This scalar is preserved for parser compatibility. When no explicit
     * `STEPSIZE` or integral `NSLICES` is provided, the sliced many-particle
     * bend tracker interprets values larger than one as a requested number of
     * rigid body slices after rounding upward.
     */
    void setSlices(double slices);

    /**
     * @brief Return the legacy bend slice-count hint.
     */
    double getSlices() const;

    /**
     * @brief Set the target path-length step used to build bend body slices.
     *
     * If positive, the many-particle bend tracker constructs
     * \f$N = \lceil L_\mathrm{support} / \Delta s \rceil\f$ rigid slices over
     * the full bend field-support interval, where \f$L_\mathrm{support}\f$
     * denotes the entry-fringe, body, and exit-fringe support combined.
     */
    void setStepsize(double stepSize);

    /**
     * @brief Return the target path-length step used for bend slice construction.
     */
    double getStepsize() const;

    /**
     * @brief Set the explicit number of rigid body slices for many-particle tracking.
     *
     * This takes precedence over the legacy floating-point `SLICES` hint but
     * is overridden by a positive `STEPSIZE`, which defines the slice count
     * geometrically from the bend support span.
     */
    void setNSlices(const std::size_t& nSlices);

    /**
     * @brief Return the explicit number of rigid body slices.
     */
    std::size_t getNSlices() const;

    /**
     * @brief Set the normalized normal quadrupole component \f$k_1\f$.
     *
     * For analytic bends this is stored alongside the dipole term in the local
     * multipole expansion and contributes only on the bend body slices unless a
     * separate thin-edge kick is applied at the faces.
     */
    void setK1(double k1);

    /**
     * @brief Return the stored normalized normal quadrupole component \f$k_1\f$.
     */
    double getK1() const;

    int getRequiredNumberOfTimeSteps() const override;

    /**
     * @brief Rigid tracking slice used by many-particle bend tracking.
     *
     * The slice coordinate system is attached to the design reference path at
     * longitudinal coordinate \f$s_\mathrm{center}\f$ measured from the bend
     * entry face. The stored transforms map either entry-frame coordinates or
     * the explicit body frame used by placed sector bends into the slice-local
     * tangent frame. The slice interval itself is
     * \f$[s_\mathrm{begin}, s_\mathrm{end}]\f$ on the same entry-based
     * reference-path coordinate.
     */
    struct TrackingSlice {
        CoordinateSystemTrafo entryToSliceLocal;
        CoordinateSystemTrafo bodyToSliceLocal;
        matrix3x3_t entryToSliceRotation;
        Vector_t<double, 3> entryOrigin;
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
     * @brief Convert a rigid entry-frame Cartesian point into the bend field chart.
     *
     * The analytic bend field is evaluated in the piecewise entry-based chart
     * \f$(x, y, s)\f$, where \f$s\f$ follows the design reference path. This
     * helper maps a point from the rigid entry frame into that chart, including
     * the straight entry and exit tangent extensions outside the body interval.
     */
    Vector_t<double, 3> convertEntryCartesianToFieldLocal(
            const Vector_t<double, 3>& entryCartesian) const;

    /**
     * @brief Convert a midpoint-body Cartesian point into the bend field chart.
     *
     * Explicitly placed bends are rendered and positioned through a rigid body
     * frame, while the analytic field is evaluated in the entry-based curved
     * reference chart. Sector bends use the reference-path midpoint transform.
     * Rectangular bends use their actual pole-face transforms because the
     * straight hardware body is not the same chart as the circular reference
     * orbit. The entrance fringe is measured in the entrance-face chart and the
     * exit fringe is measured in the exit-face chart, matching OPAL `Bend2D`:
     * \f[
     * \mathbf{r}_\mathrm{entry}
     *   = T_{b\rightarrow e}\,\mathbf{r}_\mathrm{body},
     * \qquad
     * \mathbf{r}_\mathrm{exit}
     *   = T_{b\rightarrow x}\,\mathbf{r}_\mathrm{body},
     * \qquad
     * \mathbf{r}_\mathrm{field}
     *   =
     *   \begin{cases}
     *     \mathcal{C}(\mathbf{r}_\mathrm{entry}), &
     *       z_\mathrm{exit}<-\ell_{\mathrm{fringe},x},\\
     *     (x_\mathrm{exit}, y_\mathrm{exit},
     *       L_\mathrm{ref}+z_\mathrm{exit}), &
     *       z_\mathrm{exit}\ge-\ell_{\mathrm{fringe},x},
     *   \end{cases}
     * \f]
     * where \f$\mathcal{C}\f$ is convertEntryCartesianToFieldLocal(). This
     * prevents the RBEND exit fringe from being measured from the upstream
     * circular-arc exit and removes the field discontinuity at the rectangular
     * exit plane.
     */
    Vector_t<double, 3> convertBodyCartesianToFieldLocal(
            const Vector_t<double, 3>& bodyCartesian) const;

    /**
     * @brief Rotate a body-frame Cartesian vector into the local tangent frame.
     *
     * Explicitly placed sector bends expose a rigid body frame for placement
     * and rendering, but the runtime field chart is attached to the entrance
     * tangent and extended past the body through the exit tangent. To keep
     * momentum rotation consistent with convertBodyCartesianToFieldLocal(), the
     * vector is first re-expressed in the rigid entrance or exit frame defined
     * by the actual bend geometry:
     * \f[
     * \mathbf{v}_\mathrm{entry} = R_{b\rightarrow e}\,\mathbf{v}_\mathrm{body},
     * \qquad
     * \mathbf{v}_\mathrm{exit} = R_{b\rightarrow x}\,\mathbf{v}_\mathrm{body},
     * \f]
     * and only the interior interval \f$0 < s < L\f$ is transported through
     * the curvilinear tangent basis.
     */
    Vector_t<double, 3> rotateBodyCartesianVectorToFieldLocal(
            const Vector_t<double, 3>& bodyVector, double s) const;

    /**
     * @brief Rotate a vector from the rigid entry frame into the local tangent frame.
     *
     * The analytic bend field is expressed in a tangent basis that follows the
     * reference path. If \f$s\f$ denotes the reference-path coordinate and
     * \f$h_\mathrm{geom}\f$ the geometric curvature, the local tangent frame is
     * obtained from the rigid entry frame by a rotation about \f$y\f$ through
     * \f$\phi = h_\mathrm{geom} s\f$:
     * \f[
     * \mathbf{v}_\mathrm{field} = R_y(\phi)\,\mathbf{v}_\mathrm{entry}.
     * \f]
     */
    Vector_t<double, 3> rotateEntryCartesianVectorToFieldLocal(
            const Vector_t<double, 3>& entryVector, double s) const;

    /**
     * @brief Rotate a vector from the local tangent frame back into the rigid entry frame.
     *
     * This applies the inverse of rotateEntryCartesianVectorToFieldLocal():
     * \f[
     * \mathbf{v}_\mathrm{entry} = R_y(-\phi)\,\mathbf{v}_\mathrm{field},
     * \qquad \phi = h_\mathrm{geom} s.
     * \f]
     */
    Vector_t<double, 3> rotateFieldLocalVectorToEntryCartesian(
            const Vector_t<double, 3>& fieldVector, double s) const;

    /**
     * @brief Apply the analytic bend field to particles already expressed in one slice frame.
     *
     * The particles are assumed to be in the slice-local tangent frame
     * associated with `slice`. The kernel reconstructs the corresponding
     * entry-frame position, maps it into the bend field chart
     * \f$(x, y, s)\f$, evaluates the analytic bend multipole field there, and
     * rotates the resulting field back into the slice-local frame before
     * accumulation.
     *
     * @note This path is used only for many-particle tracking. It does not
     * alter the reference-particle field evaluation.
     */
    bool applySlice(
            const std::shared_ptr<ParticleContainer_t>& pc, const TrackingSlice& slice) const;

    /**
     * @brief Return the OPAL Enge fringe support in front of the body.
     */
    double getEntryFringeSupportLength() const;

    /**
     * @brief Return the OPAL Enge fringe support behind the body.
     */
    double getExitFringeSupportLength() const;

    /**
     * @brief Return the effective integrated field length of the fringe model.
     *
     * The value is the numerical integral of the OPAL default Enge envelope
     * over the full local support. It is cached whenever the bend length, gap,
     * bend angle, or pole-face angles change. Rectangular bends use this
     * integral for field normalization; sector-bend parsers may use a
     * chord-to-arc normalization instead to match historical OPAL input.
     */
    double getEffectiveFieldLength() const;

    /**
     * @brief Return the local OPAL Enge field scale.
     *
     * The scale is zero outside the field support. At each pole face it follows
     * OPAL's default `FM1DProfile1` Enge function and equals approximately
     * \f$f(0)=0.3825\f$ at the nominal edge.
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

    /**
     * @brief Return the geometric curvature of the design reference path.
     *
     * The sliced many-particle tracker uses rigid tangent frames built on the
     * geometric reference path of the bend. That path is defined by the bend
     * angle \f$\theta\f$ and the reference-path body length \f$L_\mathrm{ref}\f$,
     * so the geometric curvature is
     * \f[
     * h_\mathrm{geom} = \frac{\theta}{L_\mathrm{ref}}.
     * \f]
     *
     * This quantity is distinct from the effective curvature
     * \f$\theta / L_\mathrm{eff}\f$ used to normalize the integrated field
     * strength when fringe support extends the effective field length.
     */
    double getReferencePathCurvature() const;

    double getSignedCurvature() const;
    /**
     * @brief Return the hard-edge horizontal entry coefficient.
     *
     * This diagnostic coefficient is the textbook \f$h\tan(E_1)\f$ value. The
     * OPAL-compatible distributed Enge field does not add it as a separate
     * \f$B_y \propto x\f$ term; the runtime field keeps the same one-dimensional
     * longitudinal field shape as OPAL `Bend2D`.
     */
    double getEntryEdgeHorizontalStrength() const;
    /**
     * @brief Return the hard-edge horizontal exit coefficient.
     *
     * This is the exit-face counterpart of
     * `getEntryEdgeHorizontalStrength()` and is retained for diagnostics and
     * possible future optics models, not applied as an extra distributed field
     * gradient in the OPAL-compatible fringe map.
     */
    double getExitEdgeHorizontalStrength() const;
    /**
     * @brief Return the signed entry-edge vertical kick strength.
     *
     * The implemented hard-edge fringe optics model uses the first-order kick
     * law
     * \f[
     *   y' \mathrel{+}= k_{y,\mathrm{entry}}\,y,
     *   \qquad
     *   k_{y,\mathrm{entry}} = -h \tan(E_1 - \psi_1),
     * \f]
     * where \f$h\f$ is the effective curvature and
     * \f[
     *   \psi_1 = h \, \mathrm{HGAP} \, \mathrm{FINT}
     *          \frac{1 + \sin^2(E_1)}{\cos(E_1)}.
     * \f]
     *
     * This quantity is a kick coefficient, not a magnetic-field gradient. The
     * distributed fringe-field model converts it into \f$B_x\f$ using the local
     * reference rigidity.
     */
    double getEntryEdgeVerticalStrength() const;
    /**
     * @brief Return the signed exit-edge vertical kick strength.
     *
     * This is the exit-face counterpart of
     * `getEntryEdgeVerticalStrength()`,
     * \f[
     *   k_{y,\mathrm{exit}} = -h \tan(E_2 - \psi_2),
     * \f]
     * with the same \f$\mathrm{HGAP}\f$/\f$\mathrm{FINT}\f$ correction model.
     * As for the entry face, the returned value is a kick coefficient that must
     * be converted to a field gradient through the local rigidity.
     */
    double getExitEdgeVerticalStrength() const;
    /**
     * @brief Return the entry-fringe coefficient for vertical edge focusing.
     *
     * The vertical edge-strength methods return the optics kick coefficient
     * \f$k_y\f$ in
     * \f[
     *   y'_\mathrm{out} = y'_\mathrm{in} + k_y y .
     * \f]
     * OPAL's analytic fringe model spreads this equivalent force over the same
     * Enge ramp used for the main dipole field.  With local field scale
     * \f$\lambda(s)\f$ and support \f$[s_-,s_+]\f$, the entry field is
     * \f[
     *   B_x^\mathrm{edge}(s,y) =
     *   \frac{B_0}{h}\,
     *   \frac{k_y}{\lambda(s_+)-\lambda(s_-)}
     *   \frac{d\lambda}{ds}\,y .
     * \f]
     * This preserves the integrated hard-edge kick while placing the force
     * over the physical fringe instead of upstream of the nominal edge only.
     */
    double getEntryEdgeVerticalFieldCoefficient() const;
    /**
     * @brief Return the exit-fringe coefficient for vertical edge focusing.
     *
     * At the exit \f$d\lambda/ds < 0\f$, so callers apply the returned
     * coefficient with the opposite derivative sign in order to preserve the
     * signed exit kick.
     */
    double getExitEdgeVerticalFieldCoefficient() const;
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
    double effectiveFieldLength_m;
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

inline void BendBase::setFullGap(double gap) {
    gap_m = std::abs(gap);
    updateFieldSupportExtent();
}

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

inline void BendBase::setFieldMapFN(std::string fileName) {
    fileName_m = std::move(fileName);
    updateFieldSupportExtent();
}

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
    return std::max(0.0, fieldEnd_m - getReferencePathLength());
}

inline double BendBase::getReferencePathLength() const { return getElementLength(); }

inline double BendBase::getEffectiveFieldLength() const { return effectiveFieldLength_m; }

inline double BendBase::getStoredExitAngle() const { return exitAngle_m; }

#endif  // OPALX_BendBase_HH
