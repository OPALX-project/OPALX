#ifndef OPALX_Multipole_HH
#define OPALX_Multipole_HH

#include "AbsBeamline/ElementBase.h"

class Fieldmap;
constexpr int MAX_MP_ORDER = 5;

/**
 * @class Multipole
 * @brief Interface for general multipole.
 *
 * Class Multipole defines the abstract interface for magnetic multipoles.
 * With normal(n) and skew(n) the stored components, the field is
 * By + i Bx = sum over n of (normal(n) - i skew(n)) (x + i y)^n:
 *
 * | n | Name       |
 * |---|------------|
 * | 0 | dipole     |
 * | 1 | quadrupole |
 * | 2 | sextupole  |
 * | 3 | octupole   |
 * | 4 | decapole   |
 *
 * Orders up to n = MAX_MP_ORDER - 1 (decapole) are used. The setters store the given
 * strength divided by n!, so normal(n) is in T/m^n. The MULTIPOLE input (OpalMultipole)
 * does not allow a dipole, because the element is placed as a straight element.
 */

class Multipole : public ElementBase {
public:
    /* ============================== Constructors ============================== */
    explicit Multipole(const std::string& name);
    Multipole();
    Multipole(const Multipole&);
    virtual ~Multipole();
    /* ========================================================================== */
    /* =========================== Multipole Expansion ========================== */

    // @note n = 0 is the dipole, n = 1 the quadrupole, up to n = 4 the decapole.

    // @returns n-th normal component of the multipole expansion
    double getNormalComponent(int n) const;

    // @brief Sets the n-th normal component
    void setNormalComponent(int n, double);

    // @brief Sets the n-th normal component with an error
    void setNormalComponent(int n, double, double);

    // @returns n-th skew component of the multipole expansion
    double getSkewComponent(int n) const;

    // @brief Sets the n-th skew component
    void setSkewComponent(int n, double);

    // @brief Sets the n-th skew component with error
    void setSkewComponent(int n, double, double);

    // @brief Get index of maximal components
    size_t getMaxNormalComponentIndex() const;
    size_t getMaxSkewComponentIndex() const;

    /* ========================================================================== */
    /* ============================== Apply Functions =========================== */
    /**
     * @brief Apply to all particles. Kernel launch moved inside the function.
     */
    virtual void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /**
     * @brief Apply to particle with position R and momentum P
     *
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     */
    virtual void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /**
     * @brief Apply to reference particle with position R and momemtum P
     *
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     *
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    /* ========================================================================== */
    /* ============================== Functions ================================= */
    // @brief Apply visitor to Multipole.
    virtual void accept(BeamlineVisitor&) const override;

    // @returns Is the n-th component focusing?
    bool isFocusing(int n) const;

    /**
     * @brief Setup, multipole goes online
     *
     * @param bunch The reference bunch
     * @param startField Where the fields start along the path
     * @param endFied Where the fields end along the path
     */
    virtual void initialise(PartBunch_t* bunch) override;

    virtual void finalise() override;

    virtual ElementType getType() const override;

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    /* ========================================================================== */
    // @brief Set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);

    // @breif Get number of slices for map tracking
    std::size_t getNSlices() const;

    // Not implemented.
    void operator=(const Multipole&);

    /* ========================================================================== */
private:
    /**
     * @brief Adds the multipole field at position R to B
     *
     * With z = x + iy and the stored coefficients kn = normal(n), ks = skew(n), every order
     * n < MAX_MP_ORDER adds Bx += kn Im(z^n) - ks Re(z^n) and By += kn Re(z^n) + ks Im(z^n).
     *
     * @note Used by both the particle kernel in apply() (device views) and computeFieldHost()
     *       (host copies), so particles and the reference particle see the same field.
     *
     * @param R Position
     * @param B Magnetic field, added to
     * @param normal View of normal components
     * @param skew View of skew components
     * @param maxNormal Number of normal components
     * @param maxSkew Number of skew components
     */
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void computeField(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& B, const ViewType& normal,
            const ViewType& skew, int maxNormal, int maxSkew);

    /**
     * @brief Computes the E and B field at position R of the reference particle
     *
     * @note Host version for the orbitthreader
     *
     * @param R Position
     * @param E Reference to electric field
     * @param B Reference to magnetic field
     */
    void computeFieldHost(
            Vector_t<double, 3> R, Vector_t<double, 3>& E, Vector_t<double, 3>& B) const;
    /* ========================================================================== */
    /* =========================== Variables ==================================== */
    Kokkos::View<double*> NormalComponents;
    Kokkos::View<double*> NormalComponentErrors;
    Kokkos::View<double*> SkewComponents;
    Kokkos::View<double*> SkewComponentErrors;
    int max_SkewComponent_m;
    int max_NormalComponent_m;
    /* =========================== Unused Variables ============================= */
    std::size_t nSlices_m;
    /* ========================================================================== */
};
/* =========================== Inline Functions ============================= */
// @note Inlined function bodies need to be in the header for libraries

// @returns n-th normal component of the multipole expansion
inline void Multipole::setNormalComponent(int n, double v) { setNormalComponent(n, v, 0.0); }
// @returns n-th skew component of the multipole expansion
inline void Multipole::setSkewComponent(int n, double v) { setSkewComponent(n, v, 0.0); }

// @returns max index of normal components
inline size_t Multipole::getMaxNormalComponentIndex() const { return NormalComponents.size(); }

// @returns max index of skew components
inline size_t Multipole::getMaxSkewComponentIndex() const { return SkewComponents.size(); }
/* ========================================================================== */
#endif  // OPALX_Multipole_HH
