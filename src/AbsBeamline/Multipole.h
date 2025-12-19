#ifndef CLASSIC_Multipole_HH
#define CLASSIC_Multipole_HH

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BMultipoleField.h"

class Fieldmap;
constexpr int MAX_MP_ORDER = 2;

/**
 * @class Multipole
 * @brief Interface for general multipole.
 *
 * Class Multipole defines the abstract interface for magnetic multipoles.
 * The order n of multipole components runs from 1 to N and is dynamically
 * adjusted. It is connected with the number of poles by the following table:
 *
 * | Order (n) | Name                |
 * |-----------|---------------------|
 * | 1         | dipole              |
 * | 2         | quadrupole          |
 * | 3         | sextupole           |
 * | 4         | octupole            |
 * | 5         | decapole            |
 *
 * Units for multipole strengths are Teslas / m^(n-1).
 */

class Multipole : public Component {
public:
/* ============================== Constructors ============================== */
    explicit Multipole(const std::string& name);
    Multipole();
    Multipole(const Multipole&);
    virtual ~Multipole();
/* ========================================================================== */
/* =========================== Multipole Expansion ========================== */

    // @note n=0 corresponds to the dipol, n=1 to the quadrupole, etc...

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
     * 
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply() override;

    /**
     * @brief Apply to particle i
     * 
     * @param i Particle index
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     * 
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply(
        const size_t& i, 
        const double& t, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;
   
    /**
     * @brief Apply to particle with position R and momentum P
     * 
     * @param R Position 
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     * 
     * @returns true if particle is out-of-bounds (lost), false otherwise
     */
    virtual bool apply(
        const Vector_t<double, 3>& R, 
        const Vector_t<double, 3>& P, 
        const double& t,
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;
    
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
        const Vector_t<double, 3>& R, 
        const Vector_t<double, 3>& P, 
        const double& t,
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;
/* ========================================================================== */
/* ============================== Functions ================================= */
    // @brief Apply visitor to Multipole.
    virtual void accept(BeamlineVisitor&) const override;

    // @brief Get multipole field.
    virtual BMultipoleField& getField() override = 0;

    // @breif Get multipole field. Version for const object.
    virtual const BMultipoleField& getField() const override = 0;
    
    // @returns Is the n-th component focusing?
    bool isFocusing(int n) const;

    /**
     * @brief Setup, multipole goes online
     * 
     * @param bunch The reference bunch
     * @param startField Where the fields start along the path
     * @param endFied Where the fields end along the path
     */
    virtual void initialise(
        PartBunch_t* bunch, 
        double& startField, 
        double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool isInside(const Vector_t<double, 3>& r) const override;
/* ========================================================================== */
/* =========================== Unused Functions ============================= */
    // @returns StraightGeometry
    virtual StraightGeometry& getGeometry() override = 0;

    // @returns StraightGeometry
    virtual const StraightGeometry& getGeometry() const override = 0;

    // @brief Set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);

    // @breif Get number of slices for map tracking 
    std::size_t getNSlices() const;
    
    // Not implemented.
    void operator=(const Multipole&);


/* ========================================================================== */
private:
    /**
     * @brief Computes the E and B field at position R
     * 
     * @note Is called inside the kernel launched by Multipole::apply()
     * 
     * @param R Position
     * @param E Reference to electric field
     * @param B Reference to magnetic field
     * @param NormalComponents View of normal components
     * @param SkewComponents View of skew components
     * @param max_NormalComponent Max order of normal components
     * @param max_SkewComponent Max order of skew components
     */
    KOKKOS_INLINE_FUNCTION
    static void computeField(
        Vector_t<double, 3> R, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B,
        const Kokkos::View<double*>& NormalComponents,
        const Kokkos::View<double*>& SkewComponents,
        int max_NormalComponent,
        int max_SkewComponent);

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
        Vector_t<double, 3> R, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) const;
/* ========================================================================== */
/* =========================== Variables ==================================== */ 
    Kokkos::View<double*> NormalComponents;
    Kokkos::View<double*> NormalComponentErrors;
    Kokkos::View<double*> SkewComponents;
    Kokkos::View<double*> SkewComponentErrors;
    unsigned int max_SkewComponent_m;
    unsigned int max_NormalComponent_m;
/* =========================== Unused Variables ============================= */ 
    std::size_t nSlices_m;
/* ========================================================================== */
};
/* =========================== Inline Functions ============================= */ 
// @note Inlined function bodies need to be in the header for libraries

// @returns n-th normal component of the multipole expansion
inline void Multipole::setNormalComponent(int n, double v) {
    setNormalComponent(n, v, 0.0);
}
// @returns n-th skew component of the multipole expansion
inline void Multipole::setSkewComponent(int n, double v) {
    setSkewComponent(n, v, 0.0);
}

inline size_t Multipole::getMaxNormalComponentIndex() const {
    return NormalComponents.size();
}

inline size_t Multipole::getMaxSkewComponentIndex() const {
    return SkewComponents.size();
}
/* ========================================================================== */
#endif  // CLASSIC_Multipole_HH
