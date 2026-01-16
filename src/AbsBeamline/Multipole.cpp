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

#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/GeneralClassicException.h"

// Unused Headers
#include "Physics/Physics.h"
#include "Fields/Fieldmap.h"

// orders
namespace {
    enum { DIPOLE, QUADRUPOLE, SEXTUPOLE, OCTUPOLE, DECAPOLE };
}

// GPU factorial
namespace {
    KOKKOS_INLINE_FUNCTION
    double factorial(unsigned int n) {
        if (n == 0) return 1.0;
        if (n == 1) return 1.0;
        if (n == 2) return 2.0;
        if (n == 3) return 6.0;
        if (n == 4) return 24.0;
        if (n == 5) return 120.0;
        if (n > 5) Kokkos::abort("factorial out of bounds");
    }
}

/* ============================== Constructors ============================== */
Multipole::Multipole() : Multipole("") {
}

/**
 * @brief Copy constructor.
 * 
 * @param right Multipole to copy
 */
Multipole::Multipole(const Multipole& right)
    : Component(right),
      NormalComponents("Multipole::Normal", 
        right.NormalComponents.extent(0)),
      NormalComponentErrors("Multipole::NormalErr", 
        right.NormalComponentErrors.extent(0)),
      SkewComponents("Multipole::Skew", 
        right.SkewComponents.extent(0)),
      SkewComponentErrors("Multipole::SkewErr", 
        right.SkewComponentErrors.extent(0)),
      
      max_SkewComponent_m(right.max_SkewComponent_m),
      max_NormalComponent_m(right.max_NormalComponent_m),
     
      nSlices_m(right.nSlices_m) 
{
    Kokkos::deep_copy(NormalComponents, right.NormalComponents);
    Kokkos::deep_copy(NormalComponentErrors, right.NormalComponentErrors);
    Kokkos::deep_copy(SkewComponents, right.SkewComponents);
    Kokkos::deep_copy(SkewComponentErrors, right.SkewComponentErrors);
}

/**
 * @brief Constructor with given name.
 * 
 * @param name Name of the Multipole element
 */
Multipole::Multipole(const std::string& name)
    : Component(name),
      NormalComponents("Multipole::Normal", 1),
      NormalComponentErrors("Multipole::NormalErr", 1),
      SkewComponents("Multipole::Skew", 1),
      SkewComponentErrors("Multipole::SkewErr", 1),
      max_SkewComponent_m(1),
      max_NormalComponent_m(1),
      nSlices_m(1) {
    
    // Initialize with 0.0
    Kokkos::deep_copy(NormalComponents, 0.0);
    Kokkos::deep_copy(SkewComponents, 0.0);
}

Multipole::~Multipole() {
}
/* ========================================================================== */
/* =========================== Multipole Expansion ========================== */

// @note n=0 corresponds to the dipol, n=1 to the quadrupole, etc...

/**
 * @brief Gets the n-th normal component of the multipole expansion 
 * 
 * Copies the device view NormalComponents n-th entry to a host
 * variable and returns it. For host-only execution the deepcopy is ignored 
 * and the n-th entry is returned directly. 
 * 
 * @param n Component
 * @returns n-th normal component of the multipole expansion
 */
double Multipole::getNormalComponent(int n) const 
{
    if (n < 0) {
        throw GeneralClassicException("Multipole::getNormalComponent", 
            "component index " + std::to_string(n) + " out of bounds");
    }
    else if (n < max_NormalComponent_m) {
        double val;
        Kokkos::deep_copy(val, Kokkos::subview(NormalComponents, n));
        return val;
    }
    else {
        return 0.0;
    }
}

/**
 * @brief Gets the n-th skew component of the multipole expansion 
 * 
 * Copies the device view SkewComponents n-th entry to a host
 * variable and returns it. For host-only execution the deepcopy is ignored 
 * and the n-th entry is returned directly.
 * 
 * @param n Component
 * @returns n-th skew component of the multipole expansion
 */
double Multipole::getSkewComponent( int n) const 
{
    if (n < 0) {
        throw GeneralClassicException("Multipole::getSkewComponent", 
            "component index " + std::to_string(n) + " out of bounds");
    }
    else if (n < max_SkewComponent_m) {
        double val;
        Kokkos::deep_copy(val, Kokkos::subview(SkewComponents, n));
        return val;
    }
    else{
        return 0.0;
    }
}

/**
 * @brief Sets the normal component of the multipole field.
 *
 * If the order is higher than the previous highest order, the view is resized
 * and the new entries (components) are 0. The value is then copied into the 
 * NormalComponents view.
 *
 * @param n The multipole order 
 * @param v The nominal value of the normal component.
 * @param vError The error associated with the normal component.
 */
void Multipole::setNormalComponent(int n, double v, double vError) 
{
    if (n < 0) {
        throw GeneralClassicException("Multipole::setNormalComponent", 
            "component index " + std::to_string(n) + " out of bounds");
    }

    if (n >= max_NormalComponent_m) {
        max_NormalComponent_m = n+1;
        Kokkos::resize(NormalComponents, max_NormalComponent_m);
        Kokkos::resize(NormalComponentErrors, max_NormalComponent_m);
    }

    double valComp, valErr;
    
    // TODO: Check the physics here (copied from OPAL)
    if (n == DIPOLE) {
        valComp = (v + vError) / 2.0;
        valErr = valComp;
    } else {
        double fact = factorial(n);
        valComp = (v + vError) / fact;
        valErr = vError / fact;
    }

    // Copy to Device
    auto sub_comp = Kokkos::subview(NormalComponents, n);
    auto sub_err  = Kokkos::subview(NormalComponentErrors, n);
    
    Kokkos::deep_copy(sub_comp, valComp);
    Kokkos::deep_copy(sub_err, valErr);
}

/**
 * @brief Sets the skew component of the multipole field.
 *
 * If the order is higher than the previous highest order, the view is resized
 * and the new entries (components) are 0. The value is then copied into the 
 * SkewComponents view.
 *
 * @param n The multipole order 
 * @param v The nominal value of the skew component.
 * @param vError The error associated with the skew component.
 */
void Multipole::setSkewComponent(int n, double v, double vError) 
{
    if (n < 0) {
        throw GeneralClassicException("Multipole::setSkewComponent", 
            "component index " + std::to_string(n) + " out of bounds");
    }

    if (n >= max_SkewComponent_m) {
        max_SkewComponent_m = n+1;
        Kokkos::resize(SkewComponents, max_SkewComponent_m);
        Kokkos::resize(SkewComponentErrors, max_SkewComponent_m);
    }

    double valComp, valErr;
    // TODO: Check physics (copied from OPAL)
    if (n == DIPOLE) {
        valComp = (v + vError) / 2.0;
        valErr = valComp;
    } else {
        double fact = factorial(n);
        valComp = (v + vError) / fact;
        valErr = vError / fact;
    }

    auto sub_comp = Kokkos::subview(SkewComponents, n);
    auto sub_err  = Kokkos::subview(SkewComponentErrors, n);

    Kokkos::deep_copy(sub_comp, valComp);
    Kokkos::deep_copy(sub_err, valErr);
}


/* ========================================================================== */
/* ============================== Apply Functions =========================== */

/**
 * @brief Applies the multipole field to all particles inside the magnet bounds
 * 
 * @note The kernel launch is moved inside this functions for GPU compatibility
 * 
 * @note TODO: Out-of-bounds check 
 *  
 * @returns true if particle is out-of-bounds (lost), false otherwise
 */
bool Multipole::apply()
{
    std::cout << "Multipole::apply() called" <<std::endl;
    // Get the particle container
    std::shared_ptr<ParticleContainer_t> pc = 
        RefPartBunch_m->getParticleContainer();

    // Get the views
    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    // Local variables that are copied into the kernel
    double elemLength = getElementLength();

    // Capture member variables by value for the kernel
    auto normalComponents = NormalComponents;
    auto skewComponents = SkewComponents;
    int maxNormal = max_NormalComponent_m;
    int maxSkew = max_SkewComponent_m;

    // Kernel launch over all particles
    Kokkos::parallel_for("Multipole::apply()", ippl::getRangePolicy(Rview), 
    KOKKOS_LAMBDA(const int i)
    {
        // Check bounds
        if (Rview(i)(2) > 0 && Rview(i)(2) <= elemLength){
            Vector_t<double,3> Ef(0.0), Bf(0.0);
            // Compute field at particle position
            computeField(Rview(i), Ef, Bf, 
                normalComponents, skewComponents,maxNormal, maxSkew);
            for(unsigned d=0; d<3; ++d){
                Eview(i)(d) += Ef(d);
                Bview(i)(d) += Bf(d);
            }
        }    
    });
    return false;
}

/**
 * @brief Applies the multipole field at particle i's position to E and B
 * 
 * @note Cannot be inside a kernel -> GPU incompatible
 * 
 * @param i Particle index
 * @param E Electric field reference
 * @param B Magnetic field reference
 * @returns true if particle is out-of-bounds (lost), false otherwise 
 */
bool Multipole::apply(
    const size_t& i, 
    const double&, 
    Vector_t<double, 3>& E, 
    Vector_t<double, 3>& B) 
{
    // Get container
    std::shared_ptr<ParticleContainer_t> pc = 
    RefPartBunch_m->getParticleContainer();
    auto Rview = pc->R.getView();
    const Vector_t<double, 3> R = Rview(i);

    // Check bounds
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    // Compute fields
    Vector_t<double, 3> Ef(0.0), Bf(0.0);
    computeFieldHost(R, Ef, Bf);

    // Apply fields
    for (unsigned int d = 0; d < 3; ++d) {
        E[d] += Ef(d);
        B[d] += Bf(d);
    }

    return false;
}

/**
 * @brief Applies the multipole field at position R to E and B
 * 
 * @note Cannot be inside a kernel -> GPU incompatible
 * 
 * @param R Position
 * @param E Electric field reference
 * @param B Magnetic field reference
 * @returns true if particle is out-of-bounds (lost), false otherwise 
 */

bool Multipole::apply(
    const Vector_t<double, 3>& R, 
    const Vector_t<double, 3>&, 
    const double&, 
    Vector_t<double, 3>& E,
    Vector_t<double, 3>& B) 
{
    // Check bounds
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    // Compute field
    computeFieldHost(R, E, B);

    return false;
}


/**
 * @brief Applies the multipole field at reference particle 
 * position R to E and B
 * 
 * @note Host-only function for the orbitthreader 
 * 
 * @param R Position
 * @param E Electric field reference
 * @param B Magnetic field reference
 * @returns true if particle is out-of-bounds (lost), false otherwise 
 */
bool Multipole::applyToReferenceParticle(
    const Vector_t<double, 3>& R, 
    const Vector_t<double, 3>&, 
    const double&, 
    Vector_t<double, 3>& E,
    Vector_t<double, 3>& B) 
{
    // Check bounds
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return true;

    /*
    for (int i = 0; i < max_NormalComponent_m; ++i) {
        NormalComponents[i] -= NormalComponentErrors[i];
    }
    for (int i = 0; i < max_SkewComponent_m; ++i) {
        SkewComponents[i] -= SkewComponentErrors[i];
    }
    */

    // Compute field
    computeFieldHost(R, E, B);

    /*
    for (int i = 0; i < max_NormalComponent_m; ++i) {
        NormalComponents[i] += NormalComponentErrors[i];
    }
    for (int i = 0; i < max_SkewComponent_m; ++i) {
        SkewComponents[i] += SkewComponentErrors[i];
    }
    */
    

    return false;
}
/* ========================================================================== */
/* ============================== Functions ================================= */
/**
 * @brief Accepts a BeamlineVisitor
 * 
 * @param visitor BeamlineVisitor reference
 */
void Multipole::accept(BeamlineVisitor& visitor) const {
    visitor.visitMultipole(*this);
}

/**
 * @brief Computes the multipole field at position R
 * 
 * @note NormalComponents and SkewComponents need to be passed in the arguments
 * because member views are not accessible inside Kokkos kernels.
 * 
 * @note Only implements up to quadrupole for now
 * 
 * @param R Position
 * @param E Electric field reference
 * @param B Magnetic field reference
 * @param NormalComponents Kokkos::View of normal components
 * @param SkewComponents Kokkos::View of skew components
 * @param max_NormalComponent Maximum normal component index
 * @param max_SkewComponent Maximum skew component index
 */
KOKKOS_INLINE_FUNCTION
void Multipole::computeField(
    Vector_t<double, 3> R, 
    Vector_t<double, 3>& E, 
    Vector_t<double, 3>& B,
    const Kokkos::View<double*>& NormalComponents,
    const Kokkos::View<double*>& SkewComponents,
    int max_NormalComponent,
    int max_SkewComponent) 
{
    // Use primitive double arrays instead of Vector_t to avoid constructor/destructor calls on GPU stack
    double Rn_x[MAX_MP_ORDER + 1];
    double Rn_y[MAX_MP_ORDER + 1];
    double fact[MAX_MP_ORDER + 1];

    // --- NORMAL COMPONENTS ---
    {
        Rn_x[0] = 1.0;
        Rn_y[0] = 1.0; // Equivalent to Vector_t(1.0)
        fact[0] = 1.0;
        
        // Use local variable for loop bound to help optimizer
        int count = (max_NormalComponent > MAX_MP_ORDER) ? 
            MAX_MP_ORDER : max_NormalComponent;

        for (int i = 0; i < count; ++i) {
            // Access Kokkos::View directly (device memory)
            double norm = NormalComponents(i);

            switch (i) {
                case DIPOLE:
                    B(1) += norm;
                    break;

                case QUADRUPOLE:
                    B(0) += norm * R(1);
                    B(1) += norm * R(0);
                    break;
            }

            Rn_x[i + 1] = Rn_x[i] * R(0);
            Rn_y[i + 1] = Rn_y[i] * R(1);
            fact[i + 1]  = fact[i] / (double)(i + 1);
        }
    }

    // --- SKEW COMPONENTS ---
    {
        Rn_x[0] = 1.0;
        Rn_y[0] = 1.0;
        fact[0] = 1.0;
        
        int count = (max_SkewComponent > MAX_MP_ORDER) ? 
            MAX_MP_ORDER : max_SkewComponent;

        for (int i = 0; i < count; ++i) {
            double skew = SkewComponents(i);

            switch (i) {
                case DIPOLE:
                    B(0) -= skew;
                    break;

                case QUADRUPOLE:
                    B(0) -= skew * R(0);
                    B(1) += skew * R(1);
                    break;
            }

            Rn_x[i + 1] = Rn_x[i] * R(0);
            Rn_y[i + 1] = Rn_y[i] * R(1);
            fact[i + 1]  = fact[i] / (double)(i + 1);
        }
    }
}

/**
 * @brief Computes the multipole field at position R (host version)
 * 
 * @note Host-only function for the orbitthreader 
 * 
 * @param R Position
 * @param E Electric field reference
 * @param B Magnetic field reference
 */
void Multipole::computeFieldHost(
    const Vector_t<double, 3> R, 
    Vector_t<double, 3>& /*E*/, 
    Vector_t<double, 3>& B) const 
{
    Vector_t<double, 3> Rn[MAX_MP_ORDER + 1];
    double fact[MAX_MP_ORDER + 1];

    auto NormalComponents_host = Kokkos::create_mirror_view(NormalComponents);
    Kokkos::deep_copy(NormalComponents_host, NormalComponents);

    auto SkewComponents_host = Kokkos::create_mirror_view(SkewComponents);
    Kokkos::deep_copy(SkewComponents_host, SkewComponents);

    // --- NORMAL COMPONENTS ---
    {
        Rn[0] = Vector_t<double, 3>(1.0);
        fact[0] = 1.0;
        
        // Use local variable for loop bound to help optimizer
        int count = (max_NormalComponent_m > MAX_MP_ORDER) ? 
            MAX_MP_ORDER : max_NormalComponent_m;

        for (int i = 0; i < count; ++i) {
            // Access Kokkos::View directly (device memory)
            double norm = NormalComponents_host(i);

            switch (i) {
                case DIPOLE:
                    B(1) += norm;
                    break;

                case QUADRUPOLE:
                    B(0) += norm * R(1);
                    B(1) += norm * R(0);
                    break;

                case SEXTUPOLE:
                    B(0) += 2 * norm * R(0) * R(1);
                    B(1) += norm * (Rn[2](0) - Rn[2](1));
                    break;

                case OCTUPOLE:
                    B(0) += norm * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                    B(1) += norm * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                    break;

                case DECAPOLE:
                    B(0) += 4 * norm * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                    B(1) += norm * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                    break;

                default: {
                    double powMinusOne = 1.0;
                    double Bx = 0.0, By = 0.0;
                    for (int j = 1; j <= (i + 1) / 2; ++j) {
                        Bx += powMinusOne * norm
                              * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] * Rn[2 * j - 1](1)
                                 * fact[2 * j - 1]);
                        By += powMinusOne * norm
                              * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] * Rn[2 * j - 2](1)
                                 * fact[2 * j - 2]);
                        powMinusOne *= -1.0;
                    }

                    if ((i + 1) / 2 == i / 2) {
                        int j = (i + 2) / 2;
                        By += powMinusOne * norm
                              * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] * Rn[2 * j - 2](1)
                                 * fact[2 * j - 2]);
                    }
                    B(0) += Bx;
                    B(1) += By;
                }
            }

            Rn[i + 1](0) = Rn[i](0) * R(0);
            Rn[i + 1](1) = Rn[i](1) * R(1);
            fact[i + 1]  = fact[i] / (double)(i + 1);
        }
    }

    // --- SKEW COMPONENTS ---
    {
        Rn[0] = Vector_t<double, 3>(1.0);
        fact[0] = 1.0;
        
        int count = (max_SkewComponent_m > MAX_MP_ORDER) ? 
            MAX_MP_ORDER : max_SkewComponent_m;

        for (int i = 0; i < count; ++i) {
            double skew = SkewComponents_host(i);

            switch (i) {
                case DIPOLE:
                    B(0) -= skew;
                    break;

                case QUADRUPOLE:
                    B(0) -= skew * R(0);
                    B(1) += skew * R(1);
                    break;

                case SEXTUPOLE:
                    B(0) -= skew * (Rn[2](0) - Rn[2](1));
                    B(1) += 2 * skew * R(0) * R(1);
                    break;

                case OCTUPOLE:
                    B(0) -= skew * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                    B(1) += skew * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                    break;

                case DECAPOLE:
                    B(0) -= skew * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                    B(1) += 4 * skew * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                    break;

                default: {
                    double powMinusOne = 1.0;
                    double Bx = 0.0, By = 0.0;
                    for (int j = 1; j <= (i + 1) / 2; ++j) {
                        Bx -= powMinusOne * skew
                              * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] * Rn[2 * j - 2](1)
                                 * fact[2 * j - 2]);
                        By += powMinusOne * skew
                              * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] * Rn[2 * j - 1](1)
                                 * fact[2 * j - 1]);
                        powMinusOne *= -1.0;
                    }

                    if ((i + 1) / 2 == i / 2) {
                        int j = (i + 2) / 2;
                        Bx -= powMinusOne * skew
                              * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] * Rn[2 * j - 2](1)
                                 * fact[2 * j - 2]);
                    }

                    B(0) += Bx;
                    B(1) += By;
                }
            }

            Rn[i + 1](0) = Rn[i](0) * R(0);
            Rn[i + 1](1) = Rn[i](1) * R(1);
            fact[i + 1]  = fact[i] / (double)(i + 1);
        }
    }
}
/**
     * @brief Setup, multipole goes online
     * 
     * @param bunch The reference bunch
     * @param startField Where the fields start along the path
     * @param endFied Where the fields end along the path
     */
void Multipole::initialise(
    PartBunch_t* bunch, 
    double& startField, 
    double& endField) 
{
    RefPartBunch_m = bunch;
    endField       = startField + getElementLength();
    online_m       = true;
}

void Multipole::finalise() {
    online_m = false;
}

bool Multipole::bends() const {
    return false;
}

void Multipole::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd   = getElementLength();
}

ElementType Multipole::getType() const {
    return ElementType::MULTIPOLE;
}

bool Multipole::isInside(const Vector_t<double, 3>& r) const {
    if (r(2) >= 0.0 && r(2) < getElementLength()) {
        return isInsideTransverse(r);
    }

    return false;
}

bool Multipole::isFocusing(int component) const {
    if (component < 0)
        throw GeneralClassicException("Multipole::isFocusing", "component negative");
    else if (component >= NormalComponents.extent(0))
        throw GeneralClassicException("Multipole::isFocusing", "component too big");

    // Fix: Use getNormalComponent() to safely retrieve the value from GPU memory
    return getNormalComponent(component) * std::pow(-1, component + 1)
               * RefPartBunch_m->getChargePerParticle() > 0.0;
}

/* ========================================================================== */
/* =========================== Unused Functions ============================= */

// set the number of slices for map tracking
void Multipole::setNSlices(const std::size_t& nSlices) {
    nSlices_m = nSlices;
}

// get the number of slices for map tracking
std::size_t Multipole::getNSlices() const {
    return nSlices_m;
}
/* ========================================================================== */


