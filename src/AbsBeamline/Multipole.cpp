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

#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/GeneralOpalException.h"

// Unused Headers
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

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
        return 0;
    }
}  // namespace

/* ============================== Constructors ============================== */
Multipole::Multipole() : Multipole("") {}

/**
 * @brief Copy constructor.
 *
 * @param right Multipole to copy
 */
Multipole::Multipole(const Multipole& right)
    : ElementBase(right),
      NormalComponents("Multipole::Normal", right.NormalComponents.extent(0)),
      NormalComponentErrors("Multipole::NormalErr", right.NormalComponentErrors.extent(0)),
      SkewComponents("Multipole::Skew", right.SkewComponents.extent(0)),
      SkewComponentErrors("Multipole::SkewErr", right.SkewComponentErrors.extent(0)),

      max_SkewComponent_m(right.max_SkewComponent_m),
      max_NormalComponent_m(right.max_NormalComponent_m),

      nSlices_m(right.nSlices_m) {
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
    : ElementBase(name),
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

Multipole::~Multipole() {}
/* ========================================================================== */
/* =========================== Multipole Expansion ========================== */

// @note n = 0 is the dipole, n = 1 the quadrupole, up to n = 4 the decapole.

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
double Multipole::getNormalComponent(int n) const {
    if (n < 0) {
        throw GeneralOpalException(
                "Multipole::getNormalComponent",
                "component index " + std::to_string(n) + " out of bounds");
    } else if (n < max_NormalComponent_m) {
        double val;
        Kokkos::deep_copy(val, Kokkos::subview(NormalComponents, n));
        return val;
    } else {
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
double Multipole::getSkewComponent(int n) const {
    if (n < 0) {
        throw GeneralOpalException(
                "Multipole::getSkewComponent",
                "component index " + std::to_string(n) + " out of bounds");
    } else if (n < max_SkewComponent_m) {
        double val;
        Kokkos::deep_copy(val, Kokkos::subview(SkewComponents, n));
        return val;
    } else {
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
void Multipole::setNormalComponent(int n, double v, double vError) {
    if (n < 0) {
        throw GeneralOpalException(
                "Multipole::setNormalComponent",
                "component index " + std::to_string(n) + " out of bounds");
    }

    if (n >= max_NormalComponent_m) {
        max_NormalComponent_m = n + 1;
        Kokkos::resize(NormalComponents, max_NormalComponent_m);
        Kokkos::resize(NormalComponentErrors, max_NormalComponent_m);
    }

    // Store v / n!, the coefficient of (x + iy)^n in the field; for the dipole (n = 0) this is v.
    const double fact    = factorial(n);
    const double valComp = (v + vError) / fact;
    const double valErr  = vError / fact;

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
void Multipole::setSkewComponent(int n, double v, double vError) {
    if (n < 0) {
        throw GeneralOpalException(
                "Multipole::setSkewComponent",
                "component index " + std::to_string(n) + " out of bounds");
    }

    if (n >= max_SkewComponent_m) {
        max_SkewComponent_m = n + 1;
        Kokkos::resize(SkewComponents, max_SkewComponent_m);
        Kokkos::resize(SkewComponentErrors, max_SkewComponent_m);
    }

    // Store v / n!, the coefficient of (x + iy)^n in the field; for the dipole (n = 0) this is v.
    const double fact    = factorial(n);
    const double valComp = (v + vError) / fact;
    const double valErr  = vError / fact;

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
 */
void Multipole::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    *gmsg << level3 << "Multipole::apply() called" << endl;
    // Get the views
    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    // Local variables that are copied into the kernel
    double elemLength               = getGeometry().getElementLength();
    const ApertureType apertureType = aperture_m.first;
    const double apertureX = aperture_m.second[0], apertureY = aperture_m.second[1];

    // Capture member variables by value for the kernel
    auto normalComponents = NormalComponents;
    auto skewComponents   = SkewComponents;
    int maxNormal         = max_NormalComponent_m;
    int maxSkew           = max_SkewComponent_m;

    // Kernel launch over all particles
    Kokkos::parallel_for(
            "Multipole::apply()", nLocal, KOKKOS_LAMBDA(const size_t i) {
                // Match isInside(): entrance-inclusive, exit-exclusive support.
                if (Rview(i)(2) >= 0 && Rview(i)(2) < elemLength
                    && ApertureHelper::isInsideAperture(
                            Rview(i)(0), Rview(i)(1), apertureType, apertureX, apertureY)) {
                    Vector_t<double, 3> Bf(0.0);
                    computeField(
                            Rview(i), Bf, normalComponents, skewComponents, maxNormal, maxSkew);
                    for (unsigned d = 0; d < 3; ++d) {
                        Bview(i)(d) += Bf(d);
                    }
                }
            });
}

/**
 * @brief Applies the multipole field at position R to E and B
 *
 * @note Cannot be inside a kernel -> GPU incompatible
 *
 * @param R Position
 * @param E Electric field reference
 * @param B Magnetic field reference
 */
void Multipole::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    // Check bounds
    if (R(2) < 0.0 || R(2) >= getGeometry().getElementLength()) return;
    if (!ApertureHelper::isInsideAperture(R, aperture_m)) {
        return;
    }

    // Compute field
    computeFieldHost(R, E, B);
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
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double&,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    // Check bounds
    if (R(2) < 0.0 || R(2) >= getGeometry().getElementLength()) return false;
    if (!ApertureHelper::isInsideAperture(R, aperture_m)) {
        return true;
    }

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
void Multipole::accept(BeamlineVisitor& visitor) const { visitor.visitMultipole(*this); }

template <class ViewType>
KOKKOS_INLINE_FUNCTION void Multipole::computeField(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& B, const ViewType& normal,
        const ViewType& skew, int maxNormal, int maxSkew) {
    int count = (maxNormal > maxSkew) ? maxNormal : maxSkew;
    if (count > MAX_MP_ORDER) {
        count = MAX_MP_ORDER;
    }

    // re + i im = (x + iy)^n, starting at n = 0
    double re = 1.0, im = 0.0;
    for (int n = 0; n < count; ++n) {
        const double kn = (n < maxNormal) ? normal(n) : 0.0;
        const double ks = (n < maxSkew) ? skew(n) : 0.0;
        B(0) += kn * im - ks * re;
        B(1) += kn * re + ks * im;

        const double nextRe = re * R(0) - im * R(1);
        im                  = re * R(1) + im * R(0);
        re                  = nextRe;
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
        const Vector_t<double, 3> R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) const {
    const auto normal = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NormalComponents);
    const auto skew   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SkewComponents);
    computeField(R, B, normal, skew, max_NormalComponent_m, max_SkewComponent_m);
}
/**
 * @brief Setup, multipole goes online
 *
 * @param bunch The reference bunch
 */
void Multipole::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    online_m       = true;
}

void Multipole::finalise() { online_m = false; }

void Multipole::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

ElementType Multipole::getType() const { return ElementType::MULTIPOLE; }

// isInside() is inherited from ElementBase (field extent [0, L) + transverse aperture).

bool Multipole::isFocusing(int component) const {
    if (component < 0)
        throw GeneralOpalException("Multipole::isFocusing", "component negative");
    else if (static_cast<unsigned long>(component) >= NormalComponents.extent(0))
        throw GeneralOpalException("Multipole::isFocusing", "component too big");

    // Fix: Use getNormalComponent() to safely retrieve the value from GPU memory
    return getNormalComponent(component) * std::pow(-1, component + 1)
                   * RefPartBunch_m->getParticleContainer()->getChargePerParticle()
           > 0.0;
}

/* ========================================================================== */
/* =========================== Unused Functions ============================= */

// set the number of slices for map tracking
void Multipole::setNSlices(const std::size_t& nSlices) { nSlices_m = nSlices; }

// get the number of slices for map tracking
std::size_t Multipole::getNSlices() const { return nSlices_m; }
/* ========================================================================== */
