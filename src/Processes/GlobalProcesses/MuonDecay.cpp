#include "Processes/GlobalProcesses/MuonDecay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"

void MuonDecay::createDaughterParticles(
        std::size_t localDestroyNum, std::size_t oldDaughterLocal,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
        const Kokkos::View<double*>& parentDt) {
    using pc_size_type = ippl::detail::size_type;

    auto dR  = daughterPC_m->R.getView();
    auto dP  = daughterPC_m->P.getView();
    auto dDt = daughterPC_m->dt.getView();

    const auto pool           = randPool_m;
    const double daughterMass = daughterMassGeV_m;

    /* Kinematics of mu -> e + nu_e + nu_mu (three-body decay).
     *
     * The two neutrinos are unobserved, so the electron energy in the muon
     * rest frame is not fixed but follows the Michel spectrum:
     *   f(x) = 2 x^2 (3 - 2x),   x = E_e / E_max,   E_max = m_mu / 2.
     *
     * Procedure per daughter:
     *   1. Sample x from the Michel spectrum via rejection sampling.
     *   2. Compute daughter energy E = x * E_max and momentum |p| = sqrt(E^2 - m_e^2).
     *   3. Pick an isotropic direction in the parent rest frame.
     *   4. Lorentz-boost the rest-frame 4-momentum to the lab frame using the
     *      parent's beta*gamma vector.
     *   5. Store the result as beta*gamma = p_lab / m_daughter. */
    const double eMax = parentMassGeV_m / 2.0;
    const double xMin = daughterMassGeV_m / eMax;
    const double fMax = 2.0;  // Michel spectrum upper bound at x = 1

    Kokkos::parallel_for(
            "MuonDecay::createDaughters", static_cast<pc_size_type>(localDestroyNum),
            KOKKOS_LAMBDA(const pc_size_type j) {
                auto gen = pool.get_state();

                /* Step 1: Sample reduced energy x from Michel spectrum (rejection). */
                double x;
                do {
                    x = xMin + (1.0 - xMin) * gen.drand(0.0, 1.0);
                } while (gen.drand(0.0, fMax) >= 2.0 * x * x * (3.0 - 2.0 * x));

                /* Step 2: Daughter energy and momentum magnitude in rest frame. */
                const double eDaughter = x * eMax;  // [GeV]
                const double pMag      = Kokkos::sqrt(
                        eDaughter * eDaughter - daughterMass * daughterMass);  // [GeV/c]

                /* Step 3: Isotropic direction in parent rest frame. */
                const double cosTheta = 2.0 * gen.drand(0.0, 1.0) - 1.0;
                const double sinTheta = Kokkos::sqrt(1.0 - cosTheta * cosTheta);
                const double phi      = 2.0 * Physics::pi * gen.drand(0.0, 1.0);
                const double pxRF     = pMag * sinTheta * Kokkos::cos(phi);
                const double pyRF     = pMag * sinTheta * Kokkos::sin(phi);
                const double pzRF     = pMag * cosTheta;

                /* Step 4: Lorentz boost to lab frame.
                 * Parent momentum is stored as beta*gamma (dimensionless).
                 * beta = P / gamma, with gamma = sqrt(1 + |P|^2).
                 * The boost formula for the spatial components is:
                 *   p_lab = p_RF + [(gamma-1)(beta . p_RF)/beta^2 + gamma*E] * beta */
                const auto& parentMom = parentP(j);
                const double bg2      = parentMom[0] * parentMom[0] + parentMom[1] * parentMom[1]
                                   + parentMom[2] * parentMom[2];
                const double gamma = Kokkos::sqrt(1.0 + bg2);
                const double betaX = parentMom[0] / gamma;
                const double betaY = parentMom[1] / gamma;
                const double betaZ = parentMom[2] / gamma;
                const double beta2 = bg2 / (gamma * gamma);

                double pxLab, pyLab, pzLab;
                if (beta2 > 0.0) {
                    const double betaDotP = betaX * pxRF + betaY * pyRF + betaZ * pzRF;
                    const double factor   = (gamma - 1.0) * betaDotP / beta2 + gamma * eDaughter;
                    pxLab                 = pxRF + factor * betaX;
                    pyLab                 = pyRF + factor * betaY;
                    pzLab                 = pzRF + factor * betaZ;
                } else {
                    pxLab = pxRF;
                    pyLab = pyRF;
                    pzLab = pzRF;
                }

                /* Step 5: Write daughter attributes.
                 * Position inherits from parent. Momentum stored as beta*gamma = p / m. */
                const pc_size_type idx = static_cast<pc_size_type>(oldDaughterLocal) + j;
                dR(idx)                = parentR(j);
                dP(idx)[0]             = pxLab / daughterMass;
                dP(idx)[1]             = pyLab / daughterMass;
                dP(idx)[2]             = pzLab / daughterMass;
                dDt(idx)               = parentDt(j);

                pool.free_state(gen);
            });
    Kokkos::fence();
}
