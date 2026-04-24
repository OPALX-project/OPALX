#include "Processes/GlobalProcesses/PionDecay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"

void PionDecay::createDaughterParticles(
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

    /* Kinematics of pi -> mu + nu_mu (two-body decay).
     *
     * With one massless neutrino, energy-momentum conservation in the pion
     * rest frame fixes the muon momentum uniquely:
     *   E_mu = (M_pi^2 + m_mu^2) / (2 M_pi)
     *   |p|  = (M_pi^2 - m_mu^2) / (2 M_pi)
     *
     * The direction is isotropic (pion is spin-0).
     *
     * Procedure per daughter:
     *   1. Compute the fixed rest-frame momentum (done once, outside the loop).
     *   2. Pick an isotropic direction in the parent rest frame.
     *   3. Lorentz-boost the rest-frame 4-momentum to the lab frame using the
     *      parent's beta*gamma vector.
     *   4. Store the result as beta*gamma = p_lab / m_daughter. */
    const double M         = parentMassGeV_m;
    const double md        = daughterMassGeV_m;
    const double pFixed    = (M * M - md * md) / (2.0 * M);
    const double eDaughter = Kokkos::sqrt(pFixed * pFixed + md * md);

    Kokkos::parallel_for(
            "PionDecay::createDaughters", static_cast<pc_size_type>(localDestroyNum),
            KOKKOS_LAMBDA(const pc_size_type j) {
                auto gen = pool.get_state();

                /* Step 2: Isotropic direction in parent rest frame. */
                const double cosTheta = 2.0 * gen.drand(0.0, 1.0) - 1.0;
                const double sinTheta = Kokkos::sqrt(1.0 - cosTheta * cosTheta);
                const double phi      = 2.0 * Physics::pi * gen.drand(0.0, 1.0);
                const double pxRF     = pFixed * sinTheta * Kokkos::cos(phi);
                const double pyRF     = pFixed * sinTheta * Kokkos::sin(phi);
                const double pzRF     = pFixed * cosTheta;

                /* Step 3: Lorentz boost to lab frame.
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

                /* Step 4: Write daughter attributes.
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
