#include "Processes/GlobalProcesses/PionDecay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"

void PionDecay::createDaughterParticles(
        std::size_t localDestroyNum, std::size_t oldDaughterLocal,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
        const Kokkos::View<double*>& parentDt,
        const Kokkos::View<ippl::Vector<double, 3>*>& /*parentPol*/) {
    using pc_size_type = ippl::detail::size_type;

    auto dR  = daughterPC_m->R.getView();
    auto dP  = daughterPC_m->P.getView();
    auto dDt = daughterPC_m->dt.getView();

    // Daughter muon polarization view (empty unallocated View when daughter has no Pol).
    const bool daughterHasSpin = daughterPC_m->hasSpin();
    using PolView_t            = Kokkos::View<ippl::Vector<double, 3>*>;
    PolView_t dPol;
    if (daughterHasSpin) {
        dPol = daughterPC_m->Pol.getView();
    }

    // V-A: in the pion rest frame the muon is fully helicity-polarized.
    //   pi-  ->  mu-  +  anti-nu_mu :  mu- spin parallel to its momentum  (h = +1)
    //   pi+  ->  mu+  +  nu_mu      :  mu+ spin anti-parallel to momentum (h = -1)
    // parentChargeSign_m is set in Decay::apply from the parent reference.
    const double helicitySign = (parentChargeSign_m < 0) ? +1.0 : -1.0;

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

                // Daughter muon spin. In the pion rest frame the spin 4-vector is h (p/m, (E/m) n),
                // n the muon direction. Boost it to the lab, then take the spin
                // in the muon rest frame reached from the lab by one pure boost (what Pol holds).
                if (daughterHasSpin) {
                    const double sT = helicitySign * pFixed / md;
                    const double sN = helicitySign * eDaughter / (md * pFixed);
                    double sTLab    = sT;
                    double sxLab    = sN * pxRF;
                    double syLab    = sN * pyRF;
                    double szLab    = sN * pzRF;
                    if (beta2 > 0.0) {
                        const double betaDotS = betaX * sxLab + betaY * syLab + betaZ * szLab;
                        const double factor   = (gamma - 1.0) * betaDotS / beta2 + gamma * sT;
                        sTLab                 = gamma * (sT + betaDotS);
                        sxLab += factor * betaX;
                        syLab += factor * betaY;
                        szLab += factor * betaZ;
                    }
                    // Spin in the muon rest frame: s = S - S^0 / (gamma_mu + 1) u, u = p_lab / m.
                    const double gammaMu = Kokkos::sqrt(
                            1.0 + dP(idx)[0] * dP(idx)[0] + dP(idx)[1] * dP(idx)[1]
                            + dP(idx)[2] * dP(idx)[2]);
                    const double k = sTLab / (gammaMu + 1.0);
                    dPol(idx)[0]   = sxLab - k * dP(idx)[0];
                    dPol(idx)[1]   = syLab - k * dP(idx)[1];
                    dPol(idx)[2]   = szLab - k * dP(idx)[2];
                }

                pool.free_state(gen);
            });
    Kokkos::fence();
}
