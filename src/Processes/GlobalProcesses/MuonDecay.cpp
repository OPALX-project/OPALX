#include "Processes/GlobalProcesses/MuonDecay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"

void MuonDecay::createDaughterParticles(
        std::size_t localDestroyNum, std::size_t oldDaughterLocal,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
        const Kokkos::View<double*>& parentDt,
        const Kokkos::View<ippl::Vector<float, 3>*>& parentPol) {
    using pc_size_type = ippl::detail::size_type;

    auto dR  = daughterPC_m->R.getView();
    auto dP  = daughterPC_m->P.getView();
    auto dDt = daughterPC_m->dt.getView();

    const auto pool           = randPool_m;
    const double daughterMass = daughterMassGeV_m;
    // V-A asymmetry sign: +1 for mu-, -1 for mu+. Set in Decay::apply from the
    // parent container's reference charge.
    const int chargeSign  = parentChargeSign_m;
    const double asymSign = (chargeSign < 0) ? +1.0 : -1.0;

    /* Kinematics of mu -> e + nu_e + nu_mu (three-body decay).
     *
     * The two neutrinos are unobserved, so the electron energy in the muon
     * rest frame is not fixed but follows the Michel spectrum:
     *   f(x) = 2 x^2 (3 - 2x),   x = E_e / E_max,   E_max = m_mu / 2.
     *
     * Procedure per daughter:
     *   1. Sample x from the Michel spectrum via rejection sampling.
     *   2. Compute daughter energy E = x * E_max and momentum |p| = sqrt(E^2 - m_e^2).
     *   3. Sample the emission direction in the parent rest frame from the
     *      polarization-dependent angular distribution
     *      1 + alpha(x) cos theta (theta measured from the parent's spin),
     *      then orient it relative to the parent's Pol vector. Reduces to
     *      isotropic when the parent is unpolarized.
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

                /* Step 3: spin-dependent emission direction in the muon rest frame.
                 *
                 * Full polarized Michel rate (V-A theory, m_e -> 0, no radiative
                 * corrections) as a joint density in reduced energy x and the
                 * emission angle theta from the muon spin, with c = cos theta:
                 *
                 *   dN/(dx dc) ~ x^2 [ (3 - 2x) + s |P| (2x - 1) c ],   s = asymSign.
                 *
                 * Marginal in x (integrate c over [-1, 1]; the c-term cancels by
                 * symmetry) is the Michel spectrum already sampled in Step 1:
                 *
                 *   p(x) ~ 2 x^2 (3 - 2x) = f(x).
                 *
                 * Conditional angle (joint / marginal) is linear in c:
                 *
                 *   p(c | x) = 1/2 [ 1 + alpha(x) c ],
                 *   alpha(x) = s |P| (2x - 1) / (3 - 2x),   |alpha| <= |P| <= 1.
                 *
                 * Inverse-CDF sampling: F(c) = 1/2 (1 + c) + alpha/4 (c^2 - 1).
                 * Setting F(c) = u (u ~ U[0,1]) and rearranging gives
                 *
                 *   (alpha/4) c^2 + (1/2) c + (1/2 - alpha/4 - u) = 0,
                 *
                 * solved below for the root in [-1, 1]. For alpha ~ 0 this reduces
                 * to the uniform c = 2u - 1 (separate branch, avoids dividing by
                 * the vanishing leading coefficient). The sampled (theta, phi) are
                 * then oriented with the polar axis along Pol-hat. */
                double polX = 0.0, polY = 0.0, polZ = 0.0;
                if (parentPol.extent(0) > 0) {
                    polX = static_cast<double>(parentPol(j)[0]);
                    polY = static_cast<double>(parentPol(j)[1]);
                    polZ = static_cast<double>(parentPol(j)[2]);
                }
                const double polMag = Kokkos::sqrt(polX * polX + polY * polY + polZ * polZ);

                const double alpha = (polMag > 0.0)
                                             ? asymSign * polMag * (2.0 * x - 1.0) / (3.0 - 2.0 * x)
                                             : 0.0;

                const double u = gen.drand(0.0, 1.0);
                double cosTheta;
                if (Kokkos::fabs(alpha) < 1.0e-12) {
                    cosTheta = 2.0 * u - 1.0;
                } else {
                    // Solve (alpha/4) c^2 + (1/2) c + (1/2 - alpha/4 - u) = 0 for c in [-1, 1].
                    const double aQ   = alpha / 4.0;
                    const double bQ   = 0.5;
                    const double cQ   = 0.5 - alpha / 4.0 - u;
                    const double disc = bQ * bQ - 4.0 * aQ * cQ;
                    const double sq   = Kokkos::sqrt(disc > 0.0 ? disc : 0.0);
                    const double r1   = (-bQ + sq) / (2.0 * aQ);
                    const double r2   = (-bQ - sq) / (2.0 * aQ);
                    // Pick the root inside [-1, 1] (numerically clamp).
                    const double pick = (r1 >= -1.0 && r1 <= 1.0) ? r1 : r2;
                    cosTheta          = pick < -1.0 ? -1.0 : (pick > 1.0 ? 1.0 : pick);
                }
                const double sinTheta = Kokkos::sqrt(
                        1.0 - cosTheta * cosTheta > 0.0 ? 1.0 - cosTheta * cosTheta : 0.0);
                const double phi = 2.0 * Physics::pi * gen.drand(0.0, 1.0);

                // Daughter momentum in the muon rest frame, expressed in the
                // spin-aligned local frame: z along Pol-hat (polar angle theta),
                // x/y an arbitrary transverse pair (azimuth phi). Rotated into
                // lab-aligned rest-frame axes just below.
                const double pLoc_x = pMag * sinTheta * Kokkos::cos(phi);
                const double pLoc_y = pMag * sinTheta * Kokkos::sin(phi);
                const double pLoc_z = pMag * cosTheta;

                // Rotate from Pol-aligned local frame (z = Pol-hat) into the lab-axis
                // rest frame. If |Pol| is 0 the choice is irrelevant.
                double pxRF, pyRF, pzRF;
                if (polMag > 1.0e-30) {
                    const double nz_x = polX / polMag;
                    const double nz_y = polY / polMag;
                    const double nz_z = polZ / polMag;
                    // Build an orthonormal basis {ex, ey, nz}. Choose ex perpendicular to nz.
                    double ex_x, ex_y, ex_z;
                    if (Kokkos::fabs(nz_x) <= Kokkos::fabs(nz_y)
                        && Kokkos::fabs(nz_x) <= Kokkos::fabs(nz_z)) {
                        const double n = Kokkos::sqrt(nz_y * nz_y + nz_z * nz_z);
                        ex_x           = 0.0;
                        ex_y           = -nz_z / n;
                        ex_z           = nz_y / n;
                    } else if (Kokkos::fabs(nz_y) <= Kokkos::fabs(nz_z)) {
                        const double n = Kokkos::sqrt(nz_x * nz_x + nz_z * nz_z);
                        ex_x           = nz_z / n;
                        ex_y           = 0.0;
                        ex_z           = -nz_x / n;
                    } else {
                        const double n = Kokkos::sqrt(nz_x * nz_x + nz_y * nz_y);
                        ex_x           = -nz_y / n;
                        ex_y           = nz_x / n;
                        ex_z           = 0.0;
                    }
                    // ey = nz x ex
                    const double ey_x = nz_y * ex_z - nz_z * ex_y;
                    const double ey_y = nz_z * ex_x - nz_x * ex_z;
                    const double ey_z = nz_x * ex_y - nz_y * ex_x;
                    pxRF              = pLoc_x * ex_x + pLoc_y * ey_x + pLoc_z * nz_x;
                    pyRF              = pLoc_x * ex_y + pLoc_y * ey_y + pLoc_z * nz_y;
                    pzRF              = pLoc_x * ex_z + pLoc_y * ey_z + pLoc_z * nz_z;
                } else {
                    pxRF = pLoc_x;
                    pyRF = pLoc_y;
                    pzRF = pLoc_z;
                }

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
