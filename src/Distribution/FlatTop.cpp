#include "FlatTop.h"
#include <cmath>
#include <memory>
#include "Distribution.h"
#include "SamplingBase.hpp"

using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;
using Dist_t        = ippl::random::NormalDistribution<double, 3>;

FlatTop::FlatTop(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist), rand_pool_m(determineRandInit()) {
    setParameters(opalDist);
}

FlatTop::FlatTop(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
        double tPulseLengthFWHM, Vector_t<double, 3> sigmaR)
    : SamplingBase(pc, fc), rand_pool_m(determineRandInit()) {
    setInternalVariables(emitting, sigmaTFall, sigmaTRise, cutoff, tPulseLengthFWHM, sigmaR);
}

FlatTop::FlatTop(
        std::shared_ptr<ParticleContainer_t> pc, bool emitting, double sigmaTFall,
        double sigmaTRise, Vector_t<double, 3> cutoff, double tPulseLengthFWHM,
        Vector_t<double, 3> sigmaR)
    : SamplingBase(pc), rand_pool_m(determineRandInit()) {
    setInternalVariables(emitting, sigmaTFall, sigmaTRise, cutoff, tPulseLengthFWHM, sigmaR);
}

void FlatTop::setWithDomainDecomp(bool withDomainDecomp) { withDomainDecomp_m = withDomainDecomp; }

size_t FlatTop::determineRandInit() {
    extern Inform* gmsg;
    size_t randInit;
    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    }
    return randInit;
}

void FlatTop::setParameters(Distribution_t* opalDist) {
    setInternalVariables(
            opalDist->emitting_m, opalDist_m->getSigmaTFall(), opalDist_m->getSigmaTRise(),
            opalDist_m->getCutoffR(), opalDist->getTPulseLengthFWHM(), opalDist_m->getSigmaR());

    opalDist_m->setTEmission(emissionTime_m);

    // make sure only z direction is decomposed
    fc_m->setDecomp({false, false, true});
}

void FlatTop::setInternalVariables(
        bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
        double tPulseLengthFWHM, Vector_t<double, 3> sigmaR) {
    emitting_m = emitting;
    // time span of fall is [0, riseTime, riseTime+flattopTime, fallTime+flattopTime+riseTime ]
    sigmaTFall_m = sigmaTFall;
    sigmaTRise_m = sigmaTRise;
    cutoffR_m    = cutoff;

    fallTime_m = sigmaTFall_m * cutoffR_m[2];  // fall is [0, fallTime]
    flattopTime_m =
            tPulseLengthFWHM - std::sqrt(2.0 * std::log(2.0)) * (sigmaTRise_m + sigmaTFall_m);
    if (flattopTime_m < 0.0) {
        flattopTime_m = 0.0;
    }
    riseTime_m = sigmaTRise_m * cutoffR_m[2];

    emissionTime_m = fallTime_m + flattopTime_m + riseTime_m;

    // These expression are taken from the old OPAL
    // I think normalizedFlankArea is int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sigma
    // Instead of int_0^{cutoff} exp(-(x/sigma)^2/2 ) / sqrt(2*pi) / sigma, which is strange!
    // So the distribution of tails are exp(-(x/sigma)^2/2 ) and not Gaussian!
    normalizedFlankArea_m =
            0.5 * std::sqrt(Physics::two_pi) * std::erf(cutoffR_m[2] / std::sqrt(2.0));
    distArea_m = flattopTime_m + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea_m;

    sigmaR_m = sigmaR;
}

void FlatTop::generateUniformDisk(size_type nlocal, size_t nNew, double dt) {
    if (nNew == 0) {
        return;
    }

    GeneratorPool rand_pool = rand_pool_m;
    view_type Rview         = pc_m->R.getView();
    view_type Pview         = pc_m->P.getView();
    auto dtView             = pc_m->dt.getView();

    double pi                  = Physics::pi;
    Vector_t<double, 3> sigmaR = sigmaR_m;
    Vector_t<double, 3> R0     = R0_m;
    Vector_t<double, 3> P0     = P0_m;

    auto range = Kokkos::RangePolicy<>(nlocal, nlocal + nNew);

    // Position sampling is shared: uniform on elliptical disk + R0 offset.
    Kokkos::parallel_for(
            "unitDisk_R", range, KOKKOS_LAMBDA(const size_t j) {
                auto generator = rand_pool.get_state();
                double r       = Kokkos::sqrt(generator.drand(0., 1.));
                double theta   = 2.0 * pi * generator.drand(0., 1.);
                double frac    = generator.drand(0., 1.);
                rand_pool.free_state(generator);

                Rview(j)[0] = r * Kokkos::cos(theta) * sigmaR[0] + R0[0];
                Rview(j)[1] = r * Kokkos::sin(theta) * sigmaR[1] + R0[1];
                Rview(j)[2] = 0.0 + R0[2];

                // Each particle is assigned a fractional timestep dt_i = f * dt where f ~ U(0,1).
                // This represents the fraction of the next integration step the particle will
                // experience, as if the particle were born at a random time within [t, t+dt]. The
                // per-particle dt is used by the Boris integrator (push/kick) and by
                // scaleDtByCharge for field deposition, so the fractional dt naturally spreads
                // particles in z and gives fractional charge contribution without needing to sample
                // Rz explicitly.
                dtView(j) = frac * dt;
            });

    // Momentum sampling depends on the emission model chosen in EMISSIONSOURCE.
    // BEAM reference momentum (avrgpz) is intentionally NOT applied:
    // for emitted beams only the thermal energy matters (old OPAL behavior).
    if (emissionModel_m == "ASTRA") {
        // ASTRA: 3D isotropic thermal emission on forward half-sphere.
        const double pTot = euclidean_norm(P0);
        Kokkos::parallel_for(
                "unitDisk_P_astra", range, KOKKOS_LAMBDA(const size_t j) {
                    auto generator = rand_pool.get_state();
                    double rand1   = generator.drand(0., 1.);
                    double rand2   = generator.drand(0., 1.);
                    rand_pool.free_state(generator);

                    double phi   = 2.0 * Kokkos::acos(Kokkos::sqrt(rand1));
                    double theta = 2.0 * pi * rand2;

                    Pview(j)[0] = pTot * Kokkos::sin(phi) * Kokkos::cos(theta);
                    Pview(j)[1] = pTot * Kokkos::sin(phi) * Kokkos::sin(theta);
                    Pview(j)[2] = pTot * Kokkos::fabs(Kokkos::cos(phi));
                });
    } else {
        // NONE: all "thermal" momentum applied in z direction.
        Kokkos::parallel_for(
                "unitDisk_P_none", range, KOKKOS_LAMBDA(const size_t j) { Pview(j) = P0; });
    }
    Kokkos::fence();

    // Correct z for the missed timeIntegration1 half-push.
    //
    // New particles are born between timeIntegration1 and timeIntegration2.
    // They only receive the timeIntegration2 push: Δz = 0.5·β_z·c·(frac·dt).
    // The inter-batch spacing is β_z·c·dt (a full step). Without correction the
    // batch spread is 0.5·β_z·c·dt — half the spacing — leaving 50% gaps that
    // appear as visible "discs" under strong acceleration.
    //
    // Adding the equivalent of the missed timeIntegration1 half-push at birth:
    //   z_init += 0.5·β_z(birth)·c·(frac·dt)
    // makes the total spread β_z·c·frac·dt ∈ [0, β_z·c·dt], exactly tiling the
    // inter-batch gap. The approximation is that β_z is evaluated at birth
    // momentum; residual error is O(F·dt²) (one-timestep field correction).
    const double c = Physics::c;
    Kokkos::parallel_for(
            "unitDisk_Zcorr", range, KOKKOS_LAMBDA(const size_t j) {
                double px     = Pview(j)[0];
                double py     = Pview(j)[1];
                double pz     = Pview(j)[2];
                double gamma  = Kokkos::sqrt(1.0 + px * px + py * py + pz * pz);
                double beta_z = pz / gamma;
                Rview(j)[2] += 0.5 * beta_z * c * dtView(j);
            });
    Kokkos::fence();
}

void FlatTop::setNr(Vector_t<double, 3> nr) { nr_m = nr; }

void FlatTop::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) {
    setNr(nr);

    // initial allocation is similar for both emitting and non-emitting cases
    allocateParticles(numberOfParticles);

    // This doesn't do anything, since the allocation is now done in TrackRun.
    /*if(emitting_m){
        // set nlocal to 0 for the very first time step, before sampling particles
        //pc_m->setLocalNum(0);

        Kokkos::View<bool*> tmp_invalid("tmp_invalid", 0);
        // \todo might be abuse of semantics: maybe think about new pc_m->setTotalNum or
    pc_m->updateTotal function instead? pc_m->destroy(tmp_invalid, pc_m->getLocalNum()); } else {
        throw OpalException(
            "FlatTop::generateParticles",
            "FlatTop does not support generating particles without emitting. Set EMITTED = true in "
            "the input file.");
    }*/
}

double FlatTop::FlatTopProfile(double t) {
    double t0;
    if (t < riseTime_m) {
        t0 = riseTime_m;
        return exp(-pow((t - t0) / sigmaTRise_m, 2) / 2.);
        //  In the old opal, tails seem to be exp(-x^2/sigma^2/2) rather than Gaussian with
        //  normalizing factor.
    } else if (t > riseTime_m && t < riseTime_m + flattopTime_m) {
        return 1.;
    } else if (t > riseTime_m + flattopTime_m && t < fallTime_m + flattopTime_m + riseTime_m) {
        t0 = fallTime_m + flattopTime_m;
        return exp(-pow((t - t0) / sigmaTFall_m, 2) / 2.);
        //  In the old opal, tails seem to be exp(-x^2/sigma^2/2) rather than Gaussian with
        //  normalizing factor.
    } else
        return 0.;
}

size_t FlatTop::computeNlocalUniformly(size_t nglobal) {
    // Use ippl::Comm so we match the communicator used for allocation in TrackRun (maxLocalNum).
    const int nranks = ippl::Comm->size();
    const int rank   = ippl::Comm->rank();

    const size_t nranks_u  = static_cast<size_t>(nranks > 0 ? nranks : 1);
    const size_type nlocal = nglobal / nranks_u;
    const size_t remaining = nglobal - nlocal * nranks_u;

    // First 'remaining' ranks get one extra particle.
    if (remaining > 0 && static_cast<size_t>(rank) < remaining) {
        return nlocal + 1;
    }
    return nlocal;
}

double FlatTop::integrateTrapezoidal(double x1, double x2, double y1, double y2) {
    return 0.5 * (y1 + y2) * fabs(x2 - x1);
}

void FlatTop::initDomainDecomp(double BoxIncr) {
    auto* mesh                 = &fc_m->getMesh();
    auto* FL                   = &fc_m->getFL();
    Vector_t<double, 3> sigmaR = sigmaR_m;
    ippl::Vector<double, 3> o;
    ippl::Vector<double, 3> e;
    double tol = 1e-15;  // enlarge grid by tol to avoid missing particles on boundaries
    o[0]       = -sigmaR[0] - tol;
    e[0]       = sigmaR[0] + tol;
    o[1]       = -sigmaR[1] - tol;
    e[1]       = sigmaR[1] + tol;
    o[2]       = 0.0 - tol;
    e[2]       = Physics::c * emissionTime_m + tol;

    ippl::Vector<double, 3> l = e - o;
    hr_m                      = (1.0 + BoxIncr / 100.) * (l / nr_m);
    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o - 0.5 * hr_m * BoxIncr / 100.);
    pc_m->getLayout().updateLayout(*FL, *mesh);
}

FlatTop::size_type FlatTop::countEnteringParticlesPerRank(double t0, double tf) {
    const double tArea = integrateTrapezoidal(t0, tf, FlatTopProfile(t0), FlatTopProfile(tf));
    // Integer-safe: floor to avoid over-counting, then clamp to totalN_m (avoids fp overflow into
    // size_type)
    const double totalNewD = std::floor(static_cast<double>(totalN_m) * tArea / distArea_m);
    size_type totalNew     = static_cast<size_type>(totalNewD);
    if (totalNew > totalN_m) {
        totalNew = totalN_m;
    }

    return pc_m ? static_cast<size_type>(computeLocalEmitCount(static_cast<size_t>(totalNew)))
                : computeNlocalUniformly(totalNew);

    // The following code is unnecessarily complicated in my opinion as a load balancer will be
    // called anyways. Apart from that it will mitigate errors in the allocated number of particles
    // per rank.
    /*
    if (totalNew > 0) {
        if(!withDomainDecomp_m){
            // the same number of particles per rank
            nlocalNew = computeNlocalUniformly(totalNew);
        }
        else{
            // select number of particles per rank using estimated domain decomposition at final
    emission time
            // find min/max of particle positions for [t,t+dt]
            Vector_t<double, 3> prmin, prmax;
            Vector_t<double, 3> sigmaR = sigmaR_m;
            prmin[0] = -sigmaR[0];
            prmax[0] =  sigmaR[0];
            prmin[1] = -sigmaR[1];
            prmax[1] =  sigmaR[1];
            prmin[2] = std::max(0.0, Physics::c * t0);
            prmax[2] = Physics::c*tf;

            double dx = prmax[0] - prmin[0];
            double dy = prmax[1] - prmin[1];
            double dz = prmax[2] - prmin[2];

            if (dx <= 0 || dy <= 0 || dz <= 0) {
                throw std::runtime_error("Invalid global particle volume: prmax must be greater than
    prmin.");
            }

            double globalpvolume = dx * dy * dz;

            // find min/max of subdomains for the current rank
            auto regions = pc_m->getLayout().getRegionLayout().gethLocalRegions();
            int rank = ippl::Comm->rank();

            if (rank < 0 || static_cast<size_t>(rank) >= regions.size()) {
                throw std::runtime_error("Invalid rank index in gethLocalRegions()");
            }

            Vector_t<double, 3> locrmin, locrmax;
            for (unsigned d = 0; d < Dim; ++d) {
               locrmax[d] = regions(rank)[d].max();
               locrmin[d] = regions(rank)[d].min();
            }

            if (prmax[0] >= locrmin[0] && prmin[0] <= locrmax[0] &&
                prmax[1] >= locrmin[1] && prmin[1] <= locrmax[1] &&
                prmax[2] >= locrmin[2] && prmin[2] <= locrmax[2]) {

                double x1 = std::max(prmin[0], locrmin[0]);
                double x2 = std::min(prmax[0], locrmax[0]);
                double y1 = std::max(prmin[1], locrmin[1]);
                double y2 = std::min(prmax[1], locrmax[1]);
                double z1 = std::max(prmin[2], locrmin[2]);
                double z2 = std::min(prmax[2], locrmax[2]);

                if (x2 >= x1 && y2 >= y1 && z2 >= z1) {
                    const double locpvolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
                    if (globalpvolume > 0) {
                        const double frac = locpvolume / globalpvolume;
                        size_type nFromFrac =
    static_cast<size_type>(std::floor(static_cast<double>(totalNew) * frac)); nlocalNew = (nFromFrac
    > totalNew) ? totalNew : nFromFrac; } else { nlocalNew = 0;
                    }
                } else {
                    nlocalNew = 0;
                }
            }
        }
    }*/
}

void FlatTop::allocateParticles(size_t numberOfParticles) {
    totalN_m = numberOfParticles;

    // Initial allocation is now handled centrally in TrackRun / PartBunch via the
    // bunch's total particle count. Here we only record the desired total number
    // of particles for this distribution. Actual per-step emission will append
    // new particles using pc_m->create(nNew), guarded by a global BEAM::NALLOC
    // limit in ParallelTracker.
}

void FlatTop::emitParticles(double t, double dt) {
    Inform msgAll = Inform("FlatTop::emitParticles", INFORM_ALL_NODES);

    // Time profile uses (t - t0) so sampling begins at t0.
    double tShift  = t - t0_m;
    double dtShift = dt;
    // Count number of new particles to be emitted in [tShift, tShift+dtShift].
    size_type nNew = countEnteringParticlesPerRank(tShift, tShift + dtShift);
    msgAll << level5 << "* " << nNew << " new particles to be emitted" << endl;
    // if (nNew == 0) { return; } // don't return, see why below

    // Current number of particles per rank.
    size_type nlocal = pc_m->getLocalNum();
    msgAll << level5 << "* " << nlocal << " particles already in the container" << endl;

    // Intended design: fill only into pre-allocated slice [nextEmitIndex, nextEmitIndex+nNew)
    // (no create), when ParticleContainer supports reserve(total) + setActiveSize. Until then
    // we extend the container each step.
    // Note that create can be safely called with nNew=0 (no-op), but it NEEDS to be called, since
    // other ranks might have != 0 leading to a MPI_Allreduce.
    pc_m->create(nNew);

    // Generate new particles on uniform disc (sample into [nlocal, nlocal+nNew)).
    // Each particle receives a fractional per-particle dt for sub-timestep spreading.
    msgAll << level3 << "* generate particles on a disc" << endl;
    generateUniformDisk(nlocal, nNew, dt);

    pc_m->markMomentsDirty();

    msgAll << level3 << "* " << nNew << " new particles emmitted" << endl;
}

void FlatTop::testNumEmitParticles(size_type nsteps, double dt) {
    size_type nNew;
    int rank, numRanks;
    double t = 0.0;
    double c = Physics::c;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    std::string filename = "timeNpart_" + std::to_string(rank) + ".txt";
    std::ofstream file(filename);

    // Check if the file opened successfully
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    for (size_type i = 0; i < nsteps; i++) {
        nNew = countEnteringParticlesPerRank(t, t + dt);

        // current number of particles per rank
        size_type nlocal = pc_m->getLocalNum();

        // extend particle container to accomodate new particles
        pc_m->create(nNew);

        // generate new particles on uniform disc
        generateUniformDisk(nlocal, nNew, dt);

        // write to a file
        auto rViewDevice = pc_m->R.getView();
        auto rView       = Kokkos::create_mirror_view(rViewDevice);
        Kokkos::deep_copy(rView, rViewDevice);

        for (size_type j = nlocal; j < nlocal + nNew; j++) {
            file << t << " " << (emissionTime_m - t) * c << " " << rView(j)[0] << " " << rView(j)[1]
                 << "\n";
        }
        file.flush();  // Ensure data is written to disk

        t = t + dt;
    }
    file.close();
    ippl::Comm->barrier();
}

void FlatTop::testEmitParticles(size_type nsteps, double dt) {
    double t = 0.0;

    for (size_type i = 0; i < nsteps; i++) {
        emitParticles(t, dt);

        t = t + dt;
    }
}
