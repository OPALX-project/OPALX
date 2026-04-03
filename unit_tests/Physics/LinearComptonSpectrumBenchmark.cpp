#include "LinearComptonSpectrumCommon.h"
#include "Distribution/MultiVariateGaussian.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Ippl.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
/**
 * @brief Command-line options for the standalone linear-Compton benchmark driver.
 *
 * The executable supports three validation modes:
 *
 * - deterministic single-electron spectra,
 * - sampled single-electron spectra,
 * - sampled finite-beam spectra driven by `MultiVariateGaussian`.
 *
 * The finite-beam parameters are intentionally simple and map directly onto the
 * Gaussian phase-space widths used in the benchmark, not onto the full OPALX
 * `BEAM` command semantics.
 */
struct CliOptions {
    std::filesystem::path output = "opalx-linear-compton-90deg-xi029-spectrum.csv";
    bool sampled = false;
    bool angular = false;
    bool finiteBeam = false;
    std::size_t samples = 250000;
    std::size_t beamParticles = 250000;
    int seed = 13579;
    double beamSigmaTransverse_m = 1.0e-3;
    double beamSigmaLongitudinal_m = 0.0;
    double beamRmsAngle_rad = 1.0e-3;
    double beamRelativeEnergySpread = 0.0;
};

/**
 * @brief Parse benchmark CLI arguments into a validated option bundle.
 *
 * The parser only performs syntactic checks here, such as missing values after
 * scalar options. Physical validation of beam sizes, angular spread, and
 * energy spread is done inside the benchmark routines so the error messages can
 * stay tied to the specific benchmark mode.
 */
CliOptions parseArguments(int argc, char** argv) {
    CliOptions options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--sampled") {
            options.sampled = true;
        } else if (arg == "--angular") {
            options.angular = true;
        } else if (arg == "--finite-beam") {
            options.finiteBeam = true;
        } else if (arg == "--samples") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --samples");
            }
            options.samples = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--beam-particles") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --beam-particles");
            }
            options.beamParticles = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--seed") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --seed");
            }
            options.seed = std::stoi(argv[++i]);
        } else if (arg == "--beam-sigma-transverse") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --beam-sigma-transverse");
            }
            options.beamSigmaTransverse_m = std::stod(argv[++i]);
        } else if (arg == "--beam-sigma-longitudinal") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --beam-sigma-longitudinal");
            }
            options.beamSigmaLongitudinal_m = std::stod(argv[++i]);
        } else if (arg == "--beam-rms-angle") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --beam-rms-angle");
            }
            options.beamRmsAngle_rad = std::stod(argv[++i]);
        } else if (arg == "--beam-relative-energy-spread") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --beam-relative-energy-spread");
            }
            options.beamRelativeEnergySpread = std::stod(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: LinearComptonSpectrumBenchmark [output.csv] [--angular] [--sampled] [--finite-beam] [--samples N] [--beam-particles N] [--seed S]\n"
                << "  default mode             : deterministic single-electron energy spectrum\n"
                << "  --angular                : write the lab polar-angle histogram instead of the energy histogram\n"
                << "  --sampled                : host-only Monte Carlo benchmark using LinearCompton::sampleEvent\n"
                << "  --finite-beam            : sample a finite-emittance electron beam with MultiVariateGaussian\n"
                << "  --samples N              : number of Monte Carlo samples in single-electron sampled mode\n"
                << "  --beam-particles N       : number of sampled beam macroparticles in finite-beam mode\n"
                << "  --beam-sigma-transverse  : transverse RMS beam size in m (default: 1e-3)\n"
                << "  --beam-sigma-longitudinal: longitudinal RMS beam size in m (default: 0.0, internally clamped)\n"
                << "  --beam-rms-angle         : transverse RMS beam angle in rad (default: 1e-3)\n"
                << "  --beam-relative-energy-spread: RMS relative energy spread sigma_E / E (default: 0.0)\n"
                << "  --seed S                 : Options::seed value used in sampled modes\n";
            std::exit(0);
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("Unknown option: " + arg);
        } else {
            options.output = arg;
        }
    }

    return options;
}

/**
 * @brief Reserve enough particle capacity for the finite-beam benchmark.
 *
 * `MultiVariateGaussian` inserts particles into an existing `ParticleContainer`.
 * The benchmark keeps allocation behavior deterministic by creating a local
 * capacity estimate up front, then clearing the initially created slots again.
 *
 * @param pc Particle container used by the finite-beam benchmark.
 * @param totalParticles Global benchmark macroparticle count.
 */
void preallocateParticleCapacity(const std::shared_ptr<ParticleContainer<double, 3>>& pc,
                                 std::size_t totalParticles) {
    const int nranks = std::max(1, ippl::Comm->size());
    const std::size_t nranksU = static_cast<std::size_t>(nranks);
    const std::size_t maxLocalNum = totalParticles / nranksU + 2 * nranksU + 1;

    pc->create(maxLocalNum);
    Kokkos::View<bool*> tmpInvalid("tmp_invalid", maxLocalNum);
    pc->destroy(tmpInvalid, maxLocalNum);
}

/**
 * @brief Create a small temporary particle container for finite-beam sampling.
 *
 * The benchmark is host-side and does not use the full OPALX beam pipeline.
 * It therefore creates a minimal particle container with a mesh domain scaled
 * to the requested RMS beam sizes. The longitudinal extent is clamped away from
 * zero so that the mesh remains well-defined even for the `sigma_z = 0`
 * benchmark cases.
 *
 * @param sigmaR RMS position widths used to size the benchmark mesh.
 * @return Temporary particle container for Gaussian electron sampling.
 */
std::shared_ptr<ParticleContainer<double, 3>> makeParticleContainer(const ippl::Vector<double, 3>& sigmaR) {
    ippl::Vector<int, 3> nr = 32;
    const double sigmaScaleXY = std::max(8.0 * sigmaR[0], 1.0e-2);
    const double sigmaScaleZ = std::max(8.0 * sigmaR[2], 1.0e-6);

    ippl::Vector<double, 3> rmin;
    rmin[0] = -sigmaScaleXY;
    rmin[1] = -sigmaScaleXY;
    rmin[2] = -sigmaScaleZ;
    ippl::Vector<double, 3> rmax;
    rmax[0] = sigmaScaleXY;
    rmax[1] = sigmaScaleXY;
    rmax[2] = sigmaScaleZ;
    ippl::Vector<double, 3> origin = rmin;
    ippl::Vector<double, 3> hr = (rmax - rmin) / ippl::Vector<double, 3>(nr);
    std::array<bool, 3> decomp = {true, true, true};

    ippl::NDIndex<3> domain;
    for (unsigned i = 0; i < 3; ++i) {
        domain[i] = ippl::Index(nr[i]);
    }

    auto fc = std::make_shared<FieldContainer_t>(hr, rmin, rmax, decomp, domain, origin, true);
    (void)fc;

    Mesh_t<3> mesh(domain, hr, origin);
    FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);
    return std::make_shared<ParticleContainer<double, 3>>(mesh, fl);
}

/**
 * @brief Sample the finite-beam photon energy spectrum with one Compton event per electron.
 *
 * Electrons are drawn from a six-dimensional Gaussian phase-space model using
 * `MultiVariateGaussian`. For each sampled electron, the local momentum vector
 * defines both the incoming direction and total energy used to build a
 * `Physics::LinearCompton::SamplingKernel`. The benchmark then draws exactly
 * one unpolarized Compton event and bins the outgoing photon energy in the lab
 * frame.
 *
 * The relative energy spread option is mapped into a longitudinal momentum RMS
 * by linearizing around the ultra-relativistic reference state,
 * sigma_pz ~= sigma_E / beta_0. This is appropriate for the
 * current weak-field, high-energy validation regime but is not yet a general
 * beam-dynamics model.
 */
LinearComptonBenchmark::SpectrumHistogram sampleFiniteBeamEnergySpectrum(
    const LinearComptonBenchmark::SpectrumConfig& config,
    const CliOptions& options) {

    if (options.beamSigmaTransverse_m <= 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires positive transverse beam size.");
    }
    if (options.beamRmsAngle_rad <= 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires positive transverse RMS angle.");
    }
    if (options.beamRelativeEnergySpread < 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires nonnegative relative energy spread.");
    }

    const double p0GeV = std::sqrt(config.electronTotalEnergyGeV * config.electronTotalEnergyGeV
        - Physics::m_e * Physics::m_e);
    const double sigmaPxGeV = p0GeV * options.beamRmsAngle_rad;
    const double beta0 = p0GeV / config.electronTotalEnergyGeV;
    const double sigmaEGeV = options.beamRelativeEnergySpread * config.electronTotalEnergyGeV;
    const double sigmaPzGeV = std::max(1.0e-15,
                                       sigmaEGeV > 0.0 ? sigmaEGeV / beta0
                                                       : 1.0e-12 * p0GeV);
    const double sigmaZ_m = std::max(1.0e-15, options.beamSigmaLongitudinal_m);

    const Vector_t<double, 3> meanR = 0.0;
    Vector_t<double, 3> meanP(0.0);
    meanP[2] = p0GeV;
    Vector_t<double, 3> sigmaR(options.beamSigmaTransverse_m,
                               options.beamSigmaTransverse_m,
                               sigmaZ_m);
    Vector_t<double, 3> sigmaP(sigmaPxGeV,
                               sigmaPxGeV,
                               sigmaPzGeV);
    const Vector_t<double, 3> cutoffR = 4.0;
    const Vector_t<double, 3> cutoffP = 4.0;

    auto pc = makeParticleContainer(sigmaR);
    preallocateParticleCapacity(pc, options.beamParticles);
    ippl::Vector<int, 3> nr = 32;

    MultiVariateGaussian sampler(pc, meanR, meanP, sigmaR, sigmaP, cutoffR, cutoffP, true, true);
    std::size_t numberOfParticles = options.beamParticles;
    sampler.generateParticles(numberOfParticles, nr);

    LinearComptonBenchmark::SpectrumHistogram histogram;
    histogram.binWidthGeV = (config.energyMaxGeV - config.energyMinGeV)
        / static_cast<double>(config.bins);
    histogram.centersGeV.resize(config.bins);
    histogram.densityPerGeV.assign(config.bins, 0.0);
    histogram.counts.assign(config.bins, 0.0);
    for (std::size_t i = 0; i < config.bins; ++i) {
        histogram.centersGeV[i] = config.energyMinGeV
            + (static_cast<double>(i) + 0.5) * histogram.binWidthGeV;
    }

    auto rView = pc->R.getView();
    auto pView = pc->P.getView();
    auto engine = Physics::LinearCompton::makeHostRandomEngine();
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(
        config.wavelength_m);

    histogram.totalWeight = static_cast<double>(pc->getLocalNum());
    for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
        (void)rView;
        Vector_t<double, 3> beamDirection(0.0);
        beamDirection[0] = pView(i)[0];
        beamDirection[1] = pView(i)[1];
        beamDirection[2] = pView(i)[2];
        const double p2 = dot(beamDirection, beamDirection);
        const double momentum = std::sqrt(p2);
        if (momentum <= 0.0) {
            continue;
        }
        beamDirection /= momentum;
        const double electronTotalEnergyGeV = std::sqrt(p2 + Physics::m_e * Physics::m_e);
        const auto kernel = Physics::LinearCompton::makeSamplingKernel(electronTotalEnergyGeV,
                                                                       laserPhotonEnergyGeV,
                                                                       beamDirection,
                                                                       config.laserDirection);
        const auto event = Physics::LinearCompton::sampleEvent(kernel, engine);
        const double energyGeV = event.scatteredPhotonEnergyLabGeV;
        if (energyGeV < config.energyMinGeV || energyGeV >= config.energyMaxGeV) {
            continue;
        }
        const std::size_t bin = static_cast<std::size_t>(
            (energyGeV - config.energyMinGeV) / histogram.binWidthGeV);
        if (bin < histogram.counts.size()) {
            histogram.counts[bin] += 1.0;
        }
    }

    double globalTotalWeight = 0.0;
    MPI_Allreduce(&histogram.totalWeight,
                  &globalTotalWeight,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  ippl::Comm->getCommunicator());
    histogram.totalWeight = globalTotalWeight;

    std::vector<double> globalCounts(histogram.counts.size(), 0.0);
    MPI_Allreduce(histogram.counts.data(),
                  globalCounts.data(),
                  static_cast<int>(histogram.counts.size()),
                  MPI_DOUBLE,
                  MPI_SUM,
                  ippl::Comm->getCommunicator());
    histogram.counts = std::move(globalCounts);

    for (std::size_t i = 0; i < histogram.densityPerGeV.size(); ++i) {
        histogram.densityPerGeV[i] = histogram.counts[i]
            / (histogram.totalWeight * histogram.binWidthGeV);
    }
    return histogram;
}

/**
 * @brief Sample the finite-beam lab-angle spectrum with one Compton event per electron.
 *
 * This uses the same electron sampling and per-particle kernel construction as
 * @ref sampleFiniteBeamEnergySpectrum, but bins the laboratory polar angle of
 * the emitted photon with respect to the benchmark beam axis. The angle is
 * therefore broadened both by Compton kinematics and by the sampled beam
 * divergence.
 */
LinearComptonBenchmark::AngleHistogram sampleFiniteBeamAngularSpectrum(
    const LinearComptonBenchmark::AngleConfig& config,
    const CliOptions& options) {

    if (options.beamSigmaTransverse_m <= 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires positive transverse beam size.");
    }
    if (options.beamRmsAngle_rad <= 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires positive transverse RMS angle.");
    }
    if (options.beamRelativeEnergySpread < 0.0) {
        throw std::runtime_error("Finite-beam benchmark requires nonnegative relative energy spread.");
    }

    const double p0GeV = std::sqrt(config.electronTotalEnergyGeV * config.electronTotalEnergyGeV
        - Physics::m_e * Physics::m_e);
    const double sigmaPxGeV = p0GeV * options.beamRmsAngle_rad;
    const double beta0 = p0GeV / config.electronTotalEnergyGeV;
    const double sigmaEGeV = options.beamRelativeEnergySpread * config.electronTotalEnergyGeV;
    const double sigmaPzGeV = std::max(1.0e-15,
                                       sigmaEGeV > 0.0 ? sigmaEGeV / beta0
                                                       : 1.0e-12 * p0GeV);
    const double sigmaZ_m = std::max(1.0e-15, options.beamSigmaLongitudinal_m);

    const Vector_t<double, 3> meanR = 0.0;
    Vector_t<double, 3> meanP(0.0);
    meanP[2] = p0GeV;
    Vector_t<double, 3> sigmaR(options.beamSigmaTransverse_m,
                               options.beamSigmaTransverse_m,
                               sigmaZ_m);
    Vector_t<double, 3> sigmaP(sigmaPxGeV,
                               sigmaPxGeV,
                               sigmaPzGeV);
    const Vector_t<double, 3> cutoffR = 4.0;
    const Vector_t<double, 3> cutoffP = 4.0;

    auto pc = makeParticleContainer(sigmaR);
    preallocateParticleCapacity(pc, options.beamParticles);
    ippl::Vector<int, 3> nr = 32;

    MultiVariateGaussian sampler(pc, meanR, meanP, sigmaR, sigmaP, cutoffR, cutoffP, true, true);
    std::size_t numberOfParticles = options.beamParticles;
    sampler.generateParticles(numberOfParticles, nr);

    LinearComptonBenchmark::AngleHistogram histogram;
    histogram.binWidthRad = (config.thetaMaxRad - config.thetaMinRad)
        / static_cast<double>(config.bins);
    histogram.centersRad.resize(config.bins);
    histogram.densityPerRad.assign(config.bins, 0.0);
    histogram.counts.assign(config.bins, 0.0);
    for (std::size_t i = 0; i < config.bins; ++i) {
        histogram.centersRad[i] = config.thetaMinRad
            + (static_cast<double>(i) + 0.5) * histogram.binWidthRad;
    }

    auto pView = pc->P.getView();
    auto engine = Physics::LinearCompton::makeHostRandomEngine();
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(
        config.wavelength_m);

    histogram.totalWeight = static_cast<double>(pc->getLocalNum());
    for (std::size_t i = 0; i < pc->getLocalNum(); ++i) {
        Vector_t<double, 3> beamDirection(0.0);
        beamDirection[0] = pView(i)[0];
        beamDirection[1] = pView(i)[1];
        beamDirection[2] = pView(i)[2];
        const double p2 = dot(beamDirection, beamDirection);
        const double momentum = std::sqrt(p2);
        if (momentum <= 0.0) {
            continue;
        }
        beamDirection /= momentum;
        const double electronTotalEnergyGeV = std::sqrt(p2 + Physics::m_e * Physics::m_e);
        const auto kernel = Physics::LinearCompton::makeSamplingKernel(electronTotalEnergyGeV,
                                                                       laserPhotonEnergyGeV,
                                                                       beamDirection,
                                                                       config.laserDirection);
        const auto event = Physics::LinearCompton::sampleEvent(kernel, engine);
        const double thetaRad = LinearComptonBenchmark::photonPolarAngleRad(
            event.scatteredPhotonDirectionLab,
            config.beamDirection);
        if (thetaRad < config.thetaMinRad || thetaRad >= config.thetaMaxRad) {
            continue;
        }
        const std::size_t bin = static_cast<std::size_t>(
            (thetaRad - config.thetaMinRad) / histogram.binWidthRad);
        if (bin < histogram.counts.size()) {
            histogram.counts[bin] += 1.0;
        }
    }

    double globalTotalWeight = 0.0;
    MPI_Allreduce(&histogram.totalWeight,
                  &globalTotalWeight,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  ippl::Comm->getCommunicator());
    histogram.totalWeight = globalTotalWeight;

    std::vector<double> globalCounts(histogram.counts.size(), 0.0);
    MPI_Allreduce(histogram.counts.data(),
                  globalCounts.data(),
                  static_cast<int>(histogram.counts.size()),
                  MPI_DOUBLE,
                  MPI_SUM,
                  ippl::Comm->getCommunicator());
    histogram.counts = std::move(globalCounts);

    for (std::size_t i = 0; i < histogram.densityPerRad.size(); ++i) {
        histogram.densityPerRad[i] = histogram.counts[i]
            / (histogram.totalWeight * histogram.binWidthRad);
    }
    return histogram;
}
}  // namespace

/**
 * @brief Entry point for the standalone linear-Compton benchmark executable.
 *
 * The executable is used in two ways:
 *
 * - directly from the command line to regenerate CSV benchmark references,
 * - indirectly from `TestLinearComptonSpectrum`, where it acts as a small
 *   regression-data generator for the finite-beam cases.
 *
 * Single-electron modes stay purely local. Finite-beam modes initialize IPPL,
 * build a temporary particle container, and reduce histogram counts across MPI
 * ranks before writing the normalized CSV output.
 */
int main(int argc, char** argv) {
    const CliOptions options = parseArguments(argc, argv);

    if (options.finiteBeam) {
        int ipplArgc = argc;
        char** ipplArgv = argv;
        ippl::initialize(ipplArgc, ipplArgv);

        const int previousSeed = Options::seed;
        Options::seed = options.seed;

        if (options.angular) {
            LinearComptonBenchmark::AngleConfig config;
            const auto histogram = sampleFiniteBeamAngularSpectrum(config, options);
            LinearComptonBenchmark::writeAngleCSV(histogram, options.output);
            if (ippl::Comm->rank() == 0) {
                const double p0GeV = std::sqrt(config.electronTotalEnergyGeV * config.electronTotalEnergyGeV
                    - Physics::m_e * Physics::m_e);
                const double sigmaPxGeV = p0GeV * options.beamRmsAngle_rad;
                const double sigmaEGeV = options.beamRelativeEnergySpread
                    * config.electronTotalEnergyGeV;
                const double geomEmittance = options.beamSigmaTransverse_m * options.beamRmsAngle_rad;
                std::cout << "Wrote OPALX finite-beam linear-Compton angular spectrum to " << options.output << '\n'
                          << "Mode = sampled finite beam (MultiVariateGaussian)\n"
                          << "Observable = theta_lab [rad]\n"
                          << "Beam sigma_x,y [m] = " << options.beamSigmaTransverse_m << '\n'
                          << "Beam sigma_z [m] = " << options.beamSigmaLongitudinal_m << " (internally clamped if zero)\n"
                          << "Beam sigma_x',y' [rad] = " << options.beamRmsAngle_rad << '\n'
                          << "Beam sigma_px,py [GeV] = " << sigmaPxGeV << '\n'
                          << "Beam sigma_E / E = " << options.beamRelativeEnergySpread << '\n'
                          << "Beam sigma_E [GeV] = " << sigmaEGeV << '\n'
                          << "Beam geometric emittance [m rad] = " << geomEmittance << '\n'
                          << "Beam macroparticles = " << options.beamParticles << '\n'
                          << "Seed = " << options.seed << '\n'
                          << "Area = " << LinearComptonBenchmark::angleHistogramArea(histogram) << '\n'
                          << "Mean angle [rad] = " << LinearComptonBenchmark::angleHistogramMeanRad(histogram) << '\n';
            }
        } else {
            LinearComptonBenchmark::SpectrumConfig config;
            const auto histogram = sampleFiniteBeamEnergySpectrum(config, options);
            LinearComptonBenchmark::writeSpectrumCSV(histogram, options.output);
            if (ippl::Comm->rank() == 0) {
                const double p0GeV = std::sqrt(config.electronTotalEnergyGeV * config.electronTotalEnergyGeV
                    - Physics::m_e * Physics::m_e);
                const double sigmaPxGeV = p0GeV * options.beamRmsAngle_rad;
                const double sigmaEGeV = options.beamRelativeEnergySpread
                    * config.electronTotalEnergyGeV;
                const double geomEmittance = options.beamSigmaTransverse_m * options.beamRmsAngle_rad;
                std::cout << "Wrote OPALX finite-beam linear-Compton spectrum to " << options.output << '\n'
                          << "Mode = sampled finite beam (MultiVariateGaussian)\n"
                          << "Observable = E_gamma [GeV]\n"
                          << "Beam sigma_x,y [m] = " << options.beamSigmaTransverse_m << '\n'
                          << "Beam sigma_z [m] = " << options.beamSigmaLongitudinal_m << " (internally clamped if zero)\n"
                          << "Beam sigma_x',y' [rad] = " << options.beamRmsAngle_rad << '\n'
                          << "Beam sigma_px,py [GeV] = " << sigmaPxGeV << '\n'
                          << "Beam sigma_E / E = " << options.beamRelativeEnergySpread << '\n'
                          << "Beam sigma_E [GeV] = " << sigmaEGeV << '\n'
                          << "Beam geometric emittance [m rad] = " << geomEmittance << '\n'
                          << "Beam macroparticles = " << options.beamParticles << '\n'
                          << "Seed = " << options.seed << '\n'
                          << "Area = " << LinearComptonBenchmark::histogramArea(histogram) << '\n'
                          << "Mean energy [GeV] = " << LinearComptonBenchmark::histogramMeanEnergyGeV(histogram) << '\n';
            }
        }

        Options::seed = previousSeed;
        ippl::finalize();
        return 0;
    }

    if (options.angular) {
        LinearComptonBenchmark::AngleConfig config;
        LinearComptonBenchmark::AngleHistogram histogram;
        if (options.sampled) {
            const int previousSeed = Options::seed;
            Options::seed = options.seed;
            histogram = LinearComptonBenchmark::sampleLabAngularSpectrum(config, options.samples);
            Options::seed = previousSeed;
        } else {
            histogram = LinearComptonBenchmark::integrateLabAngularSpectrum(config);
        }

        LinearComptonBenchmark::writeAngleCSV(histogram, options.output);

        std::cout << "Wrote OPALX linear-Compton angular spectrum to " << options.output << '\n'
                  << "Mode = " << (options.sampled ? "sampled" : "deterministic") << '\n'
                  << "Observable = theta_lab [rad]\n"
                  << "Area = " << LinearComptonBenchmark::angleHistogramArea(histogram) << '\n'
                  << "Mean angle [rad] = "
                  << LinearComptonBenchmark::angleHistogramMeanRad(histogram) << '\n';
        if (options.sampled) {
            std::cout << "Samples = " << options.samples << '\n'
                      << "Seed = " << options.seed << '\n';
        }
        return 0;
    }

    LinearComptonBenchmark::SpectrumConfig config;
    LinearComptonBenchmark::SpectrumHistogram histogram;
    if (options.sampled) {
        const int previousSeed = Options::seed;
        Options::seed = options.seed;
        histogram = LinearComptonBenchmark::sampleLabSpectrum(config, options.samples);
        Options::seed = previousSeed;
    } else {
        histogram = LinearComptonBenchmark::integrateLabSpectrum(config);
    }

    LinearComptonBenchmark::writeSpectrumCSV(histogram, options.output);

    std::cout << "Wrote OPALX linear-Compton spectrum to " << options.output << '\n'
              << "Mode = " << (options.sampled ? "sampled" : "deterministic") << '\n'
              << "Observable = E_gamma [GeV]\n"
              << "Area = " << LinearComptonBenchmark::histogramArea(histogram) << '\n'
              << "Mean energy [GeV] = "
              << LinearComptonBenchmark::histogramMeanEnergyGeV(histogram) << '\n';
    if (options.sampled) {
        std::cout << "Samples = " << options.samples << '\n'
                  << "Seed = " << options.seed << '\n';
    }
    return 0;
}
