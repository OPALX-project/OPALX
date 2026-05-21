//
// Class StatWriter
//   This class writes bunch statistics (*.stat).
//
//   One StatWriter instance corresponds to one output file. DataSink::init creates
//   statWriters_m.size() == numParticleContainers writers: for a single container the
//   file stem is the input basename (e.g. myjob -> myjob.stat); for multiple containers
//   stems are basename + "_c" + index (run_c0.stat, run_c1.stat, ...), see
//   DataSink::diagnosticStemForContainer. Each write() call must use the same
//   particleContainerIndex as that writer's slot so row data comes from
//   beam->getParticleContainer(particleContainerIndex). Shared beam-level quantities
//   (e.g. time t, dt, rmsDensity, nBins) still come from PartBunch_t regardless of index.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "StatWriter.h"

#include <Kokkos_Core.hpp>

#include "AbstractObjects/OpalData.h"
#include "BuildInfo.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

#include <sstream>

StatWriter::StatWriter(const std::string& fname, bool restart) : StatBaseWriter(fname, restart) {}

void StatWriter::fillHeader(const losses_t& losses, const std::string& species) {
    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");
    columns_m.addColumn("s", "double", "m", "Path length");
    columns_m.addColumn("numParticles", "long", "1", "Number of Macro Particles");
    columns_m.addColumn("charge", "double", "1", "Bunch Charge");
    columns_m.addColumn("energy", "double", "MeV", "Mean Bunch Energy");

    columns_m.addColumn("rms_x", "double", "m", "RMS Beamsize in x");
    columns_m.addColumn("rms_y", "double", "m", "RMS Beamsize in y");
    columns_m.addColumn("rms_s", "double", "m", "RMS Beamsize in s");

    columns_m.addColumn("rms_px", "double", "1", "RMS Normalized Momenta in x");
    columns_m.addColumn("rms_py", "double", "1", "RMS Normalized Momenta in y");
    columns_m.addColumn("rms_ps", "double", "1", "RMS Normalized Momenta in s");

    columns_m.addColumn("emit_x", "double", "m", "Normalized Emittance x");
    columns_m.addColumn("emit_y", "double", "m", "Normalized Emittance y");
    columns_m.addColumn("emit_s", "double", "m", "Normalized Emittance s");

    columns_m.addColumn("mean_x", "double", "m", "Mean Beam Position in x");
    columns_m.addColumn("mean_y", "double", "m", "Mean Beam Position in y");
    columns_m.addColumn("mean_s", "double", "m", "Mean Beam Position in s");

    columns_m.addColumn("ref_x", "double", "m", "x coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_y", "double", "m", "y coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_z", "double", "m", "z coordinate of reference particle in lab cs");

    columns_m.addColumn("ref_px", "double", "1", "x momentum of reference particle in lab cs");
    columns_m.addColumn("ref_py", "double", "1", "y momentum of reference particle in lab cs");
    columns_m.addColumn("ref_pz", "double", "1", "z momentum of reference particle in lab cs");

    columns_m.addColumn("max_x", "double", "m", "Max Beamsize in x");
    columns_m.addColumn("max_y", "double", "m", "Max Beamsize in y");
    columns_m.addColumn("max_s", "double", "m", "Max Beamsize in s");

    columns_m.addColumn("xpx", "double", "1", "Correlation xpx");
    columns_m.addColumn("ypy", "double", "1", "Correlation ypy");
    columns_m.addColumn("zpz", "double", "1", "Correlation zpz");

    columns_m.addColumn("Dx", "double", "m", "Dispersion in x");
    columns_m.addColumn("DDx", "double", "1", "Derivative of dispersion in x");
    columns_m.addColumn("Dy", "double", "m", "Dispersion in y");
    columns_m.addColumn("DDy", "double", "1", "Derivative of dispersion in y");

    columns_m.addColumn("Bx_ref", "double", "T", "Bx-Field component of ref particle");
    columns_m.addColumn("By_ref", "double", "T", "By-Field component of ref particle");
    columns_m.addColumn("Bz_ref", "double", "T", "Bz-Field component of ref particle");

    columns_m.addColumn("Ex_ref", "double", "MV/m", "Ex-Field component of ref particle");
    columns_m.addColumn("Ey_ref", "double", "MV/m", "Ey-Field component of ref particle");
    columns_m.addColumn("Ez_ref", "double", "MV/m", "Ez-Field component of ref particle");

    columns_m.addColumn("dE", "double", "MeV", "energy spread of the beam");
    columns_m.addColumn("dt", "double", "ns", "time step size");
    columns_m.addColumn("partsOutside", "double", "1", "outside n*sigma of the beam");

    columns_m.addColumn("DebyeLength", "double", "m", "Debye length in the boosted frame");
    columns_m.addColumn(
            "plasmaParameter", "double", "1",
            "Plasma parameter that gives no. of particles in a Debye sphere");
    columns_m.addColumn("temperature", "double", "K", "Temperature of the beam");
    columns_m.addColumn("rmsDensity", "double", "1", "RMS number density of the beam");

    columns_m.addColumn(
            "nBins", "int", "1",
            "Number of field solver bins, potentially after adaptive binning.");

    /// \todo Options::computePercentiles needs to be brought back
    /*
    if (Options::computePercentiles) {
        columns_m.addColumn("68_Percentile_x", "double", "m",
                            "68.27 percentile (1 sigma of normal distribution) of x-component of
    position"); columns_m.addColumn("68_Percentile_y", "double", "m", "68.27 percentile (1 sigma of
    normal distribution) of y-component of position"); columns_m.addColumn("68_Percentile_z",
    "double", "m", "68.27 percentile (1 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("95_Percentile_x", "double", "m",
                            "95.45 percentile (2 sigma of normal distribution) of x-component of
    position"); columns_m.addColumn("95_Percentile_y", "double", "m", "95.45 percentile (2 sigma of
    normal distribution) of y-component of position"); columns_m.addColumn("95_Percentile_z",
    "double", "m", "95.45 percentile (2 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("99_Percentile_x", "double", "m",
                            "99.73 percentile (3 sigma of normal distribution) of x-component of
    position"); columns_m.addColumn("99_Percentile_y", "double", "m", "99.73 percentile (3 sigma of
    normal distribution) of y-component of position"); columns_m.addColumn("99_Percentile_z",
    "double", "m", "99.73 percentile (3 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("99_99_Percentile_x", "double", "m",
                            "99.994 percentile (4 sigma of normal distribution) of x-component of
    position"); columns_m.addColumn("99_99_Percentile_y", "double", "m", "99.994 percentile (4 sigma
    of normal distribution) of y-component of position"); columns_m.addColumn("99_99_Percentile_z",
    "double", "m", "99.994 percentile (4 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("normalizedEps68Percentile_x", "double", "m",
                            "x-component of normalized emittance at 68 percentile (1 sigma of normal
    distribution)"); columns_m.addColumn("normalizedEps68Percentile_y", "double", "m", "y-component
    of normalized emittance at 68 percentile (1 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps68Percentile_z", "double", "m",
                            "z-component of normalized emittance at 68 percentile (1 sigma of normal
    distribution)");

        columns_m.addColumn("normalizedEps95Percentile_x", "double", "m",
                            "x-component of normalized emittance at 95 percentile (2 sigma of normal
    distribution)"); columns_m.addColumn("normalizedEps95Percentile_y", "double", "m", "y-component
    of normalized emittance at 95 percentile (2 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps95Percentile_z", "double", "m",
                            "z-component of normalized emittance at 95 percentile (2 sigma of normal
    distribution)");

        columns_m.addColumn("normalizedEps99Percentile_x", "double", "m",
                            "x-component of normalized emittance at 99 percentile (3 sigma of normal
    distribution)"); columns_m.addColumn("normalizedEps99Percentile_y", "double", "m", "y-component
    of normalized emittance at 99 percentile (3 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps99Percentile_z", "double", "m",
                            "z-component of normalized emittance at 99 percentile (3 sigma of normal
    distribution)");

        columns_m.addColumn("normalizedEps99_99Percentile_x", "double", "m",
                            "x-component of normalized emittance at 99.99 percentile (4 sigma of
    normal distribution)"); columns_m.addColumn("normalizedEps99_99Percentile_y", "double", "m",
                            "y-component of normalized emittance at 99.99 percentile (4 sigma of
    normal distribution)"); columns_m.addColumn("normalizedEps99_99Percentile_z", "double", "m",
                            "z-component of normalized emittance at 99.99 percentile (4 sigma of
    normal distribution)");
    }
    */
    if (OpalData::getInstance()->isInOPALCyclMode() && ippl::Comm->size() == 1) {
        columns_m.addColumn("R0_x", "double", "m", "R0 Particle position in x");
        columns_m.addColumn("R0_y", "double", "m", "R0 Particle position in y");
        columns_m.addColumn("R0_s", "double", "m", "R0 Particle position in z");

        columns_m.addColumn("P0_x", "double", "1", "R0 Particle momentum in x");
        columns_m.addColumn("P0_y", "double", "1", "R0 Particle momentum in y");
        columns_m.addColumn("P0_s", "double", "1", "R0 Particle momentum in z");
    }

    if (OpalData::getInstance()->isInOPALCyclMode()) {
        columns_m.addColumn("halo_x", "double", "1", "Halo in x");
        columns_m.addColumn("halo_y", "double", "1", "Halo in y");
        columns_m.addColumn("halo_z", "double", "1", "Halo in z");

        columns_m.addColumn("azimuth", "double", "deg", "Azimuth in global coordinates");
    }

    for (size_t i = 0; i < losses.size(); ++i) {
        columns_m.addColumn(losses[i].first, "long", "1", "Number of lost particles in element");
    }

    if (mode_m == std::ios::app) return;

    OPALTimer::Timer simtimer;
    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Statistics data '" << OpalData::getInstance()->getInputFn() << "' " << dateStr << " "
       << timeStr;

    this->addDescription(ss.str(), "stat parameters");

    std::stringstream revision;
    revision << buildinfo::project_name << " " << buildinfo::project_version << " "
             << "git rev. #" << Util::getGitRevision();

    addParameter("processors", "long", "Number of Cores used", ippl::Comm->size());
    addParameter("revision", "string", "git revision of opal", revision.str());
    addParameter("species", "string", "Particle species of container", species);

    this->addInfo("ascii", 1);
}

void StatWriter::write(
        PartBunch_t& beam, Vector_t<double, 3> FDext[], const losses_t& losses,
        const double& azimuth, const size_t npOutside, size_t particleContainerIndex) {
    using ParticleContainer_t               = ParticleContainer<T, Dim>;
    std::shared_ptr<ParticleContainer_t> pc = beam.getParticleContainer(particleContainerIndex);
    if (!pc) {
        return;
    }

    double pathLength         = pc->get_sPos();
    const std::string species = beam.getParticleName(particleContainerIndex);

    // First write to this writer's .stat file emits SDDS header via fillHeader/writeHeader.
    // File vs. container: this object was constructed with the stem for particleContainerIndex;
    // pc must be beam->getParticleContainer(particleContainerIndex) (caller responsibility).

    const size_t numParticles = pc->getTotalNum();
    double Q                  = numParticles == 0 ? 0.0 : pc->getTotalCharge();

    if (ippl::Comm->rank() != 0) {
        return;
    }

    fillHeader(losses, species);

    this->open();

    this->writeHeader();

    Vector_t<double, 3> rmsR(0.0);
    Vector_t<double, 3> rmsP(0.0);
    Vector_t<double, 3> normEmit(0.0);
    Vector_t<double, 3> meanR(0.0);
    Vector_t<double, 3> maxR(0.0);
    Vector_t<double, 3> rmsRP(0.0);
    double energy          = 0.0;
    double dE              = 0.0;
    double dx              = 0.0;
    double ddx             = 0.0;
    double dy              = 0.0;
    double ddy             = 0.0;
    double debyeLength     = 0.0;
    double plasmaParameter = 0.0;
    double temperature     = 0.0;

    if (numParticles > 0) {
        energy          = pc->getMeanKineticEnergy();
        rmsR            = pc->getRmsR();
        rmsP            = pc->getRmsP();
        normEmit        = pc->getNormEmit();
        meanR           = pc->getMeanR();
        maxR            = pc->getMaxR();
        rmsRP           = pc->getRmsRP();
        dE              = pc->getStdKineticEnergy();
        dx              = pc->getDx();
        ddx             = pc->getDDx();
        dy              = pc->getDy();
        ddy             = pc->getDDy();
        debyeLength     = pc->getDebyeLength();
        plasmaParameter = pc->getPlasmaParameter();
        temperature     = pc->getTemperature();
    }

    columns_m.addColumnValue("t", beam.getT() * Units::s2ns);        // 1
    columns_m.addColumnValue("s", pathLength);                       // 2
    columns_m.addColumnValue("numParticles", numParticles);          // 3
    columns_m.addColumnValue("charge", Q);                           // 4
    columns_m.addColumnValue("energy", energy);                      // 5

    columns_m.addColumnValue("rms_x", rmsR(0));  // 6
    columns_m.addColumnValue("rms_y", rmsR(1));  // 7
    columns_m.addColumnValue("rms_s", rmsR(2));  // 8

    columns_m.addColumnValue("rms_px", rmsP(0));  // 9
    columns_m.addColumnValue("rms_py", rmsP(1));  // 10
    columns_m.addColumnValue("rms_ps", rmsP(2));  // 11

    columns_m.addColumnValue("emit_x", normEmit(0));  // 12
    columns_m.addColumnValue("emit_y", normEmit(1));  // 13
    columns_m.addColumnValue("emit_s", normEmit(2));  // 14

    columns_m.addColumnValue("mean_x", meanR(0));  // 15
    columns_m.addColumnValue("mean_y", meanR(1));  // 16
    columns_m.addColumnValue("mean_s", meanR(2));  // 17

    columns_m.addColumnValue("ref_x", pc->getRefPartR()(0));  // 18
    columns_m.addColumnValue("ref_y", pc->getRefPartR()(1));  // 19
    columns_m.addColumnValue("ref_z", pc->getRefPartR()(2));  // 20

    columns_m.addColumnValue("ref_px", pc->getRefPartP()(0));  // 21
    columns_m.addColumnValue("ref_py", pc->getRefPartP()(1));  // 22
    columns_m.addColumnValue("ref_pz", pc->getRefPartP()(2));  // 23

    columns_m.addColumnValue("max_x", maxR(0));  // 24
    columns_m.addColumnValue("max_y", maxR(1));  // 25
    columns_m.addColumnValue("max_s", maxR(2));  // 26

    // Write out Courant Snyder parameters.
    columns_m.addColumnValue("xpx", rmsRP(0));  // 27
    columns_m.addColumnValue("ypy", rmsRP(1));  // 28
    columns_m.addColumnValue("zpz", rmsRP(2));  // 29

    // Write out dispersion.
    columns_m.addColumnValue("Dx", dx);    // 30
    columns_m.addColumnValue("DDx", ddx);  // 31
    columns_m.addColumnValue("Dy", dy);    // 32
    columns_m.addColumnValue("DDy", ddy);  // 33

    // Write head/reference particle/tail field information.
    columns_m.addColumnValue("Bx_ref", FDext[0](0));  // 34 B-ref x
    columns_m.addColumnValue("By_ref", FDext[0](1));  // 35 B-ref y
    columns_m.addColumnValue("Bz_ref", FDext[0](2));  // 36 B-ref z

    columns_m.addColumnValue("Ex_ref", FDext[1](0));  // 37 E-ref x
    columns_m.addColumnValue("Ey_ref", FDext[1](1));  // 38 E-ref y
    columns_m.addColumnValue("Ez_ref", FDext[1](2));  // 39 E-ref z

    columns_m.addColumnValue("dE", dE);                          // 40 dE energy spread
    columns_m.addColumnValue("dt", beam.getdT() * Units::s2ns);  // 41 dt time step size
    columns_m.addColumnValue("partsOutside", npOutside);  // 42 number of particles outside n*sigma

    columns_m.addColumnValue("DebyeLength", debyeLength);           // 43 Debye length
    columns_m.addColumnValue("plasmaParameter", plasmaParameter);  // 43 plasma parameter
    columns_m.addColumnValue("temperature", temperature);          // 44 Temperature
    columns_m.addColumnValue("rmsDensity", beam.get_rmsDensity());          // 45 RMS number density
    columns_m.addColumnValue("nBins", beam.getCurrentNBins());

    /*
    if (Options::computePercentiles) {
        columns_m.addColumnValue("68_Percentile_x", beam->get_68Percentile()[0]);
        columns_m.addColumnValue("68_Percentile_y", beam->get_68Percentile()[1]);
        columns_m.addColumnValue("68_Percentile_z", beam->get_68Percentile()[2]);

        columns_m.addColumnValue("95_Percentile_x", beam->get_95Percentile()[0]);
        columns_m.addColumnValue("95_Percentile_y", beam->get_95Percentile()[1]);
        columns_m.addColumnValue("95_Percentile_z", beam->get_95Percentile()[2]);

        columns_m.addColumnValue("99_Percentile_x", beam->get_99Percentile()[0]);
        columns_m.addColumnValue("99_Percentile_y", beam->get_99Percentile()[1]);
        columns_m.addColumnValue("99_Percentile_z", beam->get_99Percentile()[2]);

        columns_m.addColumnValue("99_99_Percentile_x", beam->get_99_99Percentile()[0]);
        columns_m.addColumnValue("99_99_Percentile_y", beam->get_99_99Percentile()[1]);
        columns_m.addColumnValue("99_99_Percentile_z", beam->get_99_99Percentile()[2]);

        columns_m.addColumnValue(
            "normalizedEps68Percentile_x", beam->get_normalizedEps_68Percentile()[0]);
        columns_m.addColumnValue(
            "normalizedEps68Percentile_y", beam->get_normalizedEps_68Percentile()[1]);
        columns_m.addColumnValue(
            "normalizedEps68Percentile_z", beam->get_normalizedEps_68Percentile()[2]);

        columns_m.addColumnValue(
            "normalizedEps95Percentile_x", beam->get_normalizedEps_95Percentile()[0]);
        columns_m.addColumnValue(
            "normalizedEps95Percentile_y", beam->get_normalizedEps_95Percentile()[1]);
        columns_m.addColumnValue(
            "normalizedEps95Percentile_z", beam->get_normalizedEps_95Percentile()[2]);

        columns_m.addColumnValue(
            "normalizedEps99Percentile_x", beam->get_normalizedEps_99Percentile()[0]);
        columns_m.addColumnValue(
            "normalizedEps99Percentile_y", beam->get_normalizedEps_99Percentile()[1]);
        columns_m.addColumnValue(
            "normalizedEps99Percentile_z", beam->get_normalizedEps_99Percentile()[2]);

        columns_m.addColumnValue(
            "normalizedEps99_99Percentile_x", beam->get_normalizedEps_99_99Percentile()[0]);
        columns_m.addColumnValue(
            "normalizedEps99_99Percentile_y", beam->get_normalizedEps_99_99Percentile()[1]);
        columns_m.addColumnValue(
            "normalizedEps99_99Percentile_z", beam->get_normalizedEps_99_99Percentile()[2]);
    }
    */
    if (OpalData::getInstance()->isInOPALCyclMode()) {
        if (ippl::Comm->size() == 1) {
            if (pc->getLocalNum() > 0) {
                auto rDev  = pc->R.getView();
                auto pDev  = pc->P.getView();
                auto rHost = Kokkos::create_mirror_view(rDev);
                auto pHost = Kokkos::create_mirror_view(pDev);
                Kokkos::deep_copy(rHost, rDev);
                Kokkos::deep_copy(pHost, pDev);
                columns_m.addColumnValue("R0_x", rHost(0)[0]);
                columns_m.addColumnValue("R0_y", rHost(0)[1]);
                columns_m.addColumnValue("R0_s", rHost(0)[2]);
                columns_m.addColumnValue("P0_x", pHost(0)[0]);
                columns_m.addColumnValue("P0_y", pHost(0)[1]);
                columns_m.addColumnValue("P0_s", pHost(0)[2]);
            } else {
                columns_m.addColumnValue("R0_x", 0.0);
                columns_m.addColumnValue("R0_y", 0.0);
                columns_m.addColumnValue("R0_s", 0.0);
                columns_m.addColumnValue("P0_x", 0.0);
                columns_m.addColumnValue("P0_y", 0.0);
                columns_m.addColumnValue("P0_s", 0.0);
            }
        }
        Vector_t<double, 3> halo = beam.get_halo();
        columns_m.addColumnValue("halo_x", halo(0));
        columns_m.addColumnValue("halo_y", halo(1));
        columns_m.addColumnValue("halo_z", halo(2));

        columns_m.addColumnValue("azimuth", azimuth);
    }

    for (size_t i = 0; i < losses.size(); ++i) {
        long unsigned int loss = losses[i].second;
        columns_m.addColumnValue(losses[i].first, loss);
    }

    this->writeRow();

    this->close();
}
