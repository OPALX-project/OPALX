//
// Class DataSink
//   This class acts as an observer during the calculation. It generates diagnostic
//   output of the accelerated beam such as statistical beam descriptors of particle
//   positions, momenta, beam phase space (emittance) etc. These are written to file
//   at periodic time steps during the calculation.
//
//   This class also writes the full beam phase space to an H5 file at periodic time
//   steps in the calculation (this period is different from that of the statistical
//   numbers).
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef _OPAL_DATA_SINK_H
#define _OPAL_DATA_SINK_H

#include "Structure/BinConfigWriter.h"
#include "Structure/DataSinkOptions.h"
#include "Structure/DirichletPlaneWriter.h"
#include "Structure/H5Writer.h"
#include "Structure/StatWriter.h"

#include <array>
#include <memory>
#include <span>
#include <string>
#include <vector>

class H5PartWrapper;

class DataSink {
private:
    typedef std::unique_ptr<StatWriter> statWriter_t;
    typedef std::unique_ptr<H5Writer> h5Writer_t;

public:
    /**
     * @param h5wrappers One H5 wrapper per particle container if @p options enable HDF5, none
     *        otherwise.
     * @param outputBasename Stem of the statistics and phase-space files.
     * @param inputFileName Input file named in the headers of the statistics files.
     */
    DataSink(
            const DataSinkOptions& options, const std::vector<H5PartWrapper*>& h5wrappers,
            size_t numParticleContainers, const std::string& outputBasename,
            const std::string& inputFileName, bool appendStatistics = false);

    /** Basename stem for per-container diagnostics: @c basename if @p numContainers<=1 else @c
     * basename_cN. */
    static std::string diagnosticStemForContainer(
            const std::string& inputBasename, size_t numContainers, size_t index);

    /**
     * Remove diagnostic output beyond the restored checkpoint before appending a continuation.
     * HDF5 is rewound by time, statistics by per-container path position.
     */
    void rewindToCheckpoint(PartBunch_t& beam);

    /// Moment groups the phase-space (@p psDump) and statistics (@p statDump) outputs read for a
    /// container with or without spin.
    [[nodiscard]] static opalx::moments::Request requiredMoments(
            bool psDump, bool statDump, bool hasSpin);

    /// Moment groups the phase-space (@p psDump) and statistics (@p statDump) outputs of this
    /// sink read for a container with or without spin; no phase-space groups without HDF5.
    [[nodiscard]] opalx::moments::Request momentRequest(
            bool psDump, bool statDump, bool hasSpin) const;

    /**
     * @brief Write H5 phase-space data with per-container external fields.
     * @param beam Borrowed particle bunch; DataSink does not own or retain it.
     * @param fdextPerContainer External field values indexed by particle container.
     * @param moments One snapshot per container holding H5PartWrapperForPT::requiredMoments;
     *        empty containers are skipped.
     */
    void dumpH5(
            PartBunch_t& beam,
            const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
            std::span<const opalx::moments::BeamMoments> moments) const;

    /**
     * @brief Write SDDS statistics with per-container external fields.
     * @param beam Borrowed particle bunch; DataSink does not own or retain it.
     * @param moments One snapshot per container holding StatWriter::requiredMoments; not read
     *        for empty containers.
     */
    void dumpSDDS(
            PartBunch_t& beam,
            const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
            std::span<const opalx::moments::BeamMoments> moments, int reportedBinCount) const;

    /**
     * @brief Append a binning configuration snapshot to a JSON file.
     *
     * @param step       Global tracking step index.
     * @param time       Absolute simulation time (in seconds).
     * @param preMerge   true for pre-merge histogram, false for post-merge.
     * @param binCounts  Per-bin particle counts.
     * @param binWidths  Per-bin bin widths (same length as @p binCounts).
     * @param xMin       Lower bound of the histogram coordinate.
     * @param fileName   Target JSON file name (typically from BinningCmd::getDumpBinsFileName()).
     */
    void dumpBinConfig(
            long long step, double time, bool preMerge, const std::vector<std::size_t>& binCounts,
            const std::vector<double>& binWidths, double xMin, const std::string& fileName);

    /**
     * @brief Interpolate and dump potential values on a 2D z-plane.
     *
     * The interpolation and plane extraction are performed on device by
     * DirichletPlaneWriter; only the 2D plane is copied to host for output.
     *
     * @tparam FieldType IPPL scalar field type.
     * @param step      Global tracking step index.
     * @param time      Absolute simulation time (seconds).
     * @param zPlane    Physical z location of the sampled plane (meters).
     * @param field     Scalar field to sample (`rho`/`phi` depending on solver).
     * @param solveTag  Label for the solver pass (e.g. `legacy`, `binned`).
     *
     * @return Plane diagnostics (mean/variance/sample count).
     */
    template <typename FieldType>
    DirichletPlaneWriter::PlaneDiagnostics dumpDirichletPlane(
            long long step, double time, double zPlane, const FieldType& field,
            const std::string& solveTag);

private:
    DataSink(const DataSink& ds)         = delete;
    DataSink& operator=(const DataSink&) = delete;

    DataSinkOptions options_m;

    std::vector<h5Writer_t> h5Writers_m;
    std::vector<statWriter_t> statWriters_m;

    // Writer for binning configuration JSON output (rank 0 only).
    std::unique_ptr<BinConfigWriter> binConfigWriter_m;

    // Writer for dirichlet-plane diagnostics (rank 0 only).
    std::unique_ptr<DirichletPlaneWriter> dirichletPlaneWriter_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;
};

template <typename FieldType>
DirichletPlaneWriter::PlaneDiagnostics DataSink::dumpDirichletPlane(
        long long step, double time, double zPlane, const FieldType& field,
        const std::string& solveTag) {
    if (ippl::Comm->rank() != 0) {
        return DirichletPlaneWriter::PlaneDiagnostics{};
    }

    Inform m("DataSink::dumpDirichletPlane");

    if (!dirichletPlaneWriter_m) {
        m << level4 << "Creating DirichletPlaneWriter in directory \"data/dirichletplanedumps\"."
          << endl;
        dirichletPlaneWriter_m = std::make_unique<DirichletPlaneWriter>("data/dirichletplanedumps");
    }

    return dirichletPlaneWriter_m->dumpInterpolatedPlane(step, time, zPlane, field, solveTag);
}

#endif  // DataSink_H_
