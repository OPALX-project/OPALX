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

//   Class also writes processor load balancing data to file to track parallel
//   calculation efficiency.
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

 
#include "Structure/H5Writer.h"
#include "Structure/SDDSWriter.h"
#include "Structure/StatWriter.h"
#include "Structure/BinConfigWriter.h"
#include "Structure/DirichletPlaneWriter.h"

#include <array>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

class BoundaryGeometry;
class H5PartWrapper;
class BinConfigWriter;
class DirichletPlaneWriter;

class DataSink {
private:
    typedef StatWriter::losses_t losses_t;
    typedef std::unique_ptr<StatWriter> statWriter_t;
    typedef std::unique_ptr<SDDSWriter> sddsWriter_t;
    typedef std::unique_ptr<H5Writer> h5Writer_t;

public:
    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();
    /// One H5 wrapper per particle container (when HDF5 enabled); sizes must match @p numParticleContainers.
    DataSink(const std::vector<H5PartWrapper*>& h5wrappers, bool restart, size_t numParticleContainers);
    DataSink(H5PartWrapper* h5wrapper, bool restart);
    DataSink(H5PartWrapper* h5wrapper);

    /** Basename stem for per-container diagnostics: @c basename if @p numContainers<=1 else @c basename_cN. */
    static std::string diagnosticStemForContainer(
        const std::string& inputBasename, size_t numContainers, size_t index);

    void dumpH5(const std::shared_ptr<PartBunch_t>& beam, Vector_t<double, 3> FDext[]) const;

    void dumpH5(
        const std::shared_ptr<PartBunch_t>& beam,
        const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer) const;

    int dumpH5(
        const std::shared_ptr<PartBunch_t>& beam, Vector_t<double, 3> FDext[], double meanEnergy, double refPr,
        double refPt, double refPz, double refR, double refTheta, double refZ, double azimuth,
        double elevation, bool local) const;

    void dumpSDDS(
        const std::shared_ptr<PartBunch_t>& beam,
        const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
        const double& azimuth) const;

    void dumpSDDS(
        const std::shared_ptr<PartBunch_t>& beam,
        const std::vector<std::array<Vector_t<double, 3>, 2>>& fdextPerContainer,
        const losses_t& losses,
        const double& azimuth) const;

    /** \brief Write cavity information from  H5 file
     */
    void storeCavityInformation();

    void changeH5Wrappers(const std::vector<H5PartWrapper*>& h5wrappers);

    /**
     * Write geometry points and surface triangles to vtk file
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeGeomToVtk(BoundaryGeometry& bg, std::string fn);
    // void writeGeoContourToVtk(BoundaryGeometry &bg, std::string fn);

    /**
     * Write impact number and outgoing secondaries in each time step
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeImpactStatistics(
        const PartBunch_t* beam, long long int& step, size_t& impact, double& sey_num,
        size_t numberOfFieldEmittedParticles, bool nEmissionMode, std::string fn);

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
    void dumpBinConfig(long long step,
                       double time,
                       bool preMerge,
                       const std::vector<std::size_t>& binCounts,
                       const std::vector<double>& binWidths,
                       double xMin,
                       const std::string& fileName);

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
    DirichletPlaneWriter::PlaneDiagnostics dumpDirichletPlane(long long step,
                                                               double time,
                                                               double zPlane,
                                                               const FieldType& field,
                                                               const std::string& solveTag);

private:
    DataSink(const DataSink& ds)         = delete;
    DataSink& operator=(const DataSink&) = delete;

    void rewindLines();

    void init(bool restart, const std::vector<H5PartWrapper*>& h5wrappers, size_t numParticleContainers);

    std::vector<h5Writer_t> h5Writers_m;
    std::vector<statWriter_t> statWriters_m;
    std::vector<sddsWriter_t> sddsWriter_m;

    // Writer for binning configuration JSON output (rank 0 only).
    std::unique_ptr<BinConfigWriter> binConfigWriter_m;

    // Writer for dirichlet-plane diagnostics (rank 0 only).
    std::unique_ptr<DirichletPlaneWriter> dirichletPlaneWriter_m;

    static std::string convertToString(int number, int setw = 5);

    /// needed to create index for vtk file
    unsigned int lossWrCounter_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;
};

inline std::string DataSink::convertToString(int number, int setw) {
    std::stringstream ss;
    ss << std::setw(setw) << std::setfill('0') << number;
    return ss.str();
}

template <typename FieldType>
DirichletPlaneWriter::PlaneDiagnostics DataSink::dumpDirichletPlane(
        long long step,
        double time,
        double zPlane,
        const FieldType& field,
        const std::string& solveTag) {
    if (ippl::Comm->rank() != 0) {
        return DirichletPlaneWriter::PlaneDiagnostics{};
    }

    Inform m("DataSink::dumpDirichletPlane");

    if (!dirichletPlaneWriter_m) {
        m << level4
          << "Creating DirichletPlaneWriter in directory \"data/dirichletplanedumps\"."
          << endl;
        dirichletPlaneWriter_m = std::make_unique<DirichletPlaneWriter>("data/dirichletplanedumps");
    }

    return dirichletPlaneWriter_m->dumpInterpolatedPlane(step, time, zPlane, field, solveTag);
}

#endif  // DataSink_H_
