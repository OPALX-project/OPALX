/**
 * @file FromFile.h
 * @brief Defines the FromFile class used for reading particle phase space from ASCII files.
 */

#ifndef OPALX_FROM_FILE_H
#define OPALX_FROM_FILE_H

#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "Distribution.h"
#include "Ippl.h"
#include "OPALTypes.h"
#include "SamplingBase.hpp"
#include "Utilities/OpalException.h"

/**
 * @class FromFile
 * @brief Implements the sampling method for reading particle phase space from ASCII files.
 *
 * The file must follow this structure (no variations allowed):
 * - Line 1: number of particles (single integer N).
 * - Line 2: column names (space-separated), e.g. "x px y py z pz".
 * - Remaining lines: one data row per particle (6+ columns). Blank lines and
 *   lines starting with '#' are skipped. There must be exactly N data lines.
 * If the file does not follow this structure, an error is thrown.
 *
 * Reading is collective and never holds the whole file on one rank. The constructor reads
 * the two header lines and splits the data bytes into one contiguous range per rank; each
 * rank counts the data lines that start in its range, and an allgather of the counts gives
 * every rank the global row index at each range start (and checks N). generateParticles()
 * then parses only the rows of this rank's slice, starting at the range that holds its first
 * row. The file is therefore read about twice in total, split over the ranks, and host and
 * device memory scale with the local particle count.
 */
class FromFile : public SamplingBase {
public:
    using size_type = ippl::detail::size_type;

    /**
     * @brief Constructor for FromFile.
     * @param pc Shared pointer to ParticleContainer.
     * @param opalDist Borrowed Distribution.
     */
    FromFile(std::shared_ptr<ParticleContainer_t> pc, Distribution_t* opalDist);

    /**
     * @brief Convenience constructor that takes the filename directly.
     *
     * This is primarily intended for unit tests, where constructing a full
     * OPALX Distribution object is unnecessary. It behaves like the main
     * constructor but skips querying the Distribution for the filename.
     */
    FromFile(std::shared_ptr<ParticleContainer_t> pc, const std::string& filename);

    /**
     * @brief Generates particles by reading from file.
     * @param numberOfParticles Number of particles to generate (may be overridden by file).
     * @param nr Number of grid points per direction (not used here).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Time-stepped emission hook for one-shot delayed file-based injection.
     *
     * When a source using FROMFILE has t0 > 0 and a positive NPARTDIST, the sampler
     * will read and inject particles once when [t, t+dt] crosses t0. Initial sampling
     * at t0 == 0 is handled via generateParticles.
     */
    void emitParticles(double t, double dt) override;

private:
    /**
     * @brief Reads the particle count and column header, then counts the data lines of
     *        this rank's byte range and gathers the counts of all ranks (collective).
     * @throws OpalException if the header is invalid or the number of data lines differs
     *         from the declared count (thrown on every rank).
     */
    void readFile(const std::string& filename);

    /**
     * @brief Parses data rows [first, first + count) of the file into a host view (count x 6).
     */
    Kokkos::View<double**, Kokkos::HostSpace> readRows(size_t first, size_t count) const;

    /**
     * @brief Parses the header line to get column indices. Header must contain
     *        all of x, y, z, px, py, pz (or aliases). Throws if any are missing.
     * @param headerLine The header line string.
     * @return Vector of column indices: [x_idx, y_idx, z_idx, px_idx, py_idx, pz_idx]
     */
    std::vector<size_t> parseHeader(const std::string& headerLine);

    /**
     * @brief Normalizes column name for comparison (case-insensitive, whitespace handling).
     * @param name Column name to normalize.
     * @return Normalized column name.
     */
    std::string normalizeColumnName(const std::string& name);

    /// File path to read particle data from
    std::string filename_m;

    /// Number of particles in the file
    size_t numParticles_m;

    /// Byte offset of the first byte after the header line
    std::uint64_t dataOffset_m = 0;

    /// Per rank: number of data lines starting in its byte range, and the offset of the
    /// first line starting there (end of the range if none)
    std::vector<std::uint64_t> rangeLineCounts_m;
    std::vector<std::uint64_t> rangeFirstLine_m;

    /// Column indices mapping: [x_idx, y_idx, z_idx, px_idx, py_idx, pz_idx]
    std::vector<size_t> columnIndices_m;
};

#endif  // OPALX_FROM_FILE_H
