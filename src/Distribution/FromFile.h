/**
 * @file FromFile.h
 * @brief Defines the FromFile class used for reading particle phase space from ASCII files.
 */

#ifndef OPALX_FROM_FILE_H
#define OPALX_FROM_FILE_H

#include "Distribution.h"
#include "SamplingBase.hpp"
#include "Ippl.h"
#include "OPALTypes.h"
#include "Utilities/OpalException.h"
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

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
 */
class FromFile : public SamplingBase {
public:
    using size_type = ippl::detail::size_type;

    /**
     * @brief Constructor for FromFile.
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param opalDist Borrowed Distribution.
     */
    FromFile(std::shared_ptr<ParticleContainer_t> pc,
             std::shared_ptr<FieldContainer_t> fc,
             Distribution_t* opalDist);

    /**
     * @brief Convenience constructor that takes the filename directly.
     *
     * This is primarily intended for unit tests, where constructing a full
     * OPALX Distribution object is unnecessary. It behaves like the main
     * constructor but skips querying the Distribution for the filename.
     */
    FromFile(std::shared_ptr<ParticleContainer_t> pc,
             std::shared_ptr<FieldContainer_t> fc,
             const std::string& filename);

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
     * @brief Reads and parses the particle file. Format must be: N (particle count),
     *        then header line (column names), then data lines. Comments ('#') and
     *        blank lines are skipped. Throws if format is invalid.
     * @param filename Path to the file to read.
     */
    void readFile(const std::string& filename);

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

    /// Stored particle data: each element is a 6D phase space vector [x, y, z, px, py, pz]
    std::vector<std::vector<double>> particleData_m;

    /// Number of particles in the file
    size_t numParticles_m;

    /// Column indices mapping: [x_idx, y_idx, z_idx, px_idx, py_idx, pz_idx]
    std::vector<size_t> columnIndices_m;
};

#endif // OPALX_FROM_FILE_H
