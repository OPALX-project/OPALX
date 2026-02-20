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

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;

/**
 * @class FromFile
 * @brief Implements the sampling method for reading particle phase space from ASCII files.
 *
 * This class reads particle phase space data from a simple ASCII formatted file.
 * The file format consists of:
 * - A header line with column names (e.g., "x y z px py pz")
 * - Data lines with 6 space-separated values per line
 *
 * Alternatively, the first line may be the number of particles, followed by data lines.
 */
class FromFile : public SamplingBase {
public:
    using size_type = ippl::detail::size_type;

    /**
     * @brief Constructor for FromFile.
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param opalDist Shared pointer to Distribution.
     */
    FromFile(std::shared_ptr<ParticleContainer_t> &pc,
             std::shared_ptr<FieldContainer_t> &fc,
             std::shared_ptr<Distribution_t> &opalDist);

    /**
     * @brief Generates particles by reading from file.
     * @param numberOfParticles Number of particles to generate (may be overridden by file).
     * @param nr Number of grid points per direction (not used here).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

private:
    /**
     * @brief Reads and parses the particle data file.
     * @param filename Path to the file to read.
     */
    void readFile(const std::string& filename);

    /**
     * @brief Parses a header line to identify column indices.
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
