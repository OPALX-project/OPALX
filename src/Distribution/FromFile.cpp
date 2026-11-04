#include "FromFile.h"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string_view>
#include "AbstractObjects/OpalData.h"
#include "Distribution.h"
#include "SamplingBase.hpp"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

FromFile::FromFile(std::shared_ptr<ParticleContainer_t> pc, Distribution_t* opalDist)
    : SamplingBase(pc, opalDist), numParticles_m(0) {
    // Get filename from distribution
    filename_m = opalDist->getFilename();

    if (filename_m.empty()) {
        throw OpalException(
                "FromFile::FromFile",
                "FNAME attribute must be set for FROMFILE distribution type.");
    }

    // Resolve file path (check relative to input file directory)
    namespace fs = std::filesystem;
    if (!fs::exists(filename_m)) {
        // Try relative to input file directory
        std::string inputDir = OpalData::getInstance()->getInputFn();
        if (!inputDir.empty()) {
            fs::path inputPath(inputDir);
            fs::path filePath = inputPath.parent_path() / filename_m;
            if (fs::exists(filePath)) {
                filename_m = filePath.string();
            }
        }
    }

    // Read and parse the file
    readFile(filename_m);
}

FromFile::FromFile(std::shared_ptr<ParticleContainer_t> pc, const std::string& filename)
    : SamplingBase(pc), numParticles_m(0) {
    filename_m = filename;

    if (filename_m.empty()) {
        throw OpalException("FromFile::FromFile", "Filename must not be empty.");
    }

    // Read and parse the file
    readFile(filename_m);
}

namespace {
    /// Sequential reader of newline-terminated lines from a byte offset, in large chunks.
    class LineReader {
    public:
        LineReader(const std::string& filename, std::uint64_t offset)
            : file_m(std::fopen(filename.c_str(), "rb")), buffer_m(1 << 24), bufferStart_m(offset) {
            if (file_m == nullptr || fseeko(file_m, static_cast<off_t>(offset), SEEK_SET) != 0) {
                if (file_m != nullptr) {
                    std::fclose(file_m);
                }
                throw OpalException(
                        "FromFile::LineReader", "Couldn't open file '" + filename + "'.");
            }
        }

        ~LineReader() { std::fclose(file_m); }

        LineReader(const LineReader&)            = delete;
        LineReader& operator=(const LineReader&) = delete;

        /// The next line without its newline and the byte offset of its first character.
        bool next(std::string_view& line, std::uint64_t& start) {
            for (;;) {
                const char* first = buffer_m.data() + begin_m;
                const void* nl    = std::memchr(first, '\n', end_m - begin_m);
                if (nl != nullptr || (eof_m && begin_m < end_m)) {
                    const size_t length =
                            nl != nullptr ? static_cast<const char*>(nl) - first : end_m - begin_m;
                    line  = std::string_view(first, length);
                    start = bufferStart_m + begin_m;
                    begin_m += nl != nullptr ? length + 1 : length;
                    return true;
                }
                if (eof_m) {
                    return false;
                }
                refill();
            }
        }

    private:
        void refill() {
            std::memmove(buffer_m.data(), buffer_m.data() + begin_m, end_m - begin_m);
            bufferStart_m += begin_m;
            end_m -= begin_m;
            begin_m = 0;
            if (end_m == buffer_m.size()) {
                buffer_m.resize(2 * buffer_m.size());  // a line longer than the buffer
            }
            const size_t got = std::fread(buffer_m.data() + end_m, 1, buffer_m.size() - end_m, file_m);
            end_m += got;
            eof_m = got == 0;
        }

        std::FILE* file_m;
        std::vector<char> buffer_m;
        std::uint64_t bufferStart_m;
        size_t begin_m = 0;
        size_t end_m   = 0;
        bool eof_m     = false;
    };

    /// Data lines are the lines that are neither blank nor comments ('#').
    bool isDataLine(std::string_view line) {
        const size_t first = line.find_first_not_of(" \t\r");
        return first != std::string_view::npos && line[first] != '#';
    }

    std::string trimmed(std::string_view line) {
        const size_t first = line.find_first_not_of(" \t\r");
        if (first == std::string_view::npos) {
            return {};
        }
        const size_t last = line.find_last_not_of(" \t\r");
        return std::string(line.substr(first, last - first + 1));
    }
}  // namespace

void FromFile::readFile(const std::string& filename) {
    namespace fs = std::filesystem;
    std::error_code ec;
    const std::uint64_t fileSize = fs::file_size(filename, ec);
    if (ec) {
        throw OpalException("FromFile::readFile", "Couldn't open file '" + filename + "'.");
    }

    // Required format: line 1 = particle count, line 2 = column names, then data lines
    std::string line1, headerLine;
    {
        LineReader reader(filename, 0);
        std::string_view line;
        std::uint64_t start = 0;
        auto nextDataLine   = [&]() {
            while (reader.next(line, start)) {
                if (isDataLine(line)) {
                    return true;
                }
            }
            return false;
        };
        if (!nextDataLine()) {
            throw OpalException("FromFile::readFile",
                                "File '" + filename + "' is empty or has no valid lines. "
                                "Expected format: first line = number of particles, second line = column names.");
        }
        line1 = trimmed(line);
        if (!nextDataLine()) {
            throw OpalException("FromFile::readFile",
                                "Missing header line in '" + filename + "'. "
                                "Expected format: N\\n column_names\\n data... (e.g. x px y py z pz).");
        }
        headerLine   = trimmed(line);
        dataOffset_m = std::min<std::uint64_t>(start + line.size() + 1, fileSize);
    }

    // First line must be a single integer (particle count)
    std::istringstream iss1(line1);
    long long npart = 0;
    if (!(iss1 >> npart) || npart < 0 || !iss1.eof()) {
        throw OpalException("FromFile::readFile",
                            "First non-empty line in '" + filename + "' must be the number of particles (single integer). "
                            "Got: '" + line1 + "'. Expected format: N\\n column_names\\n data...");
    }
    numParticles_m  = static_cast<size_t>(npart);
    columnIndices_m = parseHeader(headerLine);

    // Count the data lines starting in this rank's byte range. A line belongs to the range
    // holding its first byte, so the line straddling the range start belongs to the left rank.
    const int rank            = ippl::Comm->rank();
    const int nranks          = std::max(1, ippl::Comm->size());
    const std::uint64_t bytes = fileSize - dataOffset_m;
    auto rangeBegin           = [&](int r) {
        return dataOffset_m + static_cast<std::uint64_t>(
                       static_cast<long double>(bytes) * r / nranks);
    };
    const std::uint64_t begin = rangeBegin(rank);
    const std::uint64_t end   = rangeBegin(rank + 1);
    std::uint64_t local[2]    = {0, end};
    if (begin < end) {
        LineReader reader(filename, begin == dataOffset_m ? begin : begin - 1);
        std::string_view line;
        std::uint64_t start = 0;
        if (begin != dataOffset_m) {
            reader.next(line, start);  // remainder of the line holding byte begin - 1
        }
        bool first = true;
        while (reader.next(line, start) && start < end) {
            if (first) {
                local[1] = start;
                first    = false;
            }
            local[0] += isDataLine(line) ? 1 : 0;
        }
    }
    std::vector<std::uint64_t> gathered(2 * nranks);
    MPI_Allgather(local, 2, MPI_UINT64_T, gathered.data(), 2, MPI_UINT64_T,
                  ippl::Comm->getCommunicator());
    rangeLineCounts_m.resize(nranks);
    rangeFirstLine_m.resize(nranks);
    std::uint64_t dataLines = 0;
    for (int r = 0; r < nranks; ++r) {
        rangeLineCounts_m[r] = gathered[2 * r];
        rangeFirstLine_m[r]  = gathered[2 * r + 1];
        dataLines += rangeLineCounts_m[r];
    }

    if (dataLines == 0) {
        throw OpalException(
                "FromFile::readFile", "No valid particle data found in '" + filename + "'.");
    }
    if (dataLines != numParticles_m) {
        throw OpalException(
                "FromFile::readFile", "Number of data lines (" + std::to_string(dataLines)
                                              + ") does not match declared count ("
                                              + std::to_string(numParticles_m) + ") in '"
                                              + filename + "'.");
    }
}

Kokkos::View<double**, Kokkos::HostSpace> FromFile::readRows(size_t first, size_t count) const {
    Kokkos::View<double**, Kokkos::HostSpace> rows("FromFile::rows", count, 6);
    if (count == 0) {
        return rows;
    }

    // The range holding data row `first`, and the number of data lines to skip inside it
    size_t range        = 0;
    std::uint64_t index = 0;
    while (range < rangeLineCounts_m.size() && index + rangeLineCounts_m[range] <= first) {
        index += rangeLineCounts_m[range++];
    }
    if (range == rangeLineCounts_m.size()) {
        throw OpalException("FromFile::readRows", "Data row " + std::to_string(first)
                                                          + " is beyond the end of '" + filename_m
                                                          + "'.");
    }

    const size_t maxColIdx = *std::max_element(columnIndices_m.begin(), columnIndices_m.end());
    std::vector<double> values(maxColIdx + 1);
    LineReader reader(filename_m, rangeFirstLine_m[range]);
    std::string_view line;
    std::uint64_t start = 0;
    size_t row          = 0;
    while (row < count && reader.next(line, start)) {
        if (!isDataLine(line) || index++ < first) {
            continue;
        }
        // strtod skips newlines, so parse a terminated copy that cannot run into the next row
        const std::string text(line);
        const char* cursor = text.c_str();
        for (size_t col = 0; col <= maxColIdx; ++col) {
            char* next = nullptr;
            values[col] = std::strtod(cursor, &next);
            if (next == cursor) {
                throw OpalException(
                        "FromFile::readRows",
                        "Data row " + std::to_string(first + row) + " in '" + filename_m
                                + "' has fewer columns (" + std::to_string(col)
                                + ") than required (index " + std::to_string(maxColIdx) + ").");
            }
            cursor = next;
        }
        for (int j = 0; j < 6; ++j) {
            rows(row, j) = values[columnIndices_m[j]];
        }
        ++row;
    }
    if (row != count) {
        throw OpalException("FromFile::readRows", "'" + filename_m + "' ended after "
                                                          + std::to_string(first + row)
                                                          + " data rows.");
    }
    return rows;
}

std::vector<size_t> FromFile::parseHeader(const std::string& headerLine) {
    std::istringstream headerStream(headerLine);
    std::vector<std::string> columnNames;
    std::string token;

    while (headerStream >> token) {
        columnNames.push_back(token);
    }

    // Map column names to indices
    std::vector<size_t> indices(6, SIZE_MAX);

    for (size_t i = 0; i < columnNames.size(); ++i) {
        std::string normalized = normalizeColumnName(columnNames[i]);

        if (normalized == "x" && indices[0] == SIZE_MAX) {
            indices[0] = i;
        } else if (normalized == "y" && indices[1] == SIZE_MAX) {
            indices[1] = i;
        } else if (normalized == "z" && indices[2] == SIZE_MAX) {
            indices[2] = i;
        } else if ((normalized == "px" || normalized == "px/momentumx") && indices[3] == SIZE_MAX) {
            indices[3] = i;
        } else if ((normalized == "py" || normalized == "py/momentumy") && indices[4] == SIZE_MAX) {
            indices[4] = i;
        } else if ((normalized == "pz" || normalized == "pz/momentumz") && indices[5] == SIZE_MAX) {
            indices[5] = i;
        }
    }

    // Require all six phase-space columns in the header
    for (size_t i = 0; i < 6; ++i) {
        if (indices[i] == SIZE_MAX) {
            const char* names[] = {"x", "y", "z", "px", "py", "pz"};
            throw OpalException(
                    "FromFile::parseHeader",
                    "Header must contain all column names (x, y, z, px, py, pz). "
                    "Missing or unrecognized: '"
                            + std::string(names[i]) + "'.");
        }
    }

    return indices;
}

std::string FromFile::normalizeColumnName(const std::string& name) {
    std::string normalized = name;

    // Convert to lowercase
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c) {
        return std::tolower(c);
    });

    // Remove whitespace
    normalized.erase(
            std::remove_if(
                    normalized.begin(), normalized.end(),
                    [](unsigned char c) {
                        return std::isspace(c);
                    }),
            normalized.end());

    return normalized;
}

void FromFile::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
    if (emissionModel_m != "NONE")
        throw OpalException(
                "FromFile::generateParticles",
                "EMISSIONMODEL '" + emissionModel_m
                        + "' is not supported for FROMFILE distributions");

    // Only generate during initial sampling (t0 <= 0). For t0 > 0, this
    // distribution is time-independent and should not contribute here unless
    // explicitly triggered via emitParticles (which sets hasEmittedOnce_m).
    if (t0_m > 0.0 && !hasEmittedOnce_m) {
        return;
    }
    Inform m("FromFile::generateParticles");

    // Use number of particles from file if available, otherwise use requested number
    size_t totalParticles = (numParticles_m > 0) ? numParticles_m : numberOfParticles;

    // If file has fewer particles than requested, use file count
    if (numParticles_m > 0 && numParticles_m < numberOfParticles) {
        m << "Warning: File contains " << numParticles_m << " particles, but " << numberOfParticles
          << " were requested." << endl;
        m << "Using " << numParticles_m << " particles from file." << endl;
        totalParticles = numParticles_m;
    }

    // If file has more particles than requested, use requested count
    if (numParticles_m > numberOfParticles) {
        totalParticles = numberOfParticles;
    }

    numberOfParticles = totalParticles;

    // Distribute particles across MPI ranks (capacity-aware)
    const int rank        = ippl::Comm->rank();
    const int nranks      = std::max(1, ippl::Comm->size());
    const size_t nranks_u = static_cast<size_t>(nranks);
    size_t nlocal         = pc_m ? computeLocalEmitCount(totalParticles)
                                 : (totalParticles / nranks_u
                            + (static_cast<size_t>(rank) < (totalParticles % nranks_u) ? 1 : 0));
    // Use Mpi scan to get the start index for each rank (since they are all potentially different!)
    size_t startIdx = 0;
    if (pc_m && ippl::Comm->size() > 0) {
        unsigned long nlocalUL   = static_cast<unsigned long>(nlocal);
        unsigned long startIdxUL = 0;
        MPI_Exscan(
                &nlocalUL, &startIdxUL, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                ippl::Comm->getCommunicator());
        if (rank > 0) {
            startIdx = static_cast<size_t>(startIdxUL);
        }
    }

    // Allocate particles, appending after any existing ones.
    const size_t nlocalCurrent = pc_m->getLocalNum();
    pc_m->createParticles(nlocal);

    // Use subviews over the newly created slots, matching the Gaussian pattern.
    view_type RviewFull = pc_m->R.getView();
    auto Rview = Kokkos::subview(RviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));
    view_type PviewFull = pc_m->P.getView();
    auto Pview = Kokkos::subview(PviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));

    // Parse only this rank's rows [startIdx, startIdx + nlocal)
    auto hostParticleData = readRows(startIdx, nlocal);

    // Create device mirror with layout compatible for Host->Device deep_copy
    // (required when no common execution space, e.g. Host vs Cuda)
    auto deviceParticleData =
            Kokkos::create_mirror(Kokkos::DefaultExecutionSpace::memory_space(), hostParticleData);
    Kokkos::deep_copy(deviceParticleData, hostParticleData);

    // Copy particle data into the subview [0, nlocal).
    Kokkos::parallel_for(
            "FromFile::generateParticles", nlocal, KOKKOS_LAMBDA(const size_t k) {
                Rview(k)[0] = deviceParticleData(k, 0);  // x
                Rview(k)[1] = deviceParticleData(k, 1);  // y
                Rview(k)[2] = deviceParticleData(k, 2);  // z
                Pview(k)[0] = deviceParticleData(k, 3);  // px
                Pview(k)[1] = deviceParticleData(k, 4);  // py
                Pview(k)[2] = deviceParticleData(k, 5);  // pz
            });
    Kokkos::fence();

    // Apply per-emission-source offsets after all mean-fixing/corrections.
    // EMISSIONSOURCE offsets are expected to translate the generated bunch
    // without being affected by the internal "fix mean" logic.
    const Vector_t<double, 3> R0 = R0_m;
    const Vector_t<double, 3> P0 = P0_m;
    Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(const size_t k) {
                Rview(k) += R0;
                Pview(k) += P0;
            });

    Inform mALL("FromFile::generateParticles", INFORM_ALL_NODES);
    mALL << "FromFile: Loaded " << totalParticles << " particles from file '" << filename_m << "'"
         << endl;
    ippl::Comm->barrier();
    mALL << "Rank " << rank << ": " << nlocal << " local particles" << endl;
    ippl::Comm->barrier();

    fillPolarization(nlocalCurrent, nlocal);

    pc_m->markMomentsDirty();
}

void FromFile::emitParticles(double t, double dt) {
    // One-shot delayed emission for FROMFILE sources.
    const double tStart = t;
    const double tEnd   = t + dt;

    if (hasEmittedOnce_m) {
        // Don't sample again.
        return;
    }

    if (t0_m <= 0.0) {
        return;
    }

    if (!(tStart <= t0_m && t0_m < tEnd)) {
        // Not time to emit yet.
        return;
    }

    // If bound to a Distribution, respect its NPARTDIST; otherwise fall back to file count.
    size_t requested = numParticles_m;
    if (opalDist_m && opalDist_m->getNumParticles() > 0) {
        requested = opalDist_m->getNumParticles();
    }
    if (requested == 0) {
        // Short circuit: no particles to emit.
        return;
    }

    const size_t nlocalBefore = pc_m->getLocalNum();

    hasEmittedOnce_m = true;
    Vector_t<double, 3> dummyNr(0.0);
    generateParticles(requested, dummyNr);

    // Set per-particle dt for newly created particles.
    // switchToUnitlessPositions(use_dt_per_particle=true) in pushParticles scales
    // R by 1/(c * dtview(i)).  New particles come with dtview = 0 (Kokkos
    // zero-init), so that division produces inf, and the subsequent multiply by
    // c*0 in switchOffUnitlessPositions gives inf*0 = NaN positions.
    const size_t nlocalAfter = pc_m->getLocalNum();
    const size_t nNew        = nlocalAfter - nlocalBefore;
    if (nNew > 0) {
        const double fracDt = tEnd - t0_m;
        auto dtview         = pc_m->dt.getView();
        const size_t offset = nlocalBefore;
        Kokkos::parallel_for(
                "FromFile_setDt", nNew,
                KOKKOS_LAMBDA(const size_t j) { dtview(offset + j) = fracDt; });
        Kokkos::fence();
    }
}
