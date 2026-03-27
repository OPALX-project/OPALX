#include "Distribution.h"
#include "SamplingBase.hpp"
#include "FromFile.h"
#include "Utilities/OpalException.h"
#include "AbstractObjects/OpalData.h"
#include "Utility/Inform.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <filesystem>

FromFile::FromFile(std::shared_ptr<ParticleContainer_t> pc,
                   std::shared_ptr<FieldContainer_t> fc,
                   std::shared_ptr<Distribution_t> opalDist)
    : SamplingBase(pc, fc, opalDist), numParticles_m(0) {

    // Get filename from distribution
    filename_m = opalDist->getFilename();
    
    if (filename_m.empty()) {
        throw OpalException("FromFile::FromFile",
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

FromFile::FromFile(std::shared_ptr<ParticleContainer_t> pc,
                   std::shared_ptr<FieldContainer_t> fc,
                   const std::string& filename)
    : SamplingBase(pc, fc), numParticles_m(0) {

    filename_m = filename;

    if (filename_m.empty()) {
        throw OpalException("FromFile::FromFile",
                            "Filename must not be empty.");
    }

    // Read and parse the file
    readFile(filename_m);
}

void FromFile::readFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw OpalException("FromFile::readFile",
                            "Couldn't open file '" + filename + "'.");
    }

    // Helper to read next non-empty, non-comment line; returns false if none left
    auto nextLine = [&file](std::string& out) -> bool {
        std::string line;
        while (std::getline(file, line)) {
            auto first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos)
                continue;
            auto last = line.find_last_not_of(" \t\r\n");
            line = line.substr(first, last - first + 1);
            if (line.empty() || line[0] == '#')
                continue;
            out = line;
            return true;
        }
        return false;
    };

    // Required format: line 1 = particle count, line 2 = column names, then data lines
    std::string line1;
    if (!nextLine(line1)) {
        throw OpalException("FromFile::readFile",
                            "File '" + filename + "' is empty or has no valid lines. "
                            "Expected format: first line = number of particles, second line = column names.");
    }

    // First line must be a single integer (particle count)
    std::istringstream iss1(line1);
    long long npart = 0;
    if (!(iss1 >> npart) || npart < 0 || !iss1.eof()) {
        throw OpalException("FromFile::readFile",
                            "First non-empty line in '" + filename + "' must be the number of particles (single integer). "
                            "Got: '" + line1 + "'. Expected format: N\\n column_names\\n data...");
    }
    size_t expectedNumParticles = static_cast<size_t>(npart);

    std::string headerLine;
    if (!nextLine(headerLine)) {
        throw OpalException("FromFile::readFile",
                            "Missing header line in '" + filename + "'. "
                            "Expected format: N\\n column_names\\n data... (e.g. x px y py z pz).");
    }

    columnIndices_m = parseHeader(headerLine);

    // Read data lines (comments and blank lines are skipped)
    particleData_m.clear();
    size_t lineNumber = 3;

    std::string line;
    while (std::getline(file, line)) {
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            ++lineNumber;
            continue;
        }
        auto last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);
        if (line.empty() || line[0] == '#') {
            ++lineNumber;
            continue;
        }

        std::istringstream lineStream(line);
        std::vector<double> values;
        double value;
        while (lineStream >> value)
            values.push_back(value);

        size_t maxColIdx = *std::max_element(columnIndices_m.begin(), columnIndices_m.end());
        if (values.size() <= maxColIdx) {
            throw OpalException("FromFile::readFile",
                                "Line " + std::to_string(lineNumber) + " in '" + filename +
                                "' has fewer columns (" + std::to_string(values.size()) +
                                ") than required (index " + std::to_string(maxColIdx) + ").");
        }

        std::vector<double> phaseSpace(6);
        phaseSpace[0] = values[columnIndices_m[0]];
        phaseSpace[1] = values[columnIndices_m[1]];
        phaseSpace[2] = values[columnIndices_m[2]];
        phaseSpace[3] = values[columnIndices_m[3]];
        phaseSpace[4] = values[columnIndices_m[4]];
        phaseSpace[5] = values[columnIndices_m[5]];
        particleData_m.push_back(phaseSpace);
        ++lineNumber;
    }

    file.close();
    numParticles_m = particleData_m.size();

    if (numParticles_m == 0) {
        throw OpalException("FromFile::readFile",
                            "No valid particle data found in '" + filename + "'.");
    }

    if (numParticles_m != expectedNumParticles) {
        throw OpalException("FromFile::readFile",
                            "Number of data lines (" + std::to_string(numParticles_m) +
                            ") does not match declared count (" + std::to_string(expectedNumParticles) + ") in '" + filename + "'.");
    }
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
            throw OpalException("FromFile::parseHeader",
                                "Header must contain all column names (x, y, z, px, py, pz). "
                                "Missing or unrecognized: '" + std::string(names[i]) + "'.");
        }
    }

    return indices;
}

std::string FromFile::normalizeColumnName(const std::string& name) {
    std::string normalized = name;
    
    // Convert to lowercase
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    
    // Remove whitespace
    normalized.erase(std::remove_if(normalized.begin(), normalized.end(),
                                    [](unsigned char c) { return std::isspace(c); }),
                     normalized.end());
    
    return normalized;
}

void FromFile::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
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
        m << "Warning: File contains " << numParticles_m 
          << " particles, but " << numberOfParticles << " were requested." << endl;
        m << "Using " << numParticles_m << " particles from file." << endl;
        totalParticles = numParticles_m;
    }
    
    // If file has more particles than requested, use requested count
    if (numParticles_m > numberOfParticles) {
        totalParticles = numberOfParticles;
    }
    
    numberOfParticles = totalParticles;
    
    // Distribute particles across MPI ranks (capacity-aware)
    const int rank = ippl::Comm->rank();
    const int nranks = std::max(1, ippl::Comm->size());
    const size_t nranks_u = static_cast<size_t>(nranks);
    size_t nlocal = pc_m ? computeLocalEmitCount(totalParticles)
                         : (totalParticles / nranks_u
                            + (static_cast<size_t>(rank) < (totalParticles % nranks_u) ? 1 : 0));
    // Use Mpi scan to get the start index for each rank (since they are all potentially different!)
    size_t startIdx = 0;
    if (pc_m && ippl::Comm->size() > 0) {
        unsigned long nlocalUL = static_cast<unsigned long>(nlocal);
        unsigned long startIdxUL = 0;
        MPI_Exscan(&nlocalUL, &startIdxUL, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                   ippl::Comm->getCommunicator());
        if (rank > 0) {
            startIdx = static_cast<size_t>(startIdxUL);
        }
    }

    // Allocate particles, appending after any existing ones.
    const size_t nlocalCurrent = pc_m->getLocalNum();
    pc_m->create(nlocal);
    
    // Get Kokkos views
    view_type Rview = pc_m->R.getView();
    view_type Pview = pc_m->P.getView();
    
    // Create Kokkos view for particle data (host accessible)
    using HostView = Kokkos::View<double**, Kokkos::HostSpace>;
    HostView hostParticleData("hostParticleData", totalParticles, 6);
    
    // Copy data to host view
    for (size_t i = 0; i < totalParticles && i < particleData_m.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            hostParticleData(i, j) = particleData_m[i][j];
        }
    }
    
    // Create device mirror with layout compatible for Host->Device deep_copy
    // (required when no common execution space, e.g. Host vs Cuda)
    auto deviceParticleData =
        Kokkos::create_mirror(Kokkos::DefaultExecutionSpace::memory_space(), hostParticleData);
    Kokkos::deep_copy(deviceParticleData, hostParticleData);
    
    // Copy particle data into Kokkos views, appending into
    // [nlocalCurrent, nlocalCurrent + nlocal)
    Kokkos::parallel_for("FromFile::generateParticles", nlocal,
        KOKKOS_LAMBDA(const size_t k) {
            const size_t j       = nlocalCurrent + k;
            const size_t dataIdx = startIdx + k;
            Rview(j)[0] = deviceParticleData(dataIdx, 0); // x
            Rview(j)[1] = deviceParticleData(dataIdx, 1); // y
            Rview(j)[2] = deviceParticleData(dataIdx, 2); // z
            Pview(j)[0] = deviceParticleData(dataIdx, 3); // px
            Pview(j)[1] = deviceParticleData(dataIdx, 4); // py
            Pview(j)[2] = deviceParticleData(dataIdx, 5); // pz
        }
    );
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
    mALL << "FromFile: Loaded " << totalParticles << " particles from file '" 
         << filename_m << "'" << endl;
    ippl::Comm->barrier();
    mALL << "Rank " << rank << ": " << nlocal << " local particles" << endl;
    ippl::Comm->barrier();
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

    hasEmittedOnce_m = true;
    Vector_t<double, 3> dummyNr(0.0);
    generateParticles(requested, dummyNr);
}
