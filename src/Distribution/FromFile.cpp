#include "Distribution.h"
#include "SamplingBase.hpp"
#include "FromFile.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <filesystem>

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;
using Distribution_t = Distribution;

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, 3>, 1>::view_type;

FromFile::FromFile(std::shared_ptr<ParticleContainer_t> &pc,
                   std::shared_ptr<FieldContainer_t> &fc,
                   std::shared_ptr<Distribution_t> &opalDist)
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

void FromFile::readFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw OpalException("FromFile::readFile",
                            "Couldn't open file '" + filename + "'.");
    }
    
    std::string line;
    std::vector<std::string> headerTokens;
    bool headerFound = false;
    bool firstLineIsNumber = false;
    size_t expectedNumParticles = 0;
    
    // Read first non-empty line
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) {
            continue;
        }
        
        // Check if first line is a number (alternative format)
        std::istringstream firstLineStream(line);
        double firstValue;
        if (firstLineStream >> firstValue && firstLineStream.eof()) {
            // First line is a number - this is the number of particles
            expectedNumParticles = static_cast<size_t>(firstValue);
            firstLineIsNumber = true;
            continue; // Skip to next line for header/data
        }
        
        // Parse as header line
        std::istringstream headerStream(line);
        std::string token;
        while (headerStream >> token) {
            headerTokens.push_back(token);
        }
        
        if (!headerTokens.empty()) {
            headerFound = true;
            columnIndices_m = parseHeader(line);
            break;
        }
    }
    
    if (!headerFound && !firstLineIsNumber) {
        throw OpalException("FromFile::readFile",
                            "Could not find header line in file '" + filename + "'.");
    }
    
    // If first line was a number, try to read header from next line
    if (firstLineIsNumber) {
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            if (line.empty()) {
                continue;
            }
            
            // Try to parse as header
            std::istringstream headerStream(line);
            std::string token;
            std::vector<std::string> tokens;
            while (headerStream >> token) {
                tokens.push_back(token);
            }
            
            // Check if this looks like a header (contains text) or data (all numbers)
            bool looksLikeHeader = false;
            for (const auto& t : tokens) {
                bool isNumber = true;
                for (char c : t) {
                    if (!std::isdigit(c) && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E') {
                        isNumber = false;
                        break;
                    }
                }
                if (!isNumber) {
                    looksLikeHeader = true;
                    break;
                }
            }
            
            if (looksLikeHeader) {
                columnIndices_m = parseHeader(line);
                headerFound = true;
                break;
            } else {
                // This is data, rewind and parse as data
                file.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
                break;
            }
        }
    }
    
    // If no header found but first line was number, assume standard column order
    if (!headerFound && firstLineIsNumber) {
        // Assume standard order: x y z px py pz
        columnIndices_m = {0, 1, 2, 3, 4, 5};
    }
    
    // Read data lines
    particleData_m.clear();
    size_t lineNumber = (firstLineIsNumber ? 2 : 2); // Account for header line
    
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) {
            ++lineNumber;
            continue;
        }
        
        std::istringstream lineStream(line);
        std::vector<double> values;
        double value;
        
        while (lineStream >> value) {
            values.push_back(value);
        }
        
        // Check if we have enough columns
        size_t maxColIdx = *std::max_element(columnIndices_m.begin(), columnIndices_m.end());
        if (values.size() <= maxColIdx) {
            throw OpalException("FromFile::readFile",
                                "Line " + std::to_string(lineNumber) + " in file '" + filename +
                                "' has fewer columns (" + std::to_string(values.size()) +
                                ") than required (index " + std::to_string(maxColIdx) + ").");
        }
        
        // Extract phase space coordinates
        std::vector<double> phaseSpace(6);
        phaseSpace[0] = values[columnIndices_m[0]]; // x
        phaseSpace[1] = values[columnIndices_m[1]]; // y
        phaseSpace[2] = values[columnIndices_m[2]]; // z
        phaseSpace[3] = values[columnIndices_m[3]]; // px
        phaseSpace[4] = values[columnIndices_m[4]]; // py
        phaseSpace[5] = values[columnIndices_m[5]]; // pz
        
        particleData_m.push_back(phaseSpace);
        ++lineNumber;
    }
    
    file.close();
    
    numParticles_m = particleData_m.size();
    
    if (numParticles_m == 0) {
        throw OpalException("FromFile::readFile",
                            "No valid particle data found in file '" + filename + "'.");
    }
    
    // If first line specified number of particles, verify
    if (firstLineIsNumber && expectedNumParticles > 0 && numParticles_m != expectedNumParticles) {
        throw OpalException("FromFile::readFile",
                            "Number of particles in file (" + std::to_string(numParticles_m) +
                            ") does not match expected number (" + std::to_string(expectedNumParticles) + ").");
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
    
    // Check if all required columns were found
    for (size_t i = 0; i < 6; ++i) {
        if (indices[i] == SIZE_MAX) {
            // If not found, assume sequential order starting from 0
            // This handles cases where header might be missing or use different names
            indices[i] = i;
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
    extern Inform* gmsg;
    
    // Use number of particles from file if available, otherwise use requested number
    size_t totalParticles = (numParticles_m > 0) ? numParticles_m : numberOfParticles;
    
    // If file has fewer particles than requested, use file count
    if (numParticles_m > 0 && numParticles_m < numberOfParticles) {
        *gmsg << "* Warning: File contains " << numParticles_m 
              << " particles, but " << numberOfParticles << " were requested." << endl;
        *gmsg << "* Using " << numParticles_m << " particles from file." << endl;
        totalParticles = numParticles_m;
    }
    
    // If file has more particles than requested, use requested count
    if (numParticles_m > numberOfParticles) {
        totalParticles = numberOfParticles;
    }
    
    numberOfParticles = totalParticles;
    
    // Distribute particles across MPI ranks
    MPI_Comm comm = MPI_COMM_WORLD;
    int nranks, rank;
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &rank);
    
    size_t baseParticlesPerRank = totalParticles / nranks;
    size_t remaining = totalParticles - baseParticlesPerRank * nranks;
    
    // Calculate starting index for this rank
    // Rank 0 gets remaining particles, so other ranks start after that
    size_t startIdx = rank * baseParticlesPerRank + remaining;
    
    // Adjust nlocal: rank 0 gets remaining particles
    size_t nlocal = baseParticlesPerRank;
    if (remaining > 0 && rank == 0) {
        nlocal += remaining;
        startIdx = 0; // Rank 0 always starts at 0
    }
    
    // Allocate particles
    pc_m->create(nlocal);
    
    // Get Kokkos views
    view_type &Rview = pc_m->R.getView();
    view_type &Pview = pc_m->P.getView();
    
    // Create Kokkos view for particle data (host accessible)
    using HostView = Kokkos::View<double**, Kokkos::HostSpace>;
    HostView hostParticleData("hostParticleData", totalParticles, 6);
    
    // Copy data to host view
    for (size_t i = 0; i < totalParticles && i < particleData_m.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            hostParticleData(i, j) = particleData_m[i][j];
        }
    }
    
    // Create device view and copy from host
    using DeviceView = Kokkos::View<double**, Kokkos::DefaultExecutionSpace::memory_space>;
    DeviceView deviceParticleData("deviceParticleData", totalParticles, 6);
    Kokkos::deep_copy(deviceParticleData, hostParticleData);
    
    // Copy particle data into Kokkos views
    Kokkos::parallel_for(
        "FromFile::generateParticles", nlocal,
        KOKKOS_LAMBDA(const size_t k) {
            size_t dataIdx = startIdx + k;
            if (dataIdx < totalParticles) {
                Rview(k)[0] = deviceParticleData(dataIdx, 0); // x
                Rview(k)[1] = deviceParticleData(dataIdx, 1); // y
                Rview(k)[2] = deviceParticleData(dataIdx, 2); // z
                Pview(k)[0] = deviceParticleData(dataIdx, 3); // px
                Pview(k)[1] = deviceParticleData(dataIdx, 4); // py
                Pview(k)[2] = deviceParticleData(dataIdx, 5); // pz
            }
        }
    );
    Kokkos::fence();
    
    *gmsg << "* FromFile: Loaded " << totalParticles << " particles from file '" 
          << filename_m << "'" << endl;
    *gmsg << "* Rank " << rank << ": " << nlocal << " local particles" << endl;
}
