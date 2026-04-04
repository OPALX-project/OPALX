//
// Class DirichletPlaneWriter
//   Writes interpolated potential values on a 2D dirichlet plane to ASCII files.
//

#include "Structure/DirichletPlaneWriter.h"

#include "Ippl.h"
#include "Utilities/OpalException.h"

#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
std::string sanitizeTag(const std::string& tag) {
    std::string result = tag;
    for (char& c : result) {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (!std::isalnum(uc) && c != '-' && c != '_') {
            c = '_';
        }
    }
    return result.empty() ? std::string("solve") : result;
}
}  // namespace

DirichletPlaneWriter::DirichletPlaneWriter(const std::string& outputDirectory)
    : outputDirectory_m(outputDirectory) {
    Inform m("DirichletPlaneWriter::DirichletPlaneWriter");
    std::error_code ec;
    std::filesystem::create_directories(outputDirectory_m, ec);
    if (ec) {
        throw OpalException(
                "DirichletPlaneWriter::DirichletPlaneWriter",
                "Failed to create output directory \"" + outputDirectory_m.string() + "\": "
                        + ec.message());
    }

    m << level4 << "Dirichlet-plane diagnostics output directory: \""
      << outputDirectory_m.string() << "\"." << endl;
}

void DirichletPlaneWriter::writePlane(long long step,
                                      double time,
                                      double zPlane,
                                      const std::vector<double>& xCoords,
                                      const std::vector<double>& yCoords,
                                      const std::vector<double>& phiValues,
                                      std::size_t nx,
                                      std::size_t ny,
                                      const std::string& solveTag) {
    if (xCoords.size() != nx || yCoords.size() != ny || phiValues.size() != nx * ny) {
        throw OpalException(
                "DirichletPlaneWriter::writePlane",
                "Invalid plane dimensions for dirichlet-plane dump.");
    }

    const std::string safeTag = sanitizeTag(solveTag);

    std::ostringstream fileName;
    fileName << "phi_dirichlet_" << safeTag << "_step-" << std::setfill('0') << std::setw(8)
             << step << ".dat";

    const std::filesystem::path filePath = outputDirectory_m / fileName.str();
    std::ofstream out(filePath.string(), std::ios::out | std::ios::trunc);
    if (!out) {
        throw OpalException(
                "DirichletPlaneWriter::writePlane",
                "Failed to open output file \"" + filePath.string() + "\".");
    }

    out << std::setprecision(17);
    out << "# Dirichlet-plane potential dump\n";
    out << "# step=" << step << " time=" << time << " [s] zPlane=" << zPlane << " [m] tag="
        << safeTag << "\n";
    out << "# nx=" << nx << " ny=" << ny << "\n";
    out << "# columns: i j x[m] y[m] phi[V]\n";

    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) {
            const std::size_t idx = i * ny + j;
            out << std::setw(6) << i << std::setw(6) << j << " " << std::scientific
                << xCoords[i] << " " << yCoords[j] << " " << phiValues[idx] << "\n";
        }
    }
}
