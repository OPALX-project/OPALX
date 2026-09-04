/**
 * @file ReferencePath.cpp
 * @brief Implements FFT2D5 design-path loading and device transfer.
 */

#include "SpaceCharge/FFT2D5/ReferencePath.h"

#include "Utilities/OpalException.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    ReferencePath ReferencePath::load(const std::string& filename) {
        if (filename.empty()) {
            throw OpalException(
                    "ReferencePath::load", "The FFT2D5 reference-path filename is empty.");
        }
        if (!std::filesystem::exists(filename)) {
            throw OpalException("ReferencePath::load", "File does not exist: " + filename);
        }

        std::ifstream file(filename);
        if (!file) {
            throw OpalException(
                    "ReferencePath::load", "Failed to open reference path: " + filename);
        }

        std::string line;
        std::getline(file, line);
        std::vector<Vector> points;
        double length = 0.0;
        while (std::getline(file, line)) {
            std::istringstream input(line);
            double s, rx, ry, rz, px, py, pz, ex, ey, ez, bx, by, bz, ke, time;
            input >> s >> rx >> ry >> rz >> px >> py >> pz >> ex >> ey >> ez >> bx >> by >> bz >> ke
                    >> time;
            if (input.fail()) {
                throw OpalException("ReferencePath::load", "Failed to parse line: " + line);
            }
            points.emplace_back();
            points.back()[0] = rx;
            points.back()[1] = ry;
            points.back()[2] = rz;
            if (points.size() > 1) {
                const Vector segment = points.back() - points[points.size() - 2];
                length += segment.Pnorm();
            }
        }
        if (points.size() < 2 || !(length > 0.0)) {
            throw OpalException(
                    "ReferencePath::load",
                    "FFT2D5 requires at least two distinct reference-path points.");
        }

        View deviceView("FFT2D5ReferencePath", points.size());
        auto hostView = Kokkos::create_mirror_view(deviceView);
        for (std::size_t index = 0; index < points.size(); ++index) {
            hostView(index) = points[index];
        }
        Kokkos::deep_copy(deviceView, hostView);
        return ReferencePath(std::move(deviceView), length);
    }

}  // namespace opalx::spacecharge
