#include "Utilities/CollisionDynamicsAnimation.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

CollisionDynamicsAnimation::CollisionDynamicsAnimation(
    const std::string& outputPath, std::size_t width)
    : dump_m(outputPath),
      width_m(std::max<std::size_t>(width, 3)),
      firstFrame_m(true),
      frameId_m(0),
      prevLineCount_m(0) {
}

void CollisionDynamicsAnimation::render(
    double bunchCenterS, double meshBeginS, double meshEndS, double windowBeginS,
    double windowEndS, double ipCenterS, bool mirroredActive, bool fixedMesh) {
    std::vector<std::string> lines{
        renderMeshEnvelope(meshBeginS, meshEndS, windowBeginS, windowEndS),
        renderCenters(
            bunchCenterS, meshBeginS, meshEndS, windowBeginS, windowEndS, ipCenterS,
            mirroredActive),
    };

    if (!fixedMesh) {
        lines.push_back(renderMeshEnvelope(meshBeginS, meshEndS, windowBeginS, windowEndS));
    }

    lines.push_back(renderLabels(windowBeginS, windowEndS, ipCenterS, fixedMesh));

    writeToTerminal(lines);
    writeToFile(lines);
}

std::pair<double, double> CollisionDynamicsAnimation::getDisplayRange(
    double windowBeginS, double windowEndS) const {
    const double windowLength = windowEndS - windowBeginS;
    const double margin       = 0.5 * windowLength;
    return {windowBeginS - margin, windowEndS + margin};
}

bool CollisionDynamicsAnimation::isInsideDisplay(
    double z, double displayMinZ, double displayMaxZ) const {
    return z >= displayMinZ && z <= displayMaxZ;
}

std::size_t CollisionDynamicsAnimation::mapToDisplayIndex(
    double z, double displayMinZ, double displayMaxZ) const {
    const double span = displayMaxZ - displayMinZ;
    if (span <= 0.0 || width_m < 3) {
        return 1;
    }

    double u = (z - displayMinZ) / span;
    u        = std::clamp(u, 0.0, 1.0);

    std::size_t idx =
        static_cast<std::size_t>(std::llround(u * static_cast<double>(width_m - 1)));
    return std::clamp<std::size_t>(idx, 1, width_m - 2);
}

void CollisionDynamicsAnimation::placeLabel(
    std::string& line, std::size_t pos, const std::string& label) const {
    if (line.empty() || label.empty()) {
        return;
    }

    if (pos >= line.size()) {
        pos = line.size() - 1;
    }

    const std::size_t maxLen = std::min(label.size(), line.size() - pos);
    for (std::size_t i = 0; i < maxLen; ++i) {
        line[pos + i] = label[i];
    }
}

bool CollisionDynamicsAnimation::isMarkerBackground(char ch) const {
    return ch == '-' || ch == ' ' || ch == '.';
}

void CollisionDynamicsAnimation::placeMarker(
    std::string& line, std::size_t idx, char marker) const {
    if (idx >= line.size()) {
        return;
    }

    if (isMarkerBackground(line[idx])) {
        line[idx] = marker;
        return;
    }

    if (line[idx] == '[') {
        if (idx + 1 < line.size() && isMarkerBackground(line[idx + 1])) {
            line[idx + 1] = marker;
            return;
        }
        if (idx > 0 && isMarkerBackground(line[idx - 1])) {
            line[idx - 1] = marker;
            return;
        }
        return;
    }

    if (line[idx] == ']' || line[idx] == '|') {
        if (idx > 0 && isMarkerBackground(line[idx - 1])) {
            line[idx - 1] = marker;
            return;
        }
        if (idx + 1 < line.size() && isMarkerBackground(line[idx + 1])) {
            line[idx + 1] = marker;
            return;
        }
        return;
    }

    line[idx] = 'X';
}

std::string CollisionDynamicsAnimation::renderMeshEnvelope(
    double meshBeginS, double meshEndS, double windowBeginS, double windowEndS) const {
    const auto [displayMinS, displayMaxS] = getDisplayRange(windowBeginS, windowEndS);

    std::string line(width_m, ' ');

    const std::size_t meshBeginIdx = mapToDisplayIndex(meshBeginS, displayMinS, displayMaxS);
    const std::size_t meshEndIdx   = mapToDisplayIndex(meshEndS, displayMinS, displayMaxS);
    const std::size_t meshLo       = std::min(meshBeginIdx, meshEndIdx);
    const std::size_t meshHi       = std::max(meshBeginIdx, meshEndIdx);

    for (std::size_t idx = meshLo; idx <= meshHi; ++idx) {
        line[idx] = '.';
    }

    return line;
}

std::string CollisionDynamicsAnimation::renderCenters(
    double bunchCenterS, double meshBeginS, double meshEndS, double windowBeginS,
    double windowEndS, double ipCenterS, bool mirroredActive) const {
    const auto [displayMinS, displayMaxS] = getDisplayRange(windowBeginS, windowEndS);

    std::string line(width_m, '-');
    line.front() = '[';
    line.back()  = ']';

    const std::size_t startIdx = mapToDisplayIndex(windowBeginS, displayMinS, displayMaxS);
    const std::size_t endIdx   = mapToDisplayIndex(windowEndS, displayMinS, displayMaxS);

    line[startIdx] = '[';
    line[endIdx]   = ']';

    const std::size_t ipIdx = mapToDisplayIndex(ipCenterS, displayMinS, displayMaxS);
    line[ipIdx]             = '|';

    const std::size_t meshBeginIdx = mapToDisplayIndex(meshBeginS, displayMinS, displayMaxS);
    const std::size_t meshEndIdx   = mapToDisplayIndex(meshEndS, displayMinS, displayMaxS);
    const std::size_t meshLo       = std::min(meshBeginIdx, meshEndIdx);
    const std::size_t meshHi       = std::max(meshBeginIdx, meshEndIdx);

    for (std::size_t idx = meshLo; idx <= meshHi; ++idx) {
        if (line[idx] == '-') {
            line[idx] = '.';
        }
    }

    if (isInsideDisplay(bunchCenterS, displayMinS, displayMaxS)) {
        placeMarker(line, mapToDisplayIndex(bunchCenterS, displayMinS, displayMaxS), '*');
    } else if (bunchCenterS < displayMinS) {
        placeMarker(line, 1, '<');
    } else {
        placeMarker(line, width_m - 2, '>');
    }

    if (mirroredActive) {
        const double mirroredCenterS = 2.0 * ipCenterS - bunchCenterS;

        if (isInsideDisplay(mirroredCenterS, displayMinS, displayMaxS)) {
            placeMarker(
                line, mapToDisplayIndex(mirroredCenterS, displayMinS, displayMaxS), '+');
        } else if (mirroredCenterS < displayMinS) {
            placeMarker(line, 1, '<');
        } else {
            placeMarker(line, width_m - 2, '>');
        }
    }

    return line;
}

std::string CollisionDynamicsAnimation::renderLabels(
    double windowBeginS, double windowEndS, double ipCenterS, bool fixedMesh) const {
    const auto [displayMinS, displayMaxS] = getDisplayRange(windowBeginS, windowEndS);

    std::string labels(width_m, ' ');

    const std::size_t startIdx = mapToDisplayIndex(windowBeginS, displayMinS, displayMaxS);
    const std::size_t ipIdx    = mapToDisplayIndex(ipCenterS, displayMinS, displayMaxS);
    const std::size_t endIdx   = mapToDisplayIndex(windowEndS, displayMinS, displayMaxS);

    placeLabel(labels, startIdx, "S");
    placeLabel(labels, ipIdx, "IP");
    placeLabel(labels, endIdx, "E");

    if (fixedMesh) {
        for (std::size_t idx = startIdx + 1; idx < endIdx; ++idx) {
            if (labels[idx] == ' ') {
                labels[idx] = '.';
            }
        }
    }

    return labels;
}

void CollisionDynamicsAnimation::writeToTerminal(const std::vector<std::string>& lines) {
    if (firstFrame_m) {
        for (std::size_t idx = 0; idx < lines.size(); ++idx) {
            std::cout << lines[idx];
            if (idx + 1 < lines.size()) {
                std::cout << '\n';
            }
        }
        firstFrame_m = false;
    } else {
        for (std::size_t idx = 1; idx < prevLineCount_m; ++idx) {
            std::cout << "\033[1A";
        }
        for (std::size_t idx = 0; idx < prevLineCount_m; ++idx) {
            std::cout << "\r\033[2K";
            if (idx < lines.size()) {
                std::cout << lines[idx];
            }
            if (idx + 1 < prevLineCount_m) {
                std::cout << '\n';
            }
        }
    }

    prevLineCount_m = lines.size();
    std::cout << std::flush;
}

void CollisionDynamicsAnimation::writeToFile(const std::vector<std::string>& lines) {
    if (!dump_m.is_open()) {
        return;
    }

    dump_m << "FRAME " << frameId_m++ << '\n';
    for (const auto& line : lines) {
        dump_m << line << '\n';
    }
    dump_m << '\n';
    dump_m.flush();
}
