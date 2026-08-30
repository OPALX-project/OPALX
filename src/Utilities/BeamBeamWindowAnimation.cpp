#include "Utilities/BeamBeamWindowAnimation.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

BeamBeamWindowAnimation::BeamBeamWindowAnimation(const std::string& outputPath, std::size_t width)
    : dump_m(outputPath),
      width_m(std::max<std::size_t>(width, 3)),
      firstFrame_m(true),
      frameId_m(0),
      prevLineCount_m(0) {}

void BeamBeamWindowAnimation::render(
        double bunchCenterS, double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
        double beamBeamWindowEndS, double interactionPointS, bool mirroredBeamActive,
        bool useFrozenBeamBeamWindowMesh) {
    std::vector<std::string> lines{
            renderMeshEnvelope(meshBeginS, meshEndS, beamBeamWindowBeginS, beamBeamWindowEndS),
            renderCenters(
                    bunchCenterS, meshBeginS, meshEndS, beamBeamWindowBeginS, beamBeamWindowEndS,
                    interactionPointS, mirroredBeamActive),
    };

    if (!useFrozenBeamBeamWindowMesh) {
        lines.push_back(
                renderMeshEnvelope(meshBeginS, meshEndS, beamBeamWindowBeginS, beamBeamWindowEndS));
    }

    lines.push_back(renderLabels(
            beamBeamWindowBeginS, beamBeamWindowEndS, interactionPointS,
            useFrozenBeamBeamWindowMesh));

    writeToTerminal(lines);
    writeToFile(lines);
}

std::pair<double, double> BeamBeamWindowAnimation::getDisplayRange(
        double beamBeamWindowBeginS, double beamBeamWindowEndS) const {
    const double beamBeamWindowLength = beamBeamWindowEndS - beamBeamWindowBeginS;
    const double margin               = 0.5 * beamBeamWindowLength;
    return {
            beamBeamWindowBeginS - margin,
            beamBeamWindowEndS + margin,
    };
}

bool BeamBeamWindowAnimation::isInsideDisplay(
        double z, double displayMinZ, double displayMaxZ) const {
    return z >= displayMinZ && z <= displayMaxZ;
}

std::size_t BeamBeamWindowAnimation::mapToDisplayIndex(
        double z, double displayMinZ, double displayMaxZ) const {
    const double span = displayMaxZ - displayMinZ;
    if (span <= 0.0 || width_m < 3) {
        return 1;
    }

    double u = (z - displayMinZ) / span;
    u        = std::clamp(u, 0.0, 1.0);

    std::size_t idx = static_cast<std::size_t>(std::llround(u * static_cast<double>(width_m - 1)));
    return std::clamp<std::size_t>(idx, 1, width_m - 2);
}

void BeamBeamWindowAnimation::placeLabel(
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

bool BeamBeamWindowAnimation::isMarkerBackground(char ch) const {
    return ch == '-' || ch == ' ' || ch == '.';
}

void BeamBeamWindowAnimation::placeMarker(std::string& line, std::size_t idx, char marker) const {
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

std::string BeamBeamWindowAnimation::renderMeshEnvelope(
        double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
        double beamBeamWindowEndS) const {
    const auto [displayMinS, displayMaxS] =
            getDisplayRange(beamBeamWindowBeginS, beamBeamWindowEndS);

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

std::string BeamBeamWindowAnimation::renderCenters(
        double bunchCenterS, double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
        double beamBeamWindowEndS, double interactionPointS, bool mirroredBeamActive) const {
    const auto [displayMinS, displayMaxS] =
            getDisplayRange(beamBeamWindowBeginS, beamBeamWindowEndS);

    std::string line(width_m, '-');
    line.front() = '[';
    line.back()  = ']';

    const std::size_t startIdx = mapToDisplayIndex(beamBeamWindowBeginS, displayMinS, displayMaxS);
    const std::size_t endIdx   = mapToDisplayIndex(beamBeamWindowEndS, displayMinS, displayMaxS);

    line[startIdx] = '[';
    line[endIdx]   = ']';

    const std::size_t ipIdx = mapToDisplayIndex(interactionPointS, displayMinS, displayMaxS);
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

    if (mirroredBeamActive) {
        const double mirroredCenterS = 2.0 * interactionPointS - bunchCenterS;

        if (isInsideDisplay(mirroredCenterS, displayMinS, displayMaxS)) {
            placeMarker(line, mapToDisplayIndex(mirroredCenterS, displayMinS, displayMaxS), '+');
        } else if (mirroredCenterS < displayMinS) {
            placeMarker(line, 1, '<');
        } else {
            placeMarker(line, width_m - 2, '>');
        }
    }

    return line;
}

std::string BeamBeamWindowAnimation::renderLabels(
        double beamBeamWindowBeginS, double beamBeamWindowEndS, double interactionPointS,
        bool useFrozenBeamBeamWindowMesh) const {
    const auto [displayMinS, displayMaxS] =
            getDisplayRange(beamBeamWindowBeginS, beamBeamWindowEndS);

    std::string labels(width_m, ' ');

    const std::size_t startIdx = mapToDisplayIndex(beamBeamWindowBeginS, displayMinS, displayMaxS);
    const std::size_t ipIdx    = mapToDisplayIndex(interactionPointS, displayMinS, displayMaxS);
    const std::size_t endIdx   = mapToDisplayIndex(beamBeamWindowEndS, displayMinS, displayMaxS);

    placeLabel(labels, startIdx, "S");
    placeLabel(labels, ipIdx, "IP");
    placeLabel(labels, endIdx, "E");

    if (useFrozenBeamBeamWindowMesh) {
        for (std::size_t idx = startIdx + 1; idx < endIdx; ++idx) {
            if (labels[idx] == ' ') {
                labels[idx] = '.';
            }
        }
    }

    return labels;
}

void BeamBeamWindowAnimation::writeToTerminal(const std::vector<std::string>& lines) {
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

void BeamBeamWindowAnimation::writeToFile(const std::vector<std::string>& lines) {
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
