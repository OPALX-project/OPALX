//
// Class CollisionDynamicsAnimation
//   ASCII visualization helper for the collision-window prototype.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//

#ifndef OPAL_CollisionDynamicsAnimation_HH
#define OPAL_CollisionDynamicsAnimation_HH

#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

class CollisionDynamicsAnimation {
public:
    explicit CollisionDynamicsAnimation(
        const std::string& outputPath = "data/collision_ascii_frames.txt",
        std::size_t width             = 100);

    void render(
        double bunchCenterS, double meshBeginS, double meshEndS, double windowBeginS,
        double windowEndS, double ipCenterS, bool mirroredActive, bool fixedMesh);

private:
    std::pair<double, double> getDisplayRange(double windowBeginS, double windowEndS) const;
    bool isInsideDisplay(double z, double displayMinZ, double displayMaxZ) const;
    std::size_t mapToDisplayIndex(
        double z, double displayMinZ, double displayMaxZ) const;
    void placeLabel(std::string& line, std::size_t pos, const std::string& label) const;
    bool isMarkerBackground(char ch) const;
    void placeMarker(std::string& line, std::size_t idx, char marker) const;
    std::string renderMeshEnvelope(
        double meshBeginS, double meshEndS, double windowBeginS, double windowEndS) const;
    std::string renderCenters(
        double bunchCenterS, double meshBeginS, double meshEndS, double windowBeginS,
        double windowEndS, double ipCenterS, bool mirroredActive) const;
    std::string renderLabels(
        double windowBeginS, double windowEndS, double ipCenterS, bool fixedMesh) const;
    void writeToTerminal(const std::vector<std::string>& lines);
    void writeToFile(const std::vector<std::string>& lines);

    std::ofstream dump_m;
    std::size_t width_m;
    bool firstFrame_m;
    std::size_t frameId_m;
    std::size_t prevLineCount_m;
};

#endif
