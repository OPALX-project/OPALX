//
// Class BeamBeamWindowAnimation
//   ASCII visualization helper for the beam-beam-window dynamics.
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

#ifndef OPAL_BeamBeamWindowAnimation_HH
#define OPAL_BeamBeamWindowAnimation_HH

#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief ASCII visualization helper for the beam-beam-window dynamics.
 *
 * This utility renders a compact one-dimensional longitudinal view of the
 * interaction-point crossing. The rendered coordinate is the absolute path-length
 * coordinate `s`, not the bunch-local `z` coordinate.
 *
 * The visualization distinguishes two mesh regimes:
 * - moving bunch-following mesh: a short dotted envelope around the physical bunch
 * - frozen beam-beam-window mesh: a dotted envelope spanning the full
 *   beam-beam window
 *
 * The physical bunch is rendered as `*`. When the mirrored-beam model is active,
 * the virtual partner bunch is rendered as `+`, mirrored about the interaction
 * point.
 *
 * Output is written to:
 * - the terminal as an animated ASCII display
 * - a frame dump file, currently `data/collision_ascii_frames.txt`, for later GIF
 *   generation
 */
class BeamBeamWindowAnimation {
public:
    /**
     * @brief Construct the ASCII beam-beam-window animation helper.
     *
     * @param outputPath Path of the ASCII frame dump written for offline rendering.
     * @param width Number of ASCII columns used for the longitudinal display.
     */
    explicit BeamBeamWindowAnimation(
            const std::string& outputPath = "data/collision_ascii_frames.txt",
            std::size_t width             = 100);

    /**
     * @brief Render one beam-beam-window frame.
     *
     * All longitudinal positions are expressed in the absolute path-length
     * coordinate `s`.
     *
     * @param bunchCenterS Center of the physical bunch in `s`.
     * @param meshBeginS Begin of the mesh footprint shown in the frame.
     * @param meshEndS End of the mesh footprint shown in the frame.
     * @param beamBeamWindowBeginS Begin of the active beam-beam window.
     * @param beamBeamWindowEndS End of the active beam-beam window.
     * @param interactionPointS Center of the interaction point.
     * @param mirroredBeamActive True if the mirrored virtual bunch should be shown.
     * @param useFrozenBeamBeamWindowMesh True if the longitudinal mesh is fixed
     * over the beam-beam window, false if it is bunch-following.
     */
    void render(
            double bunchCenterS, double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
            double beamBeamWindowEndS, double interactionPointS, bool mirroredBeamActive,
            bool useFrozenBeamBeamWindowMesh);

private:
    /// Return the displayed `s` interval, including a symmetric margin around the beam-beam window.
    std::pair<double, double> getDisplayRange(
            double beamBeamWindowBeginS, double beamBeamWindowEndS) const;
    /// Test whether a longitudinal position lies inside the displayed `s` interval.
    bool isInsideDisplay(double z, double displayMinZ, double displayMaxZ) const;
    /// Map a longitudinal coordinate to an ASCII column index.
    std::size_t mapToDisplayIndex(double z, double displayMinZ, double displayMaxZ) const;
    /// Write a short label such as `S`, `IP`, or `E` into a line buffer.
    void placeLabel(std::string& line, std::size_t pos, const std::string& label) const;
    /// Return true if a character cell may be overwritten by a marker.
    bool isMarkerBackground(char ch) const;
    /// Place a bunch marker into a line buffer while preserving nearby frame delimiters.
    void placeMarker(std::string& line, std::size_t idx, char marker) const;
    /// Render the dotted mesh envelope shown above, and optionally below, the main axis.
    std::string renderMeshEnvelope(
            double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
            double beamBeamWindowEndS) const;
    /// Render the main axis with `S`, `IP`, `E`, the physical bunch, and the mirrored bunch.
    std::string renderCenters(
            double bunchCenterS, double meshBeginS, double meshEndS, double beamBeamWindowBeginS,
            double beamBeamWindowEndS, double interactionPointS, bool mirroredBeamActive) const;
    /// Render the label line for the beam-beam-window geometry.
    std::string renderLabels(
            double beamBeamWindowBeginS, double beamBeamWindowEndS, double interactionPointS,
            bool useFrozenBeamBeamWindowMesh) const;
    /// Update the live terminal animation.
    void writeToTerminal(const std::vector<std::string>& lines);
    /// Append one frame to the ASCII dump file used by the movie script.
    void writeToFile(const std::vector<std::string>& lines);

    /// Output stream for the persistent ASCII frame dump.
    std::ofstream dump_m;
    /// Number of ASCII columns in the rendered view.
    std::size_t width_m;
    /// True until the first frame has been written to the terminal.
    bool firstFrame_m;
    /// Sequential frame identifier used in the ASCII dump.
    std::size_t frameId_m;
    /// Number of terminal lines occupied by the previously rendered frame.
    std::size_t prevLineCount_m;
};

#endif
