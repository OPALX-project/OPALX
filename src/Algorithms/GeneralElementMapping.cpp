//
// Class GeneralElementMapping
//   Uniform lab-frame grid that maps a bunch box to the elements it can reach.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "Algorithms/GeneralElementMapping.h"

#include "BeamlineGeometry/Geometry.h"
#include "Structure/BoundingBox.h"
#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <tuple>

namespace {
    SupportBox axisAligned(const Vector_t<double, 3>& lower, const Vector_t<double, 3>& upper) {
        SupportBox box;
        box.centre = 0.5 * (lower + upper);
        box.half   = 0.5 * (upper - lower);
        for (unsigned i = 0; i < 3; ++i) {
            box.axes[i]    = Vector_t<double, 3>(0.0);
            box.axes[i](i) = 1.0;
        }
        return box;
    }
}  // namespace

GeneralElementMapping::GeneralElementMapping(
        const ElementList& elements, const std::array<unsigned, 3>& cells, const double reach)
    : cells_m(cells), meshLower_m(0.0), meshUpper_m(0.0), cellSize_m(1.0) {
    const std::size_t nCells = std::size_t(cells[0]) * cells[1] * cells[2];
    if (nCells == 0 || nCells > maximumCells) {
        throw OpalException(
                "GeneralElementMapping",
                "ELEMENTMAPPINGGRID needs at least one cell per axis "
                "and at most 2^24 cells in total.");
    }
    if (!std::isfinite(reach) || reach <= 0.0) {
        throw OpalException("GeneralElementMapping", "ELEMENTMAPPINGREACH must be positive.");
    }
    if (elements.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw OpalException("GeneralElementMapping", "Too many elements for the element grid.");
    }

    // Serial snapshot: some getters set up state lazily (ScalingFFAMagnet::getFieldExtent).
    entries_m.reserve(elements.size());
    for (const auto& element : elements) {
        Entry entry;
        entry.element    = element;
        entry.reach      = element->getReach().value_or(reach);
        entry.labToLocal = element->getMisalignment() * element->getCSTrafoGlobal2Local();
        entry.localBoxes = element->getLocalSupportBoxes(entry.reach);
        for (const SupportBox& local : entry.localBoxes)
            entry.boxes.push_back(placeInLab(local, entry.labToLocal));

        // Same window as the element's markOutsideAperture(): sector bends and curved
        // MULTIPOLET use the arc chart from 0, all others the straight body from startZ.
        const Geometry& geometry = element->getGeometry();
        const ElementType type   = element->getType();
        const bool arc           = (type == ElementType::SBEND || type == ElementType::MULTIPOLET)
                         && geometry.kind() == GeometryKind::SBend
                         && geometry.getCurvature() != 0.0;
        const auto aperture = element->getAperture();
        entry.deletes       = element->getFlagDeleteOnTransverseExit() && type != ElementType::BOX;
        entry.arcCurvature  = arc ? geometry.getCurvature() : 0.0;
        entry.bodyBegin     = arc ? 0.0 : geometry.getStartZ();
        entry.bodyEnd       = entry.bodyBegin + geometry.getElementLength();
        entry.apertureType  = aperture.first;
        entry.apertureX     = aperture.second[0];
        entry.apertureY     = aperture.second[1];
        entryOf_m.emplace(element.get(), entries_m.size());
        entries_m.push_back(entry);
    }

    cellOffsets_m.assign(nCells + 1, 0);
    if (entries_m.empty()) return;

    BoundingBox mesh;
    for (const Entry& entry : entries_m) {
        for (const SupportBox& box : entry.boxes) {
            const Vector_t<double, 3> half = labHalfExtent(box);
            mesh.enlargeToContainPosition(box.centre - half);
            mesh.enlargeToContainPosition(box.centre + half);
        }
    }
    std::tie(meshLower_m, meshUpper_m) = mesh.getCorners();
    for (unsigned d = 0; d < 3; ++d) {
        if (!(meshUpper_m(d) > meshLower_m(d))) {
            meshLower_m(d) -= 0.5;
            meshUpper_m(d) += 0.5;
        }
        cellSize_m(d) = (meshUpper_m(d) - meshLower_m(d)) / cells_m[d];
    }

    // Each element writes only its own list, so no synchronization is needed.
    std::vector<std::vector<std::uint32_t>> cellsOfEntry(entries_m.size());
    Kokkos::parallel_for(
            "GeneralElementMapping::rasterize",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, entries_m.size()),
            [&](const std::size_t e) {
                std::vector<std::uint32_t>& cells = cellsOfEntry[e];
                for (const SupportBox& box : entries_m[e].boxes) {
                    const std::vector<std::uint32_t> overlapped = overlappedCells(box);
                    cells.insert(cells.end(), overlapped.begin(), overlapped.end());
                }
                std::sort(cells.begin(), cells.end());
                cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
            });
    Kokkos::DefaultHostExecutionSpace().fence();

    // Filling in entry order keeps each row sorted, independently of the thread schedule.
    for (const auto& cellsOf : cellsOfEntry) {
        for (const std::uint32_t cell : cellsOf)
            ++cellOffsets_m[cell + 1];
    }
    std::partial_sum(cellOffsets_m.begin(), cellOffsets_m.end(), cellOffsets_m.begin());
    cellElements_m.resize(cellOffsets_m.back());
    std::vector<std::size_t> cursor(cellOffsets_m.begin(), cellOffsets_m.end() - 1);
    for (std::size_t e = 0; e < cellsOfEntry.size(); ++e) {
        for (const std::uint32_t cell : cellsOfEntry[e])
            cellElements_m[cursor[cell]++] = static_cast<std::uint32_t>(e);
    }
}

GeneralElementMapping::ElementSet GeneralElementMapping::query(
        const Vector_t<double, 3>& lower, const Vector_t<double, 3>& upper,
        const CoordinateSystemTrafo& toLab) const {
    ElementSet active;
    if (entries_m.empty()) return active;

    // The inverse of bunch-to-lab plays the role of an element's lab-to-local transform.
    const SupportBox bunch         = placeInLab(axisAligned(lower, upper), toLab.inverted());
    const Vector_t<double, 3> half = labHalfExtent(bunch);
    std::array<unsigned, 3> first{}, last{};
    if (!cellRange(bunch.centre - half, bunch.centre + half, first, last)) return active;

    std::vector<std::uint32_t> candidates;
    for (unsigned i = first[0]; i <= last[0]; ++i) {
        for (unsigned j = first[1]; j <= last[1]; ++j) {
            for (unsigned k = first[2]; k <= last[2]; ++k) {
                const std::size_t cell = flatIndex(i, j, k);
                candidates.insert(
                        candidates.end(), cellElements_m.begin() + cellOffsets_m[cell],
                        cellElements_m.begin() + cellOffsets_m[cell + 1]);
            }
        }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    for (const std::uint32_t e : candidates) {
        for (const SupportBox& box : entries_m[e].boxes) {
            if (intersect(bunch, box)) {
                active.insert(entries_m[e].element);
                break;
            }
        }
    }
    return active;
}

std::vector<GeneralElementMapping::SupportOutsideBox> GeneralElementMapping::checkReferencePoint(
        const Vector_t<double, 3>& position, const double pathLength) const {
    std::vector<SupportOutsideBox> outside;
    for (const Entry& entry : entries_m) {
        const Vector_t<double, 3> local = entry.labToLocal.transformTo(position);
        const Vector_t<double, 3> chart =
                entry.arcCurvature != 0.0
                        ? GeometryHelper::toBendArcCoords(
                                  local, entry.arcCurvature, entry.bodyEnd - entry.bodyBegin)
                        : local;
        const double offset = std::max(std::abs(chart(0)), std::abs(chart(1)));

        if (!contains(entry.localBoxes, local)) {
            if (entry.element->isInside(local))
                outside.push_back({entry.element, offset, entry.reach});
            continue;
        }

        const bool inBody = chart(2) >= entry.bodyBegin && chart(2) < entry.bodyEnd;
        if (entry.deletes && inBody
            && !ApertureHelper::isInsideAperture(
                    chart(0), chart(1), entry.apertureType, entry.apertureX, entry.apertureY)) {
            // Whether the box would still contain the point with a vanishing reach.
            const bool reachCanClear =
                    !contains(entry.element->getLocalSupportBoxes(1e-9 * entry.reach), local);
            std::ostringstream message;
            message << "The reference orbit at s = " << pathLength << " m passes element '"
                    << entry.element->getName() << "' at transverse position (" << chart(0) << ", "
                    << chart(1) << ") m, inside its box (reach " << entry.reach
                    << " m) but outside its aperture, so a particle on the reference orbit would "
                       "be deleted. ";
            if (reachCanClear) {
                message << "If the element belongs to another branch, set its REACH below "
                        << offset
                        << " m or lower OPTION ELEMENTMAPPINGREACH. Otherwise the reference orbit "
                           "hits the element's aperture.";
            } else {
                message << "REACH cannot exclude this point: the element's own extent (aperture "
                           "box, field region or bend orbit) covers it. Either the reference orbit "
                           "hits the element's aperture, or, for an element of another branch, "
                           "set DELETEONTRANSVERSEEXIT=FALSE on it.";
            }
            throw OpalException("GeneralElementMapping::checkReferencePoint", message.str());
        }
    }
    return outside;
}

bool GeneralElementMapping::coversTransversely(
        const ElementBase& element, const Vector_t<double, 3>& position) const {
    const auto found = entryOf_m.find(&element);
    if (found == entryOf_m.end()) {
        throw OpalException(
                "GeneralElementMapping::coversTransversely",
                "Element '" + element.getName() + "' is not mapped.");
    }
    // Annular sectors have no transverse axis: their box spans radius and angle, and the
    // reach does not apply to them.
    if (element.getType() == ElementType::CYCLOTRONSECTOR
        || element.getType() == ElementType::SCALINGFFAMAGNET)
        return true;
    const Entry& entry              = entries_m[found->second];
    const Vector_t<double, 3> local = entry.labToLocal.transformTo(position);
    for (const SupportBox& box : entry.localBoxes) {
        const Vector_t<double, 3> d = local - box.centre;
        if (std::abs(dot(d, box.axes[0])) <= box.half(0)
            && std::abs(dot(d, box.axes[1])) <= box.half(1))
            return true;
    }
    return false;
}

SupportBox GeneralElementMapping::placeInLab(
        const SupportBox& local, const CoordinateSystemTrafo& labToLocal) {
    SupportBox box;
    box.centre = labToLocal.transformFrom(local.centre);
    box.half   = local.half;
    for (unsigned i = 0; i < 3; ++i)
        box.axes[i] = labToLocal.rotateFrom(local.axes[i]);
    return box;
}

bool GeneralElementMapping::contains(
        const std::vector<SupportBox>& boxes, const Vector_t<double, 3>& point) {
    return std::any_of(boxes.begin(), boxes.end(), [&point](const SupportBox& box) {
        return box.contains(point);
    });
}

Vector_t<double, 3> GeneralElementMapping::labHalfExtent(const SupportBox& box) {
    Vector_t<double, 3> half(0.0);
    for (unsigned d = 0; d < 3; ++d) {
        for (unsigned i = 0; i < 3; ++i)
            half(d) += std::abs(box.axes[i](d)) * box.half(i);
    }
    return half;
}

bool GeneralElementMapping::intersect(const SupportBox& a, const SupportBox& b) {
    // Separating-axis test for two oriented boxes (Ericson, Real-Time Collision Detection,
    // section 4.4.1): the three axes of each box and their nine cross products. The epsilon
    // keeps near-parallel axes, common on straight lattices, from producing false separations.
    constexpr double epsilon = 1e-9;
    double R[3][3], absR[3][3];
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            R[i][j]    = dot(a.axes[i], b.axes[j]);
            absR[i][j] = std::abs(R[i][j]) + epsilon;
        }
    }
    const Vector_t<double, 3> d = b.centre - a.centre;
    const double t[3]           = {dot(d, a.axes[0]), dot(d, a.axes[1]), dot(d, a.axes[2])};
    const auto& ea              = a.half;
    const auto& eb              = b.half;

    for (unsigned i = 0; i < 3; ++i) {
        const double rb = eb(0) * absR[i][0] + eb(1) * absR[i][1] + eb(2) * absR[i][2];
        if (std::abs(t[i]) > ea(i) + rb) return false;
    }
    for (unsigned j = 0; j < 3; ++j) {
        const double ra = ea(0) * absR[0][j] + ea(1) * absR[1][j] + ea(2) * absR[2][j];
        if (std::abs(t[0] * R[0][j] + t[1] * R[1][j] + t[2] * R[2][j]) > ra + eb(j)) return false;
    }

    const auto separated = [](double distance, double ra, double rb) {
        return std::abs(distance) > ra + rb;
    };
    // a0 x b0, a0 x b1, a0 x b2
    if (separated(
                t[2] * R[1][0] - t[1] * R[2][0], ea(1) * absR[2][0] + ea(2) * absR[1][0],
                eb(1) * absR[0][2] + eb(2) * absR[0][1]))
        return false;
    if (separated(
                t[2] * R[1][1] - t[1] * R[2][1], ea(1) * absR[2][1] + ea(2) * absR[1][1],
                eb(0) * absR[0][2] + eb(2) * absR[0][0]))
        return false;
    if (separated(
                t[2] * R[1][2] - t[1] * R[2][2], ea(1) * absR[2][2] + ea(2) * absR[1][2],
                eb(0) * absR[0][1] + eb(1) * absR[0][0]))
        return false;
    // a1 x b0, a1 x b1, a1 x b2
    if (separated(
                t[0] * R[2][0] - t[2] * R[0][0], ea(0) * absR[2][0] + ea(2) * absR[0][0],
                eb(1) * absR[1][2] + eb(2) * absR[1][1]))
        return false;
    if (separated(
                t[0] * R[2][1] - t[2] * R[0][1], ea(0) * absR[2][1] + ea(2) * absR[0][1],
                eb(0) * absR[1][2] + eb(2) * absR[1][0]))
        return false;
    if (separated(
                t[0] * R[2][2] - t[2] * R[0][2], ea(0) * absR[2][2] + ea(2) * absR[0][2],
                eb(0) * absR[1][1] + eb(1) * absR[1][0]))
        return false;
    // a2 x b0, a2 x b1, a2 x b2
    if (separated(
                t[1] * R[0][0] - t[0] * R[1][0], ea(0) * absR[1][0] + ea(1) * absR[0][0],
                eb(1) * absR[2][2] + eb(2) * absR[2][1]))
        return false;
    if (separated(
                t[1] * R[0][1] - t[0] * R[1][1], ea(0) * absR[1][1] + ea(1) * absR[0][1],
                eb(0) * absR[2][2] + eb(2) * absR[2][0]))
        return false;
    if (separated(
                t[1] * R[0][2] - t[0] * R[1][2], ea(0) * absR[1][2] + ea(1) * absR[0][2],
                eb(0) * absR[2][1] + eb(1) * absR[2][0]))
        return false;
    return true;
}

bool GeneralElementMapping::cellRange(
        const Vector_t<double, 3>& lower, const Vector_t<double, 3>& upper,
        std::array<unsigned, 3>& first, std::array<unsigned, 3>& last) const {
    for (unsigned d = 0; d < 3; ++d) {
        if (upper(d) < meshLower_m(d) || lower(d) > meshUpper_m(d)) return false;
        // Clamp in floating point before the cast, so far-away or huge boxes stay in range.
        const double top     = cells_m[d] - 1;
        const auto cellIndex = [&](double x) {
            return static_cast<unsigned>(
                    std::clamp(std::floor((x - meshLower_m(d)) / cellSize_m(d)), 0.0, top));
        };
        first[d] = cellIndex(lower(d));
        last[d]  = cellIndex(upper(d));
    }
    return true;
}

std::size_t GeneralElementMapping::flatIndex(unsigned i, unsigned j, unsigned k) const {
    return (std::size_t(i) * cells_m[1] + j) * cells_m[2] + k;
}

std::vector<std::uint32_t> GeneralElementMapping::overlappedCells(const SupportBox& box) const {
    std::vector<std::uint32_t> overlapped;
    const Vector_t<double, 3> half = labHalfExtent(box);
    std::array<unsigned, 3> first{}, last{};
    if (!cellRange(box.centre - half, box.centre + half, first, last)) return overlapped;

    SupportBox cell;
    cell.half = 0.5 * cellSize_m;
    for (unsigned i = 0; i < 3; ++i) {
        cell.axes[i]    = Vector_t<double, 3>(0.0);
        cell.axes[i](i) = 1.0;
    }
    for (unsigned i = first[0]; i <= last[0]; ++i) {
        for (unsigned j = first[1]; j <= last[1]; ++j) {
            for (unsigned k = first[2]; k <= last[2]; ++k) {
                cell.centre = meshLower_m
                              + Vector_t<double, 3>(
                                      {(i + 0.5) * cellSize_m(0), (j + 0.5) * cellSize_m(1),
                                       (k + 0.5) * cellSize_m(2)});
                if (intersect(cell, box))
                    overlapped.push_back(static_cast<std::uint32_t>(flatIndex(i, j, k)));
            }
        }
    }
    return overlapped;
}
