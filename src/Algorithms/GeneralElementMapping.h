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
#ifndef OPALX_GENERAL_ELEMENT_MAPPING_H
#define OPALX_GENERAL_ELEMENT_MAPPING_H

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/IndexMap.h"
#include "Utilities/ElementList.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

/**
 * @brief Uniform lab-frame grid that maps a bunch box to the elements it can reach.
 *
 * Each element contributes oriented boxes: ElementBase::getLocalSupportBoxes() for its reach
 * (ElementBase::getReach(), otherwise the global default), placed in the lab by the inverse
 * of misalignment * placement, the frame in which the tracker applies the element. A bunch is
 * an oriented box as well, given by its bounds in the bunch frame and the bunch-to-lab
 * transform. An element is active iff the bunch box intersects one of its boxes (separating-axis
 * test). An active element acts on the whole bunch with its own kernels, so a bunch that
 * reaches into the box of a neighbouring branch's element receives that element entirely.
 *
 * The grid covers the lab bounding box of all element boxes with Nx x Ny x Nz cells. Each
 * cell lists the elements whose box intersects it, as compressed rows in the order of the
 * element list passed to the constructor. A query visits only the cells under the bunch's lab
 * bounding box, so its cost depends on the number of elements per cell rather than on the
 * size of the lattice.
 *
 * The grid is built once in a host Kokkos loop over the elements. Results do not depend on
 * thread scheduling, and the same element order on all ranks gives the same grid everywhere.
 * Queries run on the host.
 */
class GeneralElementMapping {
public:
    using ElementSet = IndexMap::value_t;

    /// Upper limit on Nx * Ny * Nz; the row offsets alone take 8 bytes per cell.
    static constexpr std::size_t maximumCells = std::size_t(1) << 24;

    /// A reference-orbit point inside an element's support (isInside) but outside its box.
    struct SupportOutsideBox {
        std::shared_ptr<ElementBase> element;
        /// Largest transverse offset [m] of the point from the element axis, in the arc chart
        /// for sector bends. A reach of at least this includes the point unless the element's
        /// own extent (field map, bend orbit, cavity width) is what the box misses.
        double offset;
        /// Reach [m] the element's box was built with.
        double reach;
    };

    /**
     * @param elements Placed runtime elements. Placement, misalignment and field-map extents
     *                 must be final.
     * @param cells Grid cells along lab x, y and z, each at least 1.
     * @param reach Global reach [m] for elements without their own REACH.
     * @throws OpalException If a cell count is zero or the reach is not positive.
     */
    GeneralElementMapping(
            const ElementList& elements, const std::array<unsigned, 3>& cells, double reach);

    /**
     * @brief Elements whose box intersects the bunch box.
     * @param lower Lower corner of the bunch box in the bunch frame [m].
     * @param upper Upper corner of the bunch box in the bunch frame [m].
     * @param toLab Bunch-to-lab transform: transformTo() maps a bunch-frame point to the lab.
     */
    ElementSet query(
            const Vector_t<double, 3>& lower, const Vector_t<double, 3>& upper,
            const CoordinateSystemTrafo& toLab) const;

    /**
     * @brief Check every element's reach against one point of a reference orbit.
     *
     * Throws when the point lies inside an element's box and body window but outside its
     * aperture while the element deletes particles there: a particle on the reference orbit
     * would be deleted. Either the element belongs to another branch and its reach is too
     * large, or the orbit hits the aperture.
     *
     * Returns the elements whose support contains the point but whose box does not, which a
     * bunch on this orbit does not activate. An element without an aperture in another branch
     * is expected here, so this is reported rather than thrown.
     * @param position Reference position in the lab [m].
     * @param pathLength Reference path length [m], used in the error message.
     * @throws OpalException As described above.
     */
    std::vector<SupportOutsideBox> checkReferencePoint(
            const Vector_t<double, 3>& position, double pathLength) const;

    /**
     * @brief Whether @p position lies within the element's transverse reach.
     *
     * True if the point lies inside one of the element's boxes when their longitudinal extent
     * is ignored. A bunch activates an element whose reach misses the reference only through
     * its own width, as for an element of a neighbouring branch. Always true for
     * CYCLOTRONSECTOR and SCALINGFFAMAGNET, whose boxes span radius and angle.
     * @param element An element passed to the constructor.
     * @param position Lab position [m].
     * @throws OpalException If @p element is not mapped.
     */
    bool coversTransversely(const ElementBase& element, const Vector_t<double, 3>& position) const;

    /// Number of mapped elements.
    std::size_t size() const { return entries_m.size(); }

private:
    struct Entry {
        std::shared_ptr<ElementBase> element;
        CoordinateSystemTrafo labToLocal;    ///< Misalignment after placement.
        double reach;                        ///< Reach used for the boxes [m].
        std::vector<SupportBox> localBoxes;  ///< ElementBase::getLocalSupportBoxes().
        std::vector<SupportBox> boxes;       ///< The same boxes in the lab.

        // Loss window of the element's markOutsideAperture(), for checkReferencePoint().
        bool deletes;
        double arcCurvature;  ///< Nonzero if the window is in the sector-bend arc chart.
        double bodyBegin;
        double bodyEnd;
        ApertureType apertureType;
        double apertureX;
        double apertureY;
    };

    static SupportBox placeInLab(const SupportBox& local, const CoordinateSystemTrafo& labToLocal);
    static Vector_t<double, 3> labHalfExtent(const SupportBox& box);
    static bool intersect(const SupportBox& a, const SupportBox& b);
    static bool contains(const std::vector<SupportBox>& boxes, const Vector_t<double, 3>& point);

    /// Cell index ranges under a lab bounding box; false if the box misses the mesh.
    bool cellRange(
            const Vector_t<double, 3>& lower, const Vector_t<double, 3>& upper,
            std::array<unsigned, 3>& first, std::array<unsigned, 3>& last) const;
    std::size_t flatIndex(unsigned i, unsigned j, unsigned k) const;
    std::vector<std::uint32_t> overlappedCells(const SupportBox& box) const;

    std::vector<Entry> entries_m;
    std::unordered_map<const ElementBase*, std::size_t> entryOf_m;
    std::array<unsigned, 3> cells_m;
    Vector_t<double, 3> meshLower_m;
    Vector_t<double, 3> meshUpper_m;
    Vector_t<double, 3> cellSize_m;

    /// Compressed rows: the elements of cell c are cellElements_m[cellOffsets_m[c] ..
    /// cellOffsets_m[c + 1]), in increasing entry index.
    std::vector<std::size_t> cellOffsets_m;
    std::vector<std::uint32_t> cellElements_m;
};

#endif  // OPALX_GENERAL_ELEMENT_MAPPING_H
