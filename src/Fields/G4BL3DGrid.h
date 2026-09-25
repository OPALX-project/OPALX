//
// Class G4BL3DGrid
//   Reader for G4beamline `grid` field maps.
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
#ifndef OPALX_FIELDMAPG4BL3DGRID_HH
#define OPALX_FIELDMAPG4BL3DGRID_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <string>

/**
 * @class G4BL3DGrid
 * @brief Reader for G4beamline `grid` field maps: a cartesian box of B values.
 *
 * @code
 * param current=1.
 * grid X0=-200.0 Y0=-120.0 Z0=-1380.0 nX=60 nY=25 nZ=553 dX=10.0 dY=10.0 dZ=5.0
 * data
 * <nX*nY*nZ rows of: x y z Bx By Bz>
 * @endcode
 * Positions are in mm and magnetic values in Tesla, interpolated trilinearly.
 *
 * @note Rows may also carry three more columns, x y z Bx By Bz Ex Ey Ez, which is how
 *       G4beamline writes a map that has an electric field as well. Those are in MV/m and
 *       are converted to V/m on load. The two fields are scaled independently, as in
 *       G4beamline: normB and current for the magnetic one, normE and gradient for the
 *       electric one. Storage for the electric field is only allocated once a non-zero
 *       value actually turns up, so the many nine-column maps whose Ex Ey Ez are all zero
 *       cost nothing extra.
 *
 * @note Rows are indexed by the x, y and z they carry, not by their position in the file.
 *       A missing, repeated or off-grid row is an error rather than a silent shift of the
 *       whole map. Real files do come out in a fixed order -- z fastest, then y, then x --
 *       but nothing here relies on that.
 *
 * @note Like G4BL2DMagnetoStatic and unlike every other OPALX ASCII reader, this one does
 *       @b not normalize. Values are absolute Tesla, so the element's scale is a plain
 *       multiplier and a scale of 1 reproduces a bare `fieldmap` placement.
 *
 * @note The header is parsed by the constructor and the data only by readMap(). These files
 *       are large -- the muE4 ASR61 dipole is 830k rows and 56 MB -- so unlike the 2D reader
 *       this one does not parse twice.
 */
class G4BL3DGrid : public Fieldmap {
public:
    /**
     * @brief Trilinear interpolation of any of the stored vector fields at @p R.
     *
     * The magnetic and the electric field sit on the same grid, so they share the weights
     * and this one function serves both. Adds into @p out. The caller is responsible for the
     * bounds check; a point on or past the last grid plane in any axis has no cell to
     * interpolate in.
     *
     * @note Templated on the view type, like FM2DMagnetoStatic::computeField(), so that the
     *       same code takes the device views in applyField() and the host views in
     *       getFieldstrength(). A fixed Kokkos::View<const double*> parameter lives in the
     *       default memory space, which is device memory on GPU builds, and a host view
     *       cannot be converted to it.
     */
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void interpolate(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& out, const ViewType& Vx,
            const ViewType& Vy, const ViewType& Vz, const double xbegin, const double ybegin,
            const double zbegin, const double hx, const double hy, const double hz, const int nx,
            const int ny, const int nz) {
        const double fx = (R(0) - xbegin) / hx;
        const double fy = (R(1) - ybegin) / hy;
        const double fz = (R(2) - zbegin) / hz;

        const int ix = static_cast<int>(Kokkos::floor(fx));
        const int iy = static_cast<int>(Kokkos::floor(fy));
        const int iz = static_cast<int>(Kokkos::floor(fz));

        if (ix < 0 || iy < 0 || iz < 0 || ix + 1 >= nx || iy + 1 >= ny || iz + 1 >= nz) {
            return;
        }

        const double wx = fx - ix;
        const double wy = fy - iy;
        const double wz = fz - iz;

        // z varies fastest, then y, then x -- the order the files themselves use.
        const size_t base = (static_cast<size_t>(ix) * ny + iy) * nz + iz;
        const size_t sx   = static_cast<size_t>(ny) * nz;
        const size_t sy   = static_cast<size_t>(nz);

        for (int corner = 0; corner < 8; ++corner) {
            const int cx       = (corner >> 2) & 1;
            const int cy       = (corner >> 1) & 1;
            const int cz       = corner & 1;
            const double w     = (cx ? wx : 1.0 - wx) * (cy ? wy : 1.0 - wy) * (cz ? wz : 1.0 - wz);
            const size_t index = base + cx * sx + cy * sy + cz;
            out(0) += w * Vx(index);
            out(1) += w * Vy(index);
            out(2) += w * Vz(index);
        }
    }

    virtual bool getFieldstrength(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) const override;

    /// @note Not implemented, throws.
    virtual bool getFieldDerivative(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B,
            const DiffDirection& dir) const override;

    virtual void getFieldDimensions(double& zBegin, double& zEnd) const override;

    virtual void getFieldDimensions(
            double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
            double& zFinal) const override;

    /// @note No-op. The grid format carries its own axes, so there is nothing to swap.
    virtual void swap() override;

    virtual void getInfo(Inform* msg) override;

    /// @note Not implemented, throws.
    virtual double getFrequency() const override;

    /// @note Not implemented, throws.
    virtual void setFrequency(double freq) override;

    /**
     * @brief Is @p r inside the box the map covers?
     *
     * The upper face in each axis is excluded, matching the interpolation, which needs a
     * full cell.
     */
    bool isInside(const Vector_t<double, 3>& r) const override {
        return r(0) >= xbegin_m && r(0) < xend_m && r(1) >= ybegin_m && r(1) < yend_m
               && r(2) >= zbegin_m && r(2) < zend_m;
    }

    void applyField(
            std::shared_ptr<ParticleContainer_t> pc, double scale = 1.0,
            double escale = 1.0) override;

    /// @brief Did the file actually carry a non-zero electric field?
    bool hasEField() const { return hasEField_m; }

private:
    /// @param aFilename Path to the map file.
    explicit G4BL3DGrid(std::string aFilename);
    ~G4BL3DGrid();

    void readMap() override;
    void freeMap() override;

    /// @brief Parse the header: the optional `param` line and the `grid` line.
    void readHeaderOnly();

    /// @brief Magnetic field components [T], indexed (ix * nY + iy) * nZ + iz.
    Kokkos::DualView<double*> FieldstrengthBx_m;
    Kokkos::DualView<double*> FieldstrengthBy_m;
    Kokkos::DualView<double*> FieldstrengthBz_m;

    /// @brief Electric field components [V/m], same indexing. Empty unless the file
    ///        turned out to carry a non-zero electric field.
    Kokkos::DualView<double*> FieldstrengthEx_m;
    Kokkos::DualView<double*> FieldstrengthEy_m;
    Kokkos::DualView<double*> FieldstrengthEz_m;

    /// @brief Box bounds [m].
    double xbegin_m;
    double xend_m;
    double ybegin_m;
    double yend_m;
    double zbegin_m;
    double zend_m;

    /// @brief Grid spacing [m].
    double hx_m;
    double hy_m;
    double hz_m;

    /// @brief Number of grid points.
    int num_gridpx_m;
    int num_gridpy_m;
    int num_gridpz_m;

    /// @brief Magnetic scale from the file's own `param` line, normB / current.
    double fieldScale_m;

    /// @brief Electric scale from the file's own `param` line, normE / gradient.
    double efieldScale_m;

    /// @brief Set once a non-zero electric field value is read, which is also when the
    ///        electric field storage is allocated.
    bool hasEField_m;

    friend class Fieldmap;
};

#endif  // OPALX_FIELDMAPG4BL3DGRID_HH
