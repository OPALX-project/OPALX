#ifndef OPALX_CYCLOTRON_SECTOR_FIELD_MAP_H
#define OPALX_CYCLOTRON_SECTOR_FIELD_MAP_H

#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

/**
 * @brief Shared PSI Ring median-plane magnetic map and first derivatives in SI units.
 *
 * The reader ports the layout of OPAL 2022.1 Cyclotron::getFieldFromFile_Ring:
 * radii in the file are millimetres, angles degrees, and magnetic values kilogauss.
 * Negative grid spacings denote the reciprocal of their absolute value.
 * File-supplied derivative blocks are consumed but recomputed from the base map
 * using OPAL's five-point Lagrange stencils (one-sided near either boundary).
 *
 * Both views have dimensions (nr, nt+1, 3). Channels 0, 1, 2 contain the stored
 * median-plane field F [T], dF/dr [T/m], and dF/dtheta [T/rad]. Column nt
 * copies column zero, including derivatives, to close the periodic angular seam.
 * The sign/coordinate conversion to the tracking frame belongs to CyclotronSector.
 *
 * @note Instances are shared by canonical filename and treated as read-only after
 * construction. The public Kokkos view types are writable; callers must not mutate
 * them. Changing a file while a cached instance remains alive does not reload it.
 * @see CyclotronSector
 */
class CyclotronSectorFieldMap {
public:
    using View = Kokkos::View<double***>;
    /** @brief Load or reuse a map; host-only, with Kokkos initialized.
     * @param filename Existing PSI map, resolved to a canonical filesystem path.
     * @return Shared map; its views remain alive while any sector owns the map.
     * @throws OpalException On malformed/truncated data or invalid grid dimensions.
     * @throws std::filesystem::filesystem_error If the filename cannot be resolved.
     */
    static std::shared_ptr<const CyclotronSectorFieldMap> read(const std::string& filename);
    double rmin, dr; ///< Minimum radius and positive radial spacing [m].
    double thetaMin, dtheta; ///< File angular origin and positive spacing [rad].
    int nr, nt; ///< Number of radial points and angular points excluding the seam.
    View data; ///< Device-accessible storage, populated once on the host.
    /// Host mirror for scalar reference tracking; same layout and values as data.
    decltype(Kokkos::create_mirror_view(std::declval<View>())) host;

    /** @brief Bilinear interpolation of one channel, callable on host or device.
     * @param grid View with dimensions (nr, nt+1, 3), accessible in the calling space.
     * @param channel 0=F, 1=dF/dr, 2=dF/dtheta; output units follow the channel.
     * @param u Fractional radial grid index (r-rmin)/dr, in [0,nr-1].
     * @param v Fractional angular grid index (theta-thetaMin)/dtheta, in [0,nt].
     * @param nr Number of radial grid points, at least two.
     * @param nt Number of angular intervals, at least one.
     * @return Interpolated channel value, without element scaling or sign conversion.
     * @pre The caller has checked bounds and supplied finite coordinates. This
     * routine does not wrap angles, reject negative indices, or validate channels.
     * At the upper endpoints it selects the last valid cell with weight one.
     */
    template <class V>
    KOKKOS_INLINE_FUNCTION static double interpolate(
            const V& grid, int channel, double u, double v, int nr, int nt) {
        int i = static_cast<int>(u), j = static_cast<int>(v);
        if (i >= nr - 1) i = nr - 2;
        if (j >= nt) j = nt - 1;
        const double a = u - i, b = v - j;
        return (1 - a) * (1 - b) * grid(i, j, channel) + a * (1 - b) * grid(i + 1, j, channel)
               + (1 - a) * b * grid(i, j + 1, channel) + a * b * grid(i + 1, j + 1, channel);
    }

private:
    explicit CyclotronSectorFieldMap(const std::string& filename);
};
#endif
