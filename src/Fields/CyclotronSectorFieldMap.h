#ifndef OPALX_CYCLOTRON_SECTOR_FIELD_MAP_H
#define OPALX_CYCLOTRON_SECTOR_FIELD_MAP_H

#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

/** PSI median-plane map in SI units. The final angular column is a periodic seam.
 * Channels are B, dB/dr and dB/dtheta; derivatives use the legacy five-point stencil.
 * Immutable maps are shared by canonical filename. Host and device interpolation agree.
 */
class CyclotronSectorFieldMap {
public:
    using View = Kokkos::View<double***>;
    static std::shared_ptr<const CyclotronSectorFieldMap> read(const std::string& filename);
    double rmin, dr, thetaMin, dtheta;
    int nr, nt;
    View data;
    decltype(Kokkos::create_mirror_view(std::declval<View>())) host;

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
