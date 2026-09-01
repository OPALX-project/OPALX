/**
 * @file ReferencePath.h
 * @brief Owns the validated host and device representation of an FFT2D5 design path.
 */

#ifndef OPALX_SPACE_CHARGE_PIC2D5_REFERENCE_PATH_H
#define OPALX_SPACE_CHARGE_PIC2D5_REFERENCE_PATH_H

#include "Manager/BaseManager.h"
#include "Manager/datatypes.h"

#include <Kokkos_Core.hpp>

#include <string>

namespace opalx::spacecharge {

    /** @brief Immutable reference path used for Cartesian to Frenet mapping. */
    class ReferencePath final {
    public:
        using Vector = Vector_t<double, 3>;
        using View   = Kokkos::View<Vector*>;

        [[nodiscard]] static ReferencePath load(const std::string& filename);

        ReferencePath(const ReferencePath&)            = delete;
        ReferencePath& operator=(const ReferencePath&) = delete;
        ReferencePath(ReferencePath&&)                 = default;
        ReferencePath& operator=(ReferencePath&&)      = default;

        [[nodiscard]] const View& deviceView() const { return deviceView_m; }
        [[nodiscard]] double length() const { return length_m; }

    private:
        ReferencePath(View deviceView, double length)
            : deviceView_m(std::move(deviceView)), length_m(length) {}

        View deviceView_m;
        double length_m = 0.0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC2D5_REFERENCE_PATH_H
