/**
 * @file FFT2D5FieldStorage.h
 * @brief Owns persistent 3D staging fields and the FFT2D5 slice solver array.
 */

#ifndef OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H
#define OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H

#include "Manager/BaseManager.h"
#include "Manager/datatypes.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace opalx::spacecharge {

    /**
     * @brief Persistent fields and per-slice Poisson solvers for one FFT2D5 run.
     *
     * The 3D fields stage deposition and gathering. Every transverse slice owns separate 2D IPPL
     * fields because the Field API cannot use a Kokkos subview as its storage.
     */
    class FFT2D5FieldStorage final {
    public:
        using Vector2      = Vector_t<double, 2>;
        using Vector3      = Vector_t<double, 3>;
        using Mesh2        = Mesh_t<2>;
        using Mesh3        = Mesh_t<3>;
        using Layout2      = FieldLayout_t<2>;
        using Layout3      = FieldLayout_t<3>;
        using ScalarField2 = Field_t<2>;
        using ScalarField3 = Field_t<3>;
        using VectorField2 = VField_t<double, 2>;
        using VectorField3 = VField_t<double, 3>;
        using OpenSolver2  = ippl::FFTOpenPoissonSolver<VectorField2, ScalarField2>;

        struct Slice final {
            std::unique_ptr<VectorField2> electricField;
            std::unique_ptr<ScalarField2> chargeDensity;
            std::unique_ptr<OpenSolver2> solver;
        };

        FFT2D5FieldStorage(const FFT2D5Config& config, double pathLength);

        FFT2D5FieldStorage(const FFT2D5FieldStorage&)            = delete;
        FFT2D5FieldStorage& operator=(const FFT2D5FieldStorage&) = delete;
        FFT2D5FieldStorage(FFT2D5FieldStorage&&)                 = delete;
        FFT2D5FieldStorage& operator=(FFT2D5FieldStorage&&)      = delete;

        [[nodiscard]] const Vector3& spacing() const { return spacing_m; }
        [[nodiscard]] const Vector3& origin() const { return origin_m; }
        [[nodiscard]] const Vector3& size() const { return size_m; }
        [[nodiscard]] const std::array<std::size_t, 3>& meshSize() const { return meshSize_m; }
        [[nodiscard]] Mesh2& sliceMesh() { return sliceMesh_m; }
        [[nodiscard]] Layout2& sliceLayout() { return sliceLayout_m; }
        [[nodiscard]] Layout3& layout() { return layout_m; }
        [[nodiscard]] ScalarField3& chargeDensity() { return chargeDensity_m; }
        [[nodiscard]] VectorField3& electricField() { return electricField_m; }
        [[nodiscard]] std::vector<Slice>& slices() { return slices_m; }

    private:
        std::array<std::size_t, 3> meshSize_m;
        Vector3 spacing_m;
        Vector3 origin_m;
        Vector3 size_m;
        ippl::NDIndex<3> domain_m;
        Mesh3 mesh_m;
        Layout3 layout_m;
        ScalarField3 chargeDensity_m;
        VectorField3 electricField_m;
        ippl::NDIndex<2> sliceDomain_m;
        Mesh2 sliceMesh_m;
        Layout2 sliceLayout_m;
        ippl::ParameterList solverParameters_m;
        std::vector<Slice> slices_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FFT2D5_FIELD_STORAGE_H
