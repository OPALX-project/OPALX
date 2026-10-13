#ifndef OPALX_FIELDMAP3DDYNAMIC_HH
#define OPALX_FIELDMAP3DDYNAMIC_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Cartesian RF amplitudes in the legacy OPAL 3DDynamic format.
 *
 * The header contains `3DDynamic [TRUE|FALSE]`, frequency [MHz], then one
 * `begin end intervals` record per axis [cm]. Each grid point contains
 * `Ex Ey Ez Hx Hy Hz` [MV/m, A/m], with z varying fastest, then y, then x.
 * Normalization defaults to TRUE and divides all fields by peak absolute Ez
 * interpolated at x = y = 0. Use FALSE for maps without a nonzero on-axis Ez.
 */
class FM3DDynamic : public Fieldmap {
public:
    /**
     * @brief Add the spatial field amplitudes at a position [m].
     * @param electricField Electric field [V/m].
     * @param magneticField Magnetic field [T].
     * @return true outside the map, false inside.
     */
    bool getFieldstrength(
            const Vector_t<double, 3>& position, Vector_t<double, 3>& electricField,
            Vector_t<double, 3>& magneticField) const override;

    bool getFieldDerivative(
            const Vector_t<double, 3>& position, Vector_t<double, 3>& electricField,
            Vector_t<double, 3>& magneticField, const DiffDirection& direction) const override;

    void getFieldDimensions(double& zBegin, double& zEnd) const override;

    void getFieldDimensions(
            double& xBegin, double& xEnd, double& yBegin, double& yEnd, double& zBegin,
            double& zEnd) const override;

    void swap() override;
    void getInfo(Inform* message) override;

    /// @return Angular frequency [rad/s].
    double getFrequency() const override;

    /// @param frequency Angular frequency [rad/s].
    void setFrequency(double frequency) override;

    /// @brief Host-side bounds check; lower faces are included and upper faces excluded.
    bool isInside(const Vector_t<double, 3>& position) const override {
        return position(0) >= xbegin_m && position(0) < xend_m && position(1) >= ybegin_m
               && position(1) < yend_m && position(2) >= zbegin_m && position(2) < zend_m;
    }

    /**
     * @brief Add trilinearly interpolated field amplitudes using host or device views.
     *
     * Positions, origins and spacings are in metres; fields are in V/m and T.
     * Grid counts are numbers of points. Storage has z varying fastest, then y, then x.
     * Positions outside the grid or on its upper faces leave the output unchanged.
     * RF amplitude and phase factors are applied by the caller.
     */
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void computeField(
            const Vector_t<double, 3>& position, Vector_t<double, 3>& electricField,
            Vector_t<double, 3>& magneticField, const ViewType& fieldEx, const ViewType& fieldEy,
            const ViewType& fieldEz, const ViewType& fieldBx, const ViewType& fieldBy,
            const ViewType& fieldBz, double spacingX, double spacingY, double spacingZ,
            double xBegin, double yBegin, double zBegin, double xEnd, double yEnd, double zEnd,
            int numGridPointsX, int numGridPointsY, int numGridPointsZ) {
        if (numGridPointsX < 2 || numGridPointsY < 2 || numGridPointsZ < 2
            || !(spacingX > 0.0 && spacingY > 0.0 && spacingZ > 0.0) || !Kokkos::isfinite(spacingX)
            || !Kokkos::isfinite(spacingY) || !Kokkos::isfinite(spacingZ)) {
            return;
        }

        if (!(position(0) >= xBegin && position(0) < xEnd && position(1) >= yBegin
              && position(1) < yEnd && position(2) >= zBegin && position(2) < zEnd)) {
            return;
        }

        const double gridX =
                Kokkos::fmin((position(0) - xBegin) / spacingX, double(numGridPointsX - 1));
        const double gridY =
                Kokkos::fmin((position(1) - yBegin) / spacingY, double(numGridPointsY - 1));
        const double gridZ =
                Kokkos::fmin((position(2) - zBegin) / spacingZ, double(numGridPointsZ - 1));

        const int indexX = Kokkos::min(static_cast<int>(gridX), numGridPointsX - 2);
        const int indexY = Kokkos::min(static_cast<int>(gridY), numGridPointsY - 2);
        const int indexZ = Kokkos::min(static_cast<int>(gridZ), numGridPointsZ - 2);

        const double fractionX = gridX - indexX;
        const double fractionY = gridY - indexY;
        const double fractionZ = gridZ - indexZ;

        const std::size_t strideX = static_cast<std::size_t>(numGridPointsY) * numGridPointsZ;
        const std::size_t strideY = static_cast<std::size_t>(numGridPointsZ);
        const std::size_t baseIndex =
                static_cast<std::size_t>(indexX) * strideX + indexY * strideY + indexZ;

        for (int corner = 0; corner < 8; ++corner) {
            const int offsetX   = (corner >> 2) & 1;
            const int offsetY   = (corner >> 1) & 1;
            const int offsetZ   = corner & 1;
            const double weight = (offsetX ? fractionX : 1.0 - fractionX)
                                  * (offsetY ? fractionY : 1.0 - fractionY)
                                  * (offsetZ ? fractionZ : 1.0 - fractionZ);
            const std::size_t index = baseIndex + offsetX * strideX + offsetY * strideY + offsetZ;

            electricField(0) += weight * fieldEx(index);
            electricField(1) += weight * fieldEy(index);
            electricField(2) += weight * fieldEz(index);
            magneticField(0) += weight * fieldBx(index);
            magneticField(1) += weight * fieldBy(index);
            magneticField(2) += weight * fieldBz(index);
        }
    }

    void applyField(std::shared_ptr<ParticleContainer_t> particles, double, double) override;

    /// @brief Apply separate electric and magnetic RF factors within the active z interval.
    void applyRFField(
            std::shared_ptr<ParticleContainer_t> particles, double electricScale,
            double magneticScale, double startField, double endField);

    void getOnaxisEz(std::vector<std::pair<double, double>>& field) override;

private:
    explicit FM3DDynamic(const std::string& filename);
    ~FM3DDynamic() override;

    void readMap() override;
    void freeMap() override;
    void requireLoaded() const;
    std::vector<double> sampleOnaxisEz(const Kokkos::DualView<double*>::t_host& fieldEz) const;

    Kokkos::DualView<double*> FieldstrengthEx_m;
    Kokkos::DualView<double*> FieldstrengthEy_m;
    Kokkos::DualView<double*> FieldstrengthEz_m;
    Kokkos::DualView<double*> FieldstrengthBx_m;
    Kokkos::DualView<double*> FieldstrengthBy_m;
    Kokkos::DualView<double*> FieldstrengthBz_m;

    double frequency_m = 0.0;

    double xbegin_m = 0.0;
    double xend_m   = 0.0;
    double ybegin_m = 0.0;
    double yend_m   = 0.0;
    double zbegin_m = 0.0;
    double zend_m   = 0.0;

    double hx_m = 0.0;
    double hy_m = 0.0;
    double hz_m = 0.0;

    int num_gridpx_m = 0;
    int num_gridpy_m = 0;
    int num_gridpz_m = 0;

    friend class Fieldmap;
};

#endif
