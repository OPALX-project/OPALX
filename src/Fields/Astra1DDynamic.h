#ifndef CLASSIC_AstraFIELDMAP1DDYNAMIC_HH
#define CLASSIC_AstraFIELDMAP1DDYNAMIC_HH

#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

class Astra1DDynamic: public Fieldmap {

public:
    /**
     * @brief Get the field strength at a given point.
     * 
     * @param R Position [m] relative to the field map origin.
     * @param E Output Electric field [MV/m].
     * @param B Output Magnetic field [T].
     * @return true if R is outside of the field map, false otherwise.
     */
    virtual bool getFieldstrength(
        const Vector_t<double, 3>& R, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) const override;

    /**
     * @brief Get the field derivative with respect to a direction.
     * 
     * @param R Position [m].
     * @param E Output derivative of Electric field.
     * @param B Output derivative of Magnetic field.
     * @param dir Direction of derivative (DX, DY, DZ).
     * @return true if R is outside, false otherwise.
     */
    virtual bool getFieldDerivative(
        const Vector_t<double, 3>& R, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B, 
        const DiffDirection& dir) const override;

    /**
     * @brief Get the longitudinal dimensions of the field.
     * @param zBegin Output start of field [m].
     * @param zEnd Output end of field [m].
     */
    virtual void getFieldDimensions(double& zBegin, double& zEnd) const override;

    /**
     * @brief Get the full 3D bounding box of the field.
     * 
     * @param xIni Output minimum x [m].
     * @param xFinal Output maximum x [m].
     * @param yIni Output minimum y [m].
     * @param yFinal Output maximum y [m].
     * @param zIni Output minimum z [m].
     * @param zFinal Output maximum z [m].
     */
    virtual void getFieldDimensions(
        double& xIni, double& xFinal, 
        double& yIni, double& yFinal, 
        double& zIni, double& zFinal) const override;

    /// @brief Swap coordinates
    virtual void swap() override;

    /// @brief Print info about the field map.
    virtual void getInfo(Inform *msg) override;

    /**
     * @brief Get the frequency.
     * @return Frequency [MHz].
     * @note Not implemented yet
     */
    virtual double getFrequency() const override;

    /**
     * @brief Set the frequency.
     * @param freq Frequency [MHz].
     * @note Not implemented yet
     */
    virtual void setFrequency(double freq) override;

    /**
     * @brief Checks if the given coordinate is inside the volume covered by the
     * fieldmap
     * @param r Coordinate
     * @note This cannot be called inside a GPU kernel (implicit capture of the
     * 'this' pointer not allowed on device)
     */
    bool isInside(const Vector_t<double, 3> &r) const override {
        return r(2) >= zbegin_m && r(2) < zend_m;
        // && sqrt(r(0) * r(0) + r(1) * r(1)) < rend_m;
    }

    template <class ViewType>
    KOKKOS_INLINE_FUNCTION
    static void computeField(
        const Vector_t<double, 3>& R,
        Vector_t<double, 3>& E,
        Vector_t<double, 3>& B,
        const ViewType& FourCoefs,
        double zbegin,
        double length,
        double xlrep,
        int accuracy)
    {
        const double RR2 = R(0) * R(0) + R(1) * R(1);

        const double two_pi = Physics::two_pi;
        const double pi     = Physics::pi;
        const double dk     = Physics::two_pi / length;
        const double kz     = two_pi * (R(2) - zbegin) / length + pi;
        
        double ez    = FourCoefs(0);
        double ezp   = 0.0;
        double ezpp  = 0.0;
        double ezppp = 0.0;

        int n = 1;
        for (int l = 1; l < accuracy; ++l, n += 2) {
            const double base = dk * l;

            const double coskzl = Kokkos::cos(kz * l);
            const double sinkzl = Kokkos::sin(kz * l);

            ez    += FourCoefs(n)       * coskzl - FourCoefs(n + 1) * sinkzl;
            ezp   += base               * (-FourCoefs(n) * sinkzl - FourCoefs(n + 1) * coskzl);
            ezpp  += base * base        * (-FourCoefs(n) * coskzl + FourCoefs(n + 1) * sinkzl);
            ezppp += base * base * base * ( FourCoefs(n) * sinkzl + FourCoefs(n + 1) * coskzl);
        }

        const double f  = -(ezpp  + ez  * xlrep * xlrep) / 16.0;
        const double fp = -(ezppp + ezp * xlrep * xlrep) / 16.0;

        const double EfieldR = -(ezp / 2.0 + fp * RR2);
        const double BfieldT = (ez / 2.0 + f * RR2) * xlrep / Physics::c;

        E(0) += EfieldR * R(0);
        E(1) += EfieldR * R(1);
        E(2) += ez + 4.0 * f * RR2;

        B(0) -= BfieldT * R(1);
        B(1) += BfieldT * R(0);
    }

    template <class ViewType>
    KOKKOS_INLINE_FUNCTION
    static void computeTravelingWaveField(
        const Vector_t<double, 3>& R,
        Vector_t<double, 3>& E,
        Vector_t<double, 3>& B,
        const ViewType& FourCoefs,
        double zbegin,
        double zend,
        double length,
        double xlrep,
        int accuracy,
        double entryElectricScale,
        double entryMagneticScale,
        double core1ElectricScale,
        double core1MagneticScale,
        double core2ElectricScale,
        double core2MagneticScale,
        double exitElectricScale,
        double exitMagneticScale,
        double startCoreField,
        double startExitField,
        double mappedStartExitField,
        double periodLength,
        double cellLength,
        double elementLength)
    {
        if (R(2) < -0.5 * periodLength || R(2) + 0.5 * periodLength >= elementLength) {
            return;
        }

        Vector_t<double, 3> tmpR({R(0), R(1), R(2) + 0.5 * periodLength});
        Vector_t<double, 3> tmpE(0.0), tmpB(0.0);

        if (tmpR(2) < startCoreField) {
            if (!(tmpR(2) >= zbegin && tmpR(2) < zend)) {
                return;
            }

            computeField(tmpR, tmpE, tmpB, FourCoefs, zbegin, length, xlrep, accuracy);
            E += entryElectricScale * tmpE;
            B += entryMagneticScale * tmpB;
        }
        else if (tmpR(2) < startExitField) {
            tmpR(2) -= startCoreField;
            const double z = tmpR(2);

            tmpR(2) = tmpR(2) - periodLength * Kokkos::floor(tmpR(2) / periodLength);
            tmpR(2) += startCoreField;

            if (!(tmpR(2) >= zbegin && tmpR(2) < zend)) {
                return;
            }

            computeField(tmpR, tmpE, tmpB, FourCoefs, zbegin, length, xlrep, accuracy);
            E += core1ElectricScale * tmpE;
            B += core1MagneticScale * tmpB;

            tmpE = 0.0;
            tmpB = 0.0;

            tmpR(2) = z + cellLength;
            tmpR(2) = tmpR(2) - periodLength * Kokkos::floor(tmpR(2) / periodLength);
            tmpR(2) += startCoreField;

            if (!(tmpR(2) >= zbegin && tmpR(2) < zend)) {
                return;
            }

            computeField(tmpR, tmpE, tmpB, FourCoefs, zbegin, length, xlrep, accuracy);
            E += core2ElectricScale * tmpE;
            B += core2MagneticScale * tmpB;
        }
        else {
            tmpR(2) -= mappedStartExitField;

            if (!(tmpR(2) >= zbegin && tmpR(2) < zend)) {
                return;
            }

            computeField(tmpR, tmpE, tmpB, FourCoefs, zbegin, length, xlrep, accuracy);
            E += exitElectricScale * tmpE;
            B += exitMagneticScale * tmpB;
        }
    }

    /**
     * @brief Apply the FM to all the particles
     * 
     * @param pc Particle container
     */
    void applyField(std::shared_ptr<ParticleContainer_t> pc, double) override;

    /**
     * @brief Apply the traveling-wave RF field map to all particles.
     *
     * This is the device-safe path used by `TravelingWave`. It avoids calling the
     * virtual `Fieldmap::getFieldstrength()` interface inside a Kokkos kernel and
     * reproduces the traveling-wave cavity region logic directly on device.
     *
     * The cavity is split into three longitudinal regions:
     * - entry region: direct field-map application with entry phase/scale
     * - core region: superposition of two periodically shifted field contributions
     * - exit region: remapped field-map application with exit phase/scale
     *
     * @param pc Particle container.
     *
     * @param entryElectricScale Scale factor applied to E in the entry region.
     * @param entryMagneticScale Scale factor applied to B in the entry region.
     *
     * @param core1ElectricScale Scale factor applied to the first core contribution in E.
     * @param core1MagneticScale Scale factor applied to the first core contribution in B.
     *
     * @param core2ElectricScale Scale factor applied to the second core contribution in E.
     * @param core2MagneticScale Scale factor applied to the second core contribution in B.
     *
     * @param exitElectricScale Scale factor applied to E in the exit region.
     * @param exitMagneticScale Scale factor applied to B in the exit region.
     *
     * @param startField Longitudinal beginning of the full traveling-wave structure.
     * @param startCoreField Longitudinal beginning of the periodic core region.
     * @param startExitField Longitudinal beginning of the exit region.
     * @param mappedStartExitField Shift used to remap the exit region into field-map coordinates.
     * @param periodLength Period length of the base field-map cell.
     * @param cellLength Cell-to-cell longitudinal shift in the core region.
     * @param elementLength Total longitudinal length of the traveling-wave element.
     */
    void applyTravelingWave(
        std::shared_ptr<ParticleContainer_t> pc,
        double entryElectricScale,
        double entryMagneticScale,
        double core1ElectricScale,
        double core1MagneticScale,
        double core2ElectricScale,
        double core2MagneticScale,
        double exitElectricScale,
        double exitMagneticScale,
        double startField,
        double startCoreField,
        double startExitField,
        double mappedStartExitField,
        double periodLength,
        double cellLength,
        double elementLength);

    virtual void getOnaxisEz(std::vector<std::pair<double, double>> & F) override;

private:
    Astra1DDynamic(const std::string& filename);
    ~Astra1DDynamic();

    void readMap() override;
    void freeMap() override;

    /**
     * @brief Fourier coefficients of Ez(z) (device-accessible)
     * Stored as:
     * [a0, a1, b1, a2, b2, ..., aN, bN]
     * and used in the reconstruction
     *   Ez(z) = a0 + sum_l [ a_l cos(l*kz) - b_l sin(l*kz) ].
     */
    Kokkos::DualView<double*> FourCoefs_m;

    std::vector<double> zRaw_m;
    std::vector<double> ezRaw_m;

    double frequency_m;

    /// @brief Z Bounds relative to element edge
    double zbegin_m;    /**< Start position of field map [m] */
    double zend_m;      /**< End position of field map [m] */

    /// @brief Wave number (omega / c)
    double xlrep_m;

    /// @brief Effective periodic length of the field map [m]
    double length_m;

    /// @brief Number of Fourier modes used
    int accuracy_m;

    /// @brief Number of grid points in z-direction (input sampling)
    int num_gridpz_m;

    /// @brief Allow Fieldmap factory access
    friend class Fieldmap;
};

#endif