#ifndef OPALX_FIELDMAP2DDYNAMIC_HH
#define OPALX_FIELDMAP2DDYNAMIC_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

class FM2DDynamic: public Fieldmap {

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
        const Vector_t<double, 3> &R, 
        Vector_t<double, 3> &E, 
        Vector_t<double, 3> &B) const override;
    // virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;

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
        const Vector_t<double, 3> &R, 
        Vector_t<double, 3> &E, 
        Vector_t<double, 3> &B, 
        const DiffDirection &dir) const override;
    // virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;

    /**
     * @brief Get the longitudinal dimensions of the field.
     * @param zBegin Output start of field [m].
     * @param zEnd Output end of field [m].
     */
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const override;
    // virtual void getFieldDimensions(double &zBegin, double &zEnd) const;

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
        double &xIni, double &xFinal, 
        double &yIni, double &yFinal, 
        double &zIni, double &zFinal) const override;

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
        return r(2) >= zbegin_m && r(2) < zend_m && 
        sqrt(r(0) * r(0) + r(1) * r(1)) < rend_m;
    }

    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void computeField(
        const Vector_t<double, 3>& R,
        Vector_t<double, 3>& E,
        Vector_t<double, 3>& B,
        const ViewType& Ez,
        const ViewType& Er,
        const ViewType& Bt,
        double hr_m, 
        double hz_m, 
        double zbegin_m, 
        int num_gridpr_m, 
        int num_gridpz_m) {

        const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

        const int indexr    = (int)floor(RR / hr_m);
        const double leverr = (RR / hr_m) - indexr;

        const int indexz    = (int)floor((R(2) - zbegin_m) / hz_m);
        const double leverz = (R(2) - zbegin_m) / hz_m - indexz;

        if ((indexz < 0) || (indexz + 2 > num_gridpz_m))
            return;

        if(indexr + 2 > num_gridpr_m)
            return;

        const int index1 = indexz + indexr * num_gridpz_m;
        const int index2 = index1 + num_gridpz_m;

        const double EfieldR = (1.0 - leverz)  * (1.0 - leverr) * Er(index1)
                            + leverz         * (1.0 - leverr) * Er(index1 + 1)
                            + (1.0 - leverz) * leverr         * Er(index2)
                            + leverz         * leverr         * Er(index2 + 1);

        const double EfieldZ = (1.0 - leverz)  * (1.0 - leverr) * Ez(index1)
                            + leverz         * (1.0 - leverr) * Ez(index1 + 1)
                            + (1.0 - leverz) * leverr         * Ez(index2)
                            + leverz         * leverr         * Ez(index2 + 1);

        const double BfieldT = (1.0 - leverz)  * (1.0 - leverr) * Bt(index1)
                            + leverz         * (1.0 - leverr) * Bt(index1 + 1)
                            + (1.0 - leverz) * leverr         * Bt(index2)
                            + leverz         * leverr         * Bt(index2 + 1);

        if(RR > 1e-14) {
            E(0) += EfieldR * R(0) / RR;
            E(1) += EfieldR * R(1) / RR;
            B(0) -= BfieldT * R(1) / RR;
            B(1) += BfieldT * R(0) / RR;
        }
        E(2) += EfieldZ;
        return;
    }

    /**
     * @brief Apply the FM to all the particles
     * 
     * @param pc Particle container
     */
    void applyField(std::shared_ptr<ParticleContainer_t> pc, double) override;

    /**
     * @brief Apply the RF-scaled dynamic field map to all particles.
     *
     * This is the device-safe path used by `RFCavity`. It avoids calling the
     * virtual `Fieldmap::getFieldstrength()` interface inside a Kokkos kernel.
     *
     * @param pc Particle container.
     * @param electricScale Scale factor applied to the electric field.
     * @param magneticScale Scale factor applied to the magnetic field.
     * @param startField Begin of the active cavity region.
     * @param endField End of the active cavity region.
     */
    void applyRFField(
        std::shared_ptr<ParticleContainer_t> pc,
        double electricScale,
        double magneticScale,
        double startField,
        double endField);
 
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F) override;

private:
    FM2DDynamic(const std::string& filename);
    ~FM2DDynamic();

    void readMap() override;
    void freeMap() override;

    /// @brief Fieldstrengths
    Kokkos::DualView<double*> FieldstrengthEz_m;    /**< 2D array with Ez, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
    Kokkos::DualView<double*> FieldstrengthEr_m;    /**< 2D array with Er, read in like Ez*/
    Kokkos::DualView<double*> FieldstrengthBt_m;    /**< 2D array with Er, read in like Ez*/

    double frequency_m;

    /// @brief Radius Bounds
    double rbegin_m;
    double rend_m;

    /// @brief Z Bounds relative to element edge
    double zbegin_m;
    double zend_m;

    /// @brief Grid
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    double hr_m;                   /**< length between points in grid, r-direction, m*/
    int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

#endif
