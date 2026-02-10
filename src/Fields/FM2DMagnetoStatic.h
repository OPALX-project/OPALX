#ifndef CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH
#define CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

class FM2DMagnetoStatic: public Fieldmap {

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

    /**
     * @brief Get the longitudinal dimensions of the field.
     * @param zBegin Output start of field [m].
     * @param zEnd Output end of field [m].
     */
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const override;

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
     */
    virtual double getFrequency() const override;

    /**
     * @brief Set the frequency.
     * @param freq Frequency [MHz].
     */
    virtual void setFrequency(double freq) override;

    /**
     * @brief Checks if the given coordinate is inside the volume covered by the
     * fieldmap
     * @param r Coordinate
     * @note This cannot be called inside a GPU kernel (implicit capture of the
     * 'this' pointer not allowed on device)
     */
    bool isInside(const Vector_t<double, 3> &r)const{
        return r(2) >= zbegin_m && r(2) < zend_m && 
        sqrt(r(0)*r(0) + r(1)*r(1)) < rend_m;
    }

    /**
     * @brief Computes the magnetic field B at the position R by interpolating
     * from the fieldmap specified by Bz, Br, hr, hz, zbegin, num_gridpr, and
     * num_gridpz.
     *
     * @param R Position [m] relative to the element edge
     * @param B Output magnetic field [T]
     * @param Bz Longitudinal value of the fieldmap [T]
     * @param Br Radial value of the fieldmap [T]
     * @param hr Radial grid spacing [m]
     * @param hz Longitudinal grid spacing [m]
     * @param zbegin Start of the fieldmap relative to the element edge [m]
     * @param num_gridpr Number of radial grid points in the fieldmap
     * @param num_gridpz Number of longitudinal grid points in the fieldmap
     * 
     * @note We need to pass the fieldmap data as arguments to this function 
     * because the 'this' pointer cannot be captured in a GPU kernel, 
     * and thus we cannot access member variables directly.
     */
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static bool computeField(
        const Vector_t<double, 3>& R, 
        Vector_t<double, 3>& B,
        const ViewType& Bz, 
        const ViewType& Br,
        double hr, 
        double hz, 
        double zbegin, 
        int num_gridpr, 
        int num_gridpz) {
        
        double RR = sqrt(R(0) * R(0) + R(1) * R(1));

        int indexr    = (int)floor(RR / hr);
        double leverr = (RR / hr) - indexr;

        int indexz    = (int)floor((R(2) - zbegin) / hz);
        double leverz = (R(2) - zbegin) / hz - indexz;

        if ((indexz < 0) || (indexz + 2 > num_gridpz))
            return false;

        if (indexr + 2 > num_gridpr)
            return true;

        int index1     = indexz + indexr * num_gridpz;
        int index2     = index1 + num_gridpz;
        double BfieldR = (1.0 - leverz) * (1.0 - leverr) * Br(index1)
                               + leverz * (1.0 - leverr) * Br(index1 + 1)
                               + (1.0 - leverz) * leverr * Br(index2)
                               + leverz * leverr * Br(index2 + 1);

        double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * Bz(index1)
                               + leverz * (1.0 - leverr) * Bz(index1 + 1)
                               + (1.0 - leverz) * leverr * Bz(index2)
                               + leverz * leverr * Bz(index2 + 1);

        if (RR != 0) {
            B(0) += BfieldR * R(0) / RR;
            B(1) += BfieldR * R(1) / RR;
        }
        B(2) += BfieldZ;
        return false;
    }

    /**
     * @brief Compute the derivative of the magnetic field B at the position R.
     * 
     * @param R Position [m].
     * @param B Output derivative of the magnetic field.
     * @param Bz Longitudinal value of the fieldmap [T].
     * @param Br Radial value of the fieldmap [T].
     * @param hr Radial grid spacing [m].
     * @param hz Longitudinal grid spacing [m].
     * @param zbegin Start of the fieldmap relative to the element edge [m].
     * @param num_gridpr Number of radial grid points in the fieldmap.
     * @param num_gridpz Number of longitudinal grid points in the fieldmap.
     * @param dir Direction of the derivative (DX, DY, DZ).
     * 
     * @note We need to pass the fieldmap data as arguments to this function 
     * because the 'this' pointer cannot be captured in a GPU kernel, 
     * and thus we cannot access member variables directly.
     */
    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static bool computeFieldDerivative(
        const Vector_t<double, 3>& R,
        Vector_t<double, 3>& B,
        const ViewType& Bz,
        const ViewType& Br,
        double hr,
        double hz,
        double /* zbegin */,
        int num_gridpr,
        int num_gridpz,
        const DiffDirection& dir)
    {
        double BfieldR = 0;
        double BfieldZ = 0;

        double RR = sqrt(R(0) * R(0) + R(1) * R(1));

        int indexr    = (int)floor(RR / hr);
        double leverr = (RR / hr) - indexr;

        int indexz    = (int)floor((R(2)) / hz);
        
        // Unused:
        //double leverz = (R(2) / hz) - indexz;

        // test

        if ((indexz < 0) || (indexz + 2 > num_gridpz))
            return false;

        if (indexr + 2 > num_gridpr)
            return true;

        int index1 = indexz + indexr * num_gridpz;
        int index2 = index1 + num_gridpz;
        if (dir == DZ) {
            if (indexz > 0) {
                if (indexz < num_gridpz - 1) {
                    BfieldR = (1.0 - leverr)
                                * ((Br(index1 - 1) - 2. * Br(index1)
                                    + Br(index1 + 1))
                                        * (R(2) - indexz * hz) / (hz * hz)
                                    + (Br(index1 + 1) - Br(index1 - 1))
                                        / (2. * hz))
                            + leverr
                                    * ((Br(index2 - 1) - 2. * Br(index2)
                                        + Br(index2 + 1))
                                        * (R(2) - indexz * hz) / (hz * hz)
                                    + (Br(index2 + 1) - Br(index2 - 1))
                                            / (2. * hz));

                    BfieldZ = (1.0 - leverr)
                                * ((Bz(index1 - 1) - 2. * Bz(index1)
                                    + Bz(index1 + 1))
                                        * (R(2) - indexz * hz) / (hz * hz)
                                    + (Bz(index1 + 1) - Bz(index1 - 1))
                                        / (2. * hz))
                            + leverr
                                    * ((Bz(index2 - 1) - 2. * Bz(index2)
                                        + Bz(index2 + 1))
                                        * (R(2) - indexz * hz) / (hz * hz)
                                    + (Bz(index2 + 1) - Bz(index2 - 1))
                                            / (2. * hz));
                } else {
                    BfieldR =
                        (1.0 - leverr)
                            * ((Br(index1) - 2. * Br(index1 - 1)
                                + Br(index1 - 2))
                                * (R(2) - indexz * hz) / (hz * hz)
                            - (3. * Br(index1) - 4. * Br(index1 - 1)
                                + Br(index1 - 2))
                                    / (2. * hz))
                        + leverr
                            * ((Br(index2) - 2. * Br(index2 - 1)
                                + Br(index2 - 2))
                                    * (R(2) - indexz * hz) / (hz * hz)
                                - (3. * Br(index2) - 4. * Br(index2 - 1)
                                    + Br(index2 - 2))
                                    / (2. * hz));

                    BfieldZ =
                        (1.0 - leverr)
                            * ((Bz(index1) - 2. * Bz(index1 - 1)
                                + Bz(index1 - 2))
                                * (R(2) - indexz * hz) / (hz * hz)
                            - (3. * Bz(index1) - 4. * Bz(index1 - 1)
                                + Bz(index1 - 2))
                                    / (2. * hz))
                        + leverr
                            * ((Bz(index2) - 2. * Bz(index2 - 1)
                                + Bz(index2 - 2))
                                    * (R(2) - indexz * hz) / (hz * hz)
                                - (3. * Bz(index2) - 4. * Bz(index2 - 1)
                                    + Bz(index2 - 2))
                                    / (2. * hz));
                }
            } else {
                BfieldR =
                    (1.0 - leverr)
                        * ((Br(index1) - 2. * Br(index1 + 1)
                            + Br(index1 + 2))
                            * (R(2) - indexz * hz) / (hz * hz)
                        - (3. * Br(index1) - 4. * Br(index1 + 1)
                            + Br(index1 + 2))
                                / (2. * hz))
                    + leverr
                        * ((Br(index2) - 2. * Br(index2 + 1)
                            + Br(index2 + 2))
                                * (R(2) - indexz * hz) / (hz * hz)
                            - (3. * Br(index2) - 4. * Br(index2 + 1)
                                + Br(index2 + 2))
                                / (2. * hz));

                BfieldZ =
                    (1.0 - leverr)
                        * ((Bz(index1) - 2. * Bz(index1 + 1)
                            + Bz(index1 + 2))
                            * (R(2) - indexz * hz) / (hz * hz)
                        - (3. * Bz(index1) - 4. * Bz(index1 + 1)
                            + Bz(index1 + 2))
                                / (2. * hz))
                    + leverr
                        * ((Bz(index2) - 2. * Bz(index2 + 1)
                            + Bz(index2 + 2))
                            * (R(2) - indexz * hz) / (hz * hz)
                        - (3. * Bz(index2) - 4. * Bz(index2 + 1)
                            + Bz(index2 + 2))
                                / (2. * hz));
            }
        } else {
             // For directions other than DZ, currently not implemented or fall through? 
             // Original logic only handled DZ condition fully in the snippet.
             // If dir != DZ, return false?
             // Original code had `if (dir == DZ)` block, but what if not?
             // The snippet ended abruptly.
             // Assuming no other derivatives supported for magnetostatic map normally (only d/dz implies fringe?) or maybe radial.
             // For now, I'll close the function. If dir != DZ, no B changes.
        }

        if (RR != 0) {
            B(0) += BfieldR * R(0) / RR;
            B(1) += BfieldR * R(1) / RR;
        }
        B(2) += BfieldZ;
        return false;
    }

    /**
     * @brief Apply the FM to all the particles
     * 
     * @param pc Particle container
     */
    void applyField(std::shared_ptr<ParticleContainer_t> pc);

private:
    FM2DMagnetoStatic(std::string aFilename);
    ~FM2DMagnetoStatic();

    virtual void readMap();
    virtual void freeMap();

    /// @brief Fieldstrengths 
    Kokkos::DualView<double*> FieldstrengthBz_m;    
    Kokkos::DualView<double*> FieldstrengthBr_m;    

    /// @brief Radius Bounds
    double rbegin_m;   
    double rend_m;

    /// @brief Z Bounds relative to element edge
    double zbegin_m;
    double zend_m;

    /// @brief Grid 
    double hz_m;                   
    double hr_m;                   
    int num_gridpr_m;              
    int num_gridpz_m;              

    bool swap_m;
    friend class Fieldmap;
};


#endif