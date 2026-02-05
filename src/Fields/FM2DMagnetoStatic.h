#ifndef CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH
#define CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

class FM2DMagnetoStatic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t<double, 3> &R, Vector_t<double, 3> &E, Vector_t<double, 3> &B) const;
    virtual bool getFieldDerivative(const Vector_t<double, 3> &R, Vector_t<double, 3> &E, Vector_t<double, 3> &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    KOKKOS_INLINE_FUNCTION
    virtual bool isInside(const Vector_t<double, 3> &r) const;

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

    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static bool computeFieldDerivative(
        const Vector_t<double, 3>& R,   // Position
        Vector_t<double, 3>& B,         // Derivative of B
        const ViewType& Bz,             // Longitudinal Value of FM 
        const ViewType& Br,             // Radial Value of FM
        double hr,                      // Radial Grid Spacing 
        double hz,                      // Longitudinal Grid Spacing 
        double /* zbegin */,            // Start Relative to Element Edge 
        int num_gridpr,                 // # Radial Gridpoints 
        int num_gridpz,                 // # Longitudinal Gridpoints
        const DiffDirection& dir)       // Direction of Derivative 
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
     * @param Rview View of particle positions
     * @param Eview View of E-field at particle positions
     * @param Bview View of B-field at particle positions
     */
    void applyField(std::shared_ptr<ParticleContainer_t> pc)
    {
        // local variables to copy to the kernel
        double zbegin = zbegin_m;
        double zend = zend_m;
        double rend = rend_m;
        double hr = hr_m;
        double hz = hz_m;
        int num_gridpr = num_gridpr_m;
        int num_gridpz = num_gridpz_m;

        // capture device views
        auto Bz_device = FieldstrengthBz_m.view_device();
        auto Br_device = FieldstrengthBr_m.view_device();

        // capture views
        auto Rview = pc->R.getView();
        auto Bview = pc->B.getView();

        Kokkos::parallel_for("FM2DMagnetoStatic::applyField",
        ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const int i)
        {
            // Check bounds
            if(Rview(i)(2) >= zbegin &&
                Rview(i)(2) < zend &&
                sqrt(Rview(i)(0)*Rview(i)(0) + Rview(i)(1)*Rview(i)(1)) < rend) 
            {
                computeField(Rview(i), 
                    Bview(i), 
                    Bz_device,           // Use device view
                    Br_device,           // Use device view
                    hr, hz, zbegin, num_gridpr, num_gridpz);
            }
        });
    }

private:
    FM2DMagnetoStatic(std::string aFilename);
    ~FM2DMagnetoStatic();

    virtual void readMap();
    virtual void freeMap();

    Kokkos::DualView<double*> FieldstrengthBz_m;    /**< 2D array with Ez, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
    Kokkos::DualView<double*> FieldstrengthBr_m;    /**< 2D array with Er, read in like Ez*/

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    double hr_m;                   /**< length between points in grid, r-direction, m*/
    int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

/**
 * @brief Checks if position r is inside the bounds
 * 
 * @todo This does not work on gpu since zbegin_m, zend_m and rend_m is
 * not accessible on GPU
 */
KOKKOS_INLINE_FUNCTION
bool FM2DMagnetoStatic::isInside(const Vector_t<double, 3> &r) const
{
    return r(2) >= zbegin_m && r(2) < zend_m && sqrt(r(0)*r(0) + r(1)*r(1)) < rend_m;
}

#endif