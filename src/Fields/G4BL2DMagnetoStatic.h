#ifndef OPALX_FIELDMAPG4BL2DMAGNETOSTATIC_HH
#define OPALX_FIELDMAPG4BL2DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <string>

/**
 * @class G4BL2DMagnetoStatic
 * @brief Reader for G4beamline `cylinder` field maps.
 *
 * The file holds an axisymmetric (r, z) magnetic field:
 * @code
 * param normB=1. current=1.
 * cylinder Z0=-1500.0000 nZ=301 dZ=10.0000 nR=51 dR=10.0000
 * Bz
 * <nZ rows of nR values>
 * Br
 * <nZ rows of nR values>
 * @endcode
 * Positions are in mm, field values in Tesla. Within a block the row index is
 * the z index and the column index is the r index; column 0 is on axis. r
 * always starts at 0, so the `cylinder` line carries no R0.
 *
 * @note Unlike every other OPALX ASCII reader this one does @b not normalize.
 * The values are absolute Tesla, so `KS` on the element is a plain multiplier.
 * G4beamline evaluates B * normB * current_deck / current_file, where
 * current_file is the `param current=` of the file and current_deck the
 * `current=` of the `fieldmap` command. Both halves that live in the file,
 * normB and current_file, are folded in on load, so `KS` is exactly
 * G4beamline's deck-side `current=` and `KS = 1` matches a bare
 * `fieldmap ... ` placement.
 *
 * @note The grid layout and the interpolation are identical to
 * FM2DMagnetoStatic, so FM2DMagnetoStatic::computeField() is reused rather
 * than duplicated. Only the parser is new.
 */
class G4BL2DMagnetoStatic : public Fieldmap {
public:
    /**
     * @brief Get the field strength at a given point.
     *
     * @param R Position [m] relative to the field map origin.
     * @param E Output Electric field [MV/m], untouched.
     * @param B Output Magnetic field [T].
     * @return true if R is outside of the field map, false otherwise.
     */
    virtual bool getFieldstrength(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) const override;

    /**
     * @brief Get the field derivative with respect to a direction.
     * @note Not implemented, throws.
     */
    virtual bool getFieldDerivative(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B,
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
            double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
            double& zFinal) const override;

    /**
     * @brief Swap coordinates.
     * @note No-op. The G4beamline cylinder format has a fixed block layout, so
     * there is no orientation to swap.
     */
    virtual void swap() override;

    /// @brief Print info about the field map.
    virtual void getInfo(Inform* msg) override;

    /**
     * @brief Get the frequency.
     * @note Not implemented, throws.
     */
    virtual double getFrequency() const override;

    /**
     * @brief Set the frequency.
     * @note Not implemented, throws.
     */
    virtual void setFrequency(double freq) override;

    /**
     * @brief Checks if the given coordinate is inside the volume covered by the
     * fieldmap
     * @param r Coordinate
     * @note This cannot be called inside a GPU kernel (implicit capture of the
     * 'this' pointer not allowed on device)
     */
    bool isInside(const Vector_t<double, 3>& r) const override {
        return r(2) >= zbegin_m && r(2) < zend_m && sqrt(r(0) * r(0) + r(1) * r(1)) < rend_m;
    }

    /**
     * @brief Apply the FM to all the particles
     *
     * @param pc Particle container
     * @param scale Scaling factor applied to the field
     */
    void applyField(std::shared_ptr<ParticleContainer_t> pc, double scale = 1.0) override;

private:
    /**
     * @param aFilename Path to the .g4blmap file
     * @param zReverse Mirror the map in z and negate Bz on load, which is what
     * G4beamline's `place ... rotation=Y180` does to an axisymmetric map
     */
    G4BL2DMagnetoStatic(std::string aFilename, bool zReverse = false);
    ~G4BL2DMagnetoStatic();

    void readMap() override;
    void freeMap() override;

    /**
     * @brief Read the whole file, validating as it goes.
     *
     * Called twice: from the constructor with @p loadData false, to parse the
     * header and check the file is well formed, and from readMap() with
     * @p loadData true, to fill the views. One function so the two passes
     * cannot drift apart.
     *
     * @param loadData If true, store the values into FieldstrengthB*_m, which
     * must already be allocated.
     */
    void parseFile(bool loadData);

    /// @brief Fieldstrengths [T], indexed indexz + indexr * num_gridpz_m
    Kokkos::DualView<double*> FieldstrengthBz_m;
    Kokkos::DualView<double*> FieldstrengthBr_m;

    /// @brief Radius bounds [m]
    double rbegin_m;
    double rend_m;

    /// @brief Z bounds relative to element edge [m], after any reversal
    double zbegin_m;
    double zend_m;

    /// @brief Grid spacing [m]
    double hz_m;
    double hr_m;

    /// @brief Number of grid points
    int num_gridpr_m;
    int num_gridpz_m;

    /// @brief Field scale factor from the file's `param` line, normB / current
    double fieldScale_m;

    friend class Fieldmap;
};

#endif
