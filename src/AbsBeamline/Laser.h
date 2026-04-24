#ifndef OPALX_Laser_HH
#define OPALX_Laser_HH

#include "AbsBeamline/Component.h"

/**
 * @brief Passive OPALX laser element.
 *
 * The laser currently acts as a beamline component that stores an analytic laser
 * pulse definition without applying fields or transport kicks. Besides the basic
 * geometric/component interface, it exposes a small set of helper functions for
 * the first linear Compton benchmarks used in the unit tests.
 *
 * Stored quantities follow the current OPALX conventions:
 * - wavelength in m
 * - pulse energy in J
 * - pulse length in s
 * - waists in m
 * - direction as a normalized laboratory-frame unit vector
 * - Stokes vector as a dimensionless polarization vector
 *
 * The Compton helpers operate on electron total energy in GeV and return photon
 * energies in GeV.
 */
class Laser : public Component {
public:
    /** @brief Construct an unnamed laser element. */
    Laser();

    /** @brief Construct a laser element with the given name. */
    explicit Laser(const std::string& name);

    /** @brief Copy-construct a laser element. */
    Laser(const Laser& right);

    /** @brief Destroy the laser element. */
    ~Laser() override;

    /** @brief Accept a beamline visitor through the generic component interface. */
    void accept(BeamlineVisitor&) const override;

    /**
     * @brief Initialize the element for a tracking pass.
     * @param bunch Referenced bunch.
     * @param startField Element entrance position in the current lattice traversal.
     * @param endField Filled with the element exit position.
     */
    void initialise(PartBunch_t* bunch, double& startField, double& endField) override;

    /** @brief Finalize the element after a tracking pass. */
    void finalise() override;

    /** @brief Return false because the element has straight reference geometry. */
    bool bends() const override;

    /** @brief Get the longitudinal element extent in the current lattice traversal. */
    void getFieldExtend(double& zBegin, double& zEnd) const override;

    /** @brief Return the OPALX element type identifier. */
    ElementType getType() const override;

    /** @brief Return the minimum number of time steps required by this passive element. */
    int getRequiredNumberOfTimeSteps() const override;

    /** @brief Set the laser wavelength in m. */
    void setWavelength(double wavelength);

    /** @brief Get the laser wavelength in m. */
    double getWavelength() const;

    /** @brief Set the total pulse energy in J. */
    void setPulseEnergy(double pulseEnergy);

    /** @brief Get the total pulse energy in J. */
    double getPulseEnergy() const;

    /** @brief Set the pulse duration in s. */
    void setPulseLength(double pulseLength);

    /** @brief Get the pulse duration in s. */
    double getPulseLength() const;

    /** @brief Set the horizontal waist in m. */
    void setWaistX(double waistX);

    /** @brief Get the horizontal waist in m. */
    double getWaistX() const;

    /** @brief Set the vertical waist in m. */
    void setWaistY(double waistY);

    /** @brief Get the vertical waist in m. */
    double getWaistY() const;

    /** @brief Set the laboratory-frame laser propagation direction. */
    void setDirection(const Vector_t<double, 3>& direction);

    /** @brief Get the laboratory-frame laser propagation direction. */
    const Vector_t<double, 3>& getDirection() const;

    /** @brief Set the normalized Stokes polarization vector. */
    void setStokes(const Vector_t<double, 3>& stokes);

    /** @brief Get the normalized Stokes polarization vector. */
    const Vector_t<double, 3>& getStokes() const;

    /**
     * @brief Convert the stored laser wavelength to a single-photon energy.
     * @return Laser photon energy in GeV.
     * @throws OpalException If the stored wavelength is non-positive.
     */
    double getPhotonEnergyGeV() const;

    /**
     * @brief Compute the linear Compton invariant x for a given beam direction.
     *
     * This follows the CAIN linear-Compton rest-frame mapping
     * @f$\omega_1 = \gamma \omega_L (1 - \beta \cos \alpha)@f$ and returns
     * @f$x = 2 \omega_1 / m_e@f$.
     *
     * @param electronTotalEnergyGeV Electron total energy in GeV.
     * @param beamDirection Laboratory-frame unit direction of the incoming electron beam.
     * @return Dimensionless linear Compton invariant @f$x@f$.
     * @throws OpalException If the electron energy is below rest energy or if either
     *         direction vector is zero.
     */
    double getLinearComptonInvariantX(
            double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection) const;

    /**
     * @brief Compute the forward scattered photon energy for linear Compton scattering.
     *
     * The helper uses the same CAIN-inspired sequence as the first OPALX benchmark:
     * boost the laser photon into the electron rest frame, apply the Compton recoil,
     * and boost the forward photon back to the laboratory frame.
     *
     * @param electronTotalEnergyGeV Electron total energy in GeV.
     * @param beamDirection Laboratory-frame unit direction of the incoming electron beam.
     * @return Forward scattered photon energy in GeV.
     * @throws OpalException If the electron energy is below rest energy or if either
     *         direction vector is zero.
     */
    double getLinearComptonForwardPhotonEnergyGeV(
            double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection) const;

private:
    double startField_m;
    double wavelength_m;
    double pulseEnergy_m;
    double pulseLength_m;
    double waistX_m;
    double waistY_m;
    Vector_t<double, 3> direction_m;
    Vector_t<double, 3> stokes_m;
};

inline int Laser::getRequiredNumberOfTimeSteps() const { return 1; }

#endif  // OPALX_Laser_HH
