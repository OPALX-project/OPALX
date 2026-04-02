#include "AbsBeamline/Laser.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>

Laser::Laser()
    : Laser("") {
}

Laser::Laser(const std::string& name)
    : Component(name),
      startField_m(0.0),
      wavelength_m(0.0),
      pulseEnergy_m(0.0),
      pulseLength_m(0.0),
      waistX_m(0.0),
      waistY_m(0.0),
      direction_m(0.0),
      stokes_m(0.0) {
}

Laser::Laser(const Laser& right)
    : Component(right),
      startField_m(right.startField_m),
      wavelength_m(right.wavelength_m),
      pulseEnergy_m(right.pulseEnergy_m),
      pulseLength_m(right.pulseLength_m),
      waistX_m(right.waistX_m),
      waistY_m(right.waistY_m),
      direction_m(right.direction_m),
      stokes_m(right.stokes_m) {
}

Laser::~Laser() {
}

namespace {
    // Normalize a direction-like vector and fail fast on invalid zero input so the
    // public Compton helpers report a clear argument error.
    Vector_t<double, 3> normalizeDirection(const Vector_t<double, 3>& direction,
                                           const char* where,
                                           const char* argument) {
        const double norm = std::sqrt(dot(direction, direction));
        if (norm <= 0.0) {
            throw OpalException(where, std::string("\"") + argument + "\" must be a non-zero vector.");
        }

        return direction / norm;
    }

    // Convert total electron energy to beta using OPALX GeV units and reject
    // unphysical energies below the electron rest energy.
    double getElectronBeta(double electronTotalEnergyGeV, const char* where) {
        if (electronTotalEnergyGeV <= Physics::m_e) {
            throw OpalException(where, "Electron total energy must be larger than the electron rest energy.");
        }

        const double gamma = electronTotalEnergyGeV / Physics::m_e;
        return std::sqrt(1.0 - 1.0 / (gamma * gamma));
    }

    // Clamp roundoff-sensitive dot products before they enter scattering-angle formulas.
    double clampCosine(double value) {
        return std::max(-1.0, std::min(1.0, value));
    }
}

void Laser::accept(BeamlineVisitor& visitor) const {
    visitor.visitComponent(*this);
}

void Laser::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}

void Laser::finalise() {
}

bool Laser::bends() const {
    return false;
}

void Laser::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementType Laser::getType() const {
    return ElementType::LASER;
}

void Laser::setWavelength(double wavelength) {
    wavelength_m = wavelength;
}

double Laser::getWavelength() const {
    return wavelength_m;
}

void Laser::setPulseEnergy(double pulseEnergy) {
    pulseEnergy_m = pulseEnergy;
}

double Laser::getPulseEnergy() const {
    return pulseEnergy_m;
}

void Laser::setPulseLength(double pulseLength) {
    pulseLength_m = pulseLength;
}

double Laser::getPulseLength() const {
    return pulseLength_m;
}

void Laser::setWaistX(double waistX) {
    waistX_m = waistX;
}

double Laser::getWaistX() const {
    return waistX_m;
}

void Laser::setWaistY(double waistY) {
    waistY_m = waistY;
}

double Laser::getWaistY() const {
    return waistY_m;
}

void Laser::setDirection(const Vector_t<double, 3>& direction) {
    direction_m = direction;
}

const Vector_t<double, 3>& Laser::getDirection() const {
    return direction_m;
}

void Laser::setStokes(const Vector_t<double, 3>& stokes) {
    stokes_m = stokes;
}

double Laser::getPhotonEnergyGeV() const {
    constexpr const char* where = "Laser::getPhotonEnergyGeV()";

    if (wavelength_m <= 0.0) {
        throw OpalException(where, "Laser wavelength must be greater than 0.");
    }

    return Physics::two_pi * Physics::h_bar * Physics::c / wavelength_m;
}

double Laser::getLinearComptonInvariantX(double electronTotalEnergyGeV,
                                         const Vector_t<double, 3>& beamDirection) const {
    constexpr const char* where = "Laser::getLinearComptonInvariantX()";

    const double beta = getElectronBeta(electronTotalEnergyGeV, where);
    const double gamma = electronTotalEnergyGeV / Physics::m_e;
    const Vector_t<double, 3> normalizedBeamDirection =
        normalizeDirection(beamDirection, where, "beamDirection");
    const Vector_t<double, 3> normalizedLaserDirection =
        normalizeDirection(direction_m, where, "direction");

    // CAIN lncpgn.f uses W1 = gamma * W * (1 - beta * cos(alpha)) in the electron rest frame.
    const double cosAlpha = clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));
    const double laserPhotonEnergy = getPhotonEnergyGeV();
    const double restFramePhotonEnergy = gamma * laserPhotonEnergy * (1.0 - beta * cosAlpha);

    return 2.0 * restFramePhotonEnergy / Physics::m_e;
}

double Laser::getLinearComptonForwardPhotonEnergyGeV(
    double electronTotalEnergyGeV,
    const Vector_t<double, 3>& beamDirection) const {
    constexpr const char* where = "Laser::getLinearComptonForwardPhotonEnergyGeV()";

    const double beta = getElectronBeta(electronTotalEnergyGeV, where);
    const double gamma = electronTotalEnergyGeV / Physics::m_e;
    const Vector_t<double, 3> normalizedBeamDirection =
        normalizeDirection(beamDirection, where, "beamDirection");
    const Vector_t<double, 3> normalizedLaserDirection =
        normalizeDirection(direction_m, where, "direction");

    const double cosAlpha = clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));
    const double laserPhotonEnergy = getPhotonEnergyGeV();

    // Follow the CAIN linear Compton kinematics in lncpgn.f: boost the incoming laser photon
    // into the electron rest frame, apply the Compton energy shift there, then boost the
    // forward-scattered photon back to the lab frame.
    const double restFramePhotonEnergy = gamma * laserPhotonEnergy * (1.0 - beta * cosAlpha);
    const double scatteringCosine = (cosAlpha - beta) / (1.0 - beta * cosAlpha);
    const double recoilFactor =
        1.0 + restFramePhotonEnergy / Physics::m_e * (1.0 - scatteringCosine);
    const double restFrameScatteredPhotonEnergy = restFramePhotonEnergy / recoilFactor;

    return gamma * (1.0 + beta) * restFrameScatteredPhotonEnergy;
}

const Vector_t<double, 3>& Laser::getStokes() const {
    return stokes_m;
}
