#include "AbsBeamline/Laser.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/LinearCompton.h"

Laser::Laser() : Laser("") {}

Laser::Laser(const std::string& name)
    : Component(name),
      startField_m(0.0),
      wavelength_m(0.0),
      pulseEnergy_m(0.0),
      pulseLength_m(0.0),
      waistX_m(0.0),
      waistY_m(0.0),
      direction_m(0.0),
      stokes_m(0.0) {}

Laser::Laser(const Laser& right)
    : Component(right),
      startField_m(right.startField_m),
      wavelength_m(right.wavelength_m),
      pulseEnergy_m(right.pulseEnergy_m),
      pulseLength_m(right.pulseLength_m),
      waistX_m(right.waistX_m),
      waistY_m(right.waistY_m),
      direction_m(right.direction_m),
      stokes_m(right.stokes_m) {}

Laser::~Laser() {}

void Laser::accept(BeamlineVisitor& visitor) const { visitor.visitLaser(*this); }

void Laser::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    endField       = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m   = startField;
}

void Laser::finalise() {}

bool Laser::bends() const { return false; }

void Laser::getFieldExtend(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = startField_m + getElementLength();
}

ElementType Laser::getType() const { return ElementType::LASER; }

void Laser::setWavelength(double wavelength) { wavelength_m = wavelength; }

double Laser::getWavelength() const { return wavelength_m; }

void Laser::setPulseEnergy(double pulseEnergy) { pulseEnergy_m = pulseEnergy; }

double Laser::getPulseEnergy() const { return pulseEnergy_m; }

void Laser::setPulseLength(double pulseLength) { pulseLength_m = pulseLength; }

double Laser::getPulseLength() const { return pulseLength_m; }

void Laser::setWaistX(double waistX) { waistX_m = waistX; }

double Laser::getWaistX() const { return waistX_m; }

void Laser::setWaistY(double waistY) { waistY_m = waistY; }

double Laser::getWaistY() const { return waistY_m; }

void Laser::setDirection(const Vector_t<double, 3>& direction) { direction_m = direction; }

const Vector_t<double, 3>& Laser::getDirection() const { return direction_m; }

void Laser::setStokes(const Vector_t<double, 3>& stokes) { stokes_m = stokes; }

double Laser::getPhotonEnergyGeV() const {
    return Physics::LinearCompton::photonEnergyFromWavelengthGeV(wavelength_m);
}

double Laser::getLinearComptonInvariantX(
        double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection) const {
    return Physics::LinearCompton::invariantX(
            electronTotalEnergyGeV, getPhotonEnergyGeV(), beamDirection, direction_m);
}

double Laser::getLinearComptonForwardPhotonEnergyGeV(
        double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection) const {
    return Physics::LinearCompton::labForwardPhotonEnergyGeV(
            electronTotalEnergyGeV, getPhotonEnergyGeV(), beamDirection, direction_m);
}

const Vector_t<double, 3>& Laser::getStokes() const { return stokes_m; }
