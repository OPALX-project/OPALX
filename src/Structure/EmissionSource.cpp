#include "Structure/EmissionSource.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

#include <cmath>

EmissionSource::EmissionSource()
    : Definition(
              SIZE, "EMISSIONSOURCE",
              "The EMISSIONSOURCE statement defines a particle injector with "
              "position R0, momentum P0, start time t0, and a linked Distribution.") {
    itsAttr[DISTRIBUTION] =
            Attributes::makeString("DISTRIBUTION", "Name of the Distribution to sample from.");

    itsAttr[R0X] = Attributes::makeReal("R0X", "Position offset x added after sampling [m].", 0.0);
    itsAttr[R0Y] = Attributes::makeReal("R0Y", "Position offset y added after sampling [m].", 0.0);
    itsAttr[R0Z] = Attributes::makeReal("R0Z", "Position offset z added after sampling [m].", 0.0);

    itsAttr[P0X] = Attributes::makeReal(
            "P0X", "Momentum offset x added after sampling [beta*gamma].", 0.0);
    itsAttr[P0Y] = Attributes::makeReal(
            "P0Y", "Momentum offset y added after sampling [beta*gamma].", 0.0);
    itsAttr[P0Z] = Attributes::makeReal(
            "P0Z", "Momentum offset z added after sampling [beta*gamma].", 0.0);

    itsAttr[T0] = Attributes::makeReal("T0", "Start time when sampling begins [s].", 0.0);

    itsAttr[ZEROFACE_R0Z] = Attributes::makeBool(
            "ZEROFACE_R0Z", "Set Dirichlet boundary conditions (0 potential) in xy plane at R0Z.",
            false);

    itsAttr[SHIFTED_GREENS_FUNCTION] = Attributes::makeBool(
            "SHIFTED_GREENS_FUNCTION",
            "Enforce Dirichlet BCs (0 potential) in xy plane at R0Z via a shifted Green's "
            "function instead of explicit image charges. Mutually exclusive with ZEROFACE_R0Z. "
            "Requires the OPEN field solver. Active for ZEROFACE_MAXSTEPS steps (0 = unlimited).",
            false);

    itsAttr[ZEROFACEPLANEDUMP] = Attributes::makeReal(
            "ZEROFACEPLANEDUMP",
            "Dump interpolated potential on the ZEROFACE_R0Z plane every n-th global "
            "timestep (0 disables dumping, since it's an expensive operation).",
            0.0);

    itsAttr[ZEROFACE_MAXSTEPS] = Attributes::makeReal(
            "ZEROFACE_MAXSTEPS",
            "Number of timesteps for which image charges are active. "
            "After this many steps, the solver continues without image charges. "
            "0 means unlimited (image charges always active).",
            0.0);

    itsAttr[EMISSIONMODEL] = Attributes::makePredefinedString(
            "EMISSIONMODEL",
            "Emission model for cathode emission. "
            "NONE: thermal momentum applied in z only. "
            "ASTRA: 3D isotropic thermal emission on forward half-sphere.",
            {"NONE", "ASTRA"}, "NONE");

    itsAttr[EKIN] = Attributes::makeReal(
            "EKIN",
            "Kinetic emission energy [eV]. The normalized source samples are stretched to "
            "this beta*gamma magnitude before P0 is applied. For EMISSIONMODEL=NONE this "
            "is a longitudinal momentum increment; for EMISSIONMODEL=ASTRA this is the "
            "thermal momentum magnitude used to sample the forward half-sphere.",
            0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

EmissionSource::EmissionSource(const std::string& name, EmissionSource* parent)
    : Definition(name, parent) {}

EmissionSource::~EmissionSource() {}

bool EmissionSource::canReplaceBy(Object* object) {
    return dynamic_cast<EmissionSource*>(object) != nullptr;
}

EmissionSource* EmissionSource::clone(const std::string& name) {
    return new EmissionSource(name, this);
}

void EmissionSource::execute() {}

EmissionSource* EmissionSource::find(const std::string& name) {
    Object* obj = OpalData::getInstance()->find(name);
    auto* es    = dynamic_cast<EmissionSource*>(obj);
    if (!es) {
        throw OpalException("EmissionSource::find()", "EmissionSource \"" + name + "\" not found.");
    }
    return es;
}

std::string EmissionSource::getDistributionName() const {
    return Attributes::getString(itsAttr[DISTRIBUTION]);
}

ippl::Vector<double, 3> EmissionSource::getR0() const {
    return ippl::Vector<double, 3>(
            Attributes::getReal(itsAttr[R0X]), Attributes::getReal(itsAttr[R0Y]),
            Attributes::getReal(itsAttr[R0Z]));
}

ippl::Vector<double, 3> EmissionSource::getP0() const {
    return ippl::Vector<double, 3>(
            Attributes::getReal(itsAttr[P0X]), Attributes::getReal(itsAttr[P0Y]),
            Attributes::getReal(itsAttr[P0Z]));
}

double EmissionSource::getT0() const { return Attributes::getReal(itsAttr[T0]); }

bool EmissionSource::getZeroFaceR0Z() const { return Attributes::getBool(itsAttr[ZEROFACE_R0Z]); }

bool EmissionSource::getShiftedGreensFunction() const {
    return Attributes::getBool(itsAttr[SHIFTED_GREENS_FUNCTION]);
}

int EmissionSource::getZeroFacePlaneDumpFrequency() const {
    const double rawFrequency = Attributes::getReal(itsAttr[ZEROFACEPLANEDUMP]);
    const int frequency       = static_cast<int>(rawFrequency);
    if (rawFrequency < 0.0 || std::floor(rawFrequency) != rawFrequency) {
        throw OpalException(
                "EmissionSource::getZeroFacePlaneDumpFrequency",
                "ZEROFACEPLANEDUMP must be a non-negative integer value.");
    }
    return frequency;
}

int EmissionSource::getZerofaceMaxSteps() const {
    const double rawValue = Attributes::getReal(itsAttr[ZEROFACE_MAXSTEPS]);
    const int value       = static_cast<int>(rawValue);
    if (rawValue < 0.0 || std::floor(rawValue) != rawValue) {
        throw OpalException(
                "EmissionSource::getZerofaceMaxSteps",
                "ZEROFACE_MAXSTEPS must be a non-negative integer value.");
    }
    return value;
}

std::string EmissionSource::getEmissionModel() const {
    return Attributes::getString(itsAttr[EMISSIONMODEL]);
}

double EmissionSource::getKineticEnergy() const {
    return std::abs(Attributes::getReal(itsAttr[EKIN]));
}
