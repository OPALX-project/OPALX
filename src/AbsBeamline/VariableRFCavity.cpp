//
// Class VariableRFCavity
//   Defines the abstract interface for a RF Cavity
//   with Time Dependent Parameters.
//
// Copyright (c) 2014 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"

VariableRFCavity::VariableRFCavity(const std::string& name) : Component(name) {
    initNull();  // initialise everything to nullptr
}

VariableRFCavity::VariableRFCavity() {
    initNull();  // initialise everything to nullptr
}

VariableRFCavity::VariableRFCavity(const VariableRFCavity& var) : Component(var) {
    initNull();  // initialise everything to nullptr
    *this = var;
}

VariableRFCavity& VariableRFCavity::operator=(const VariableRFCavity& rhs) {
    if (&rhs == this) {
        return *this;
    }
    setName(rhs.getName());
    setPhaseModel(nullptr);
    setAmplitudeModel(nullptr);
    setFrequencyModel(nullptr);
    if (rhs.phaseTD_m != nullptr) {
        setPhaseModel(std::shared_ptr<AbstractTimeDependence>(rhs.phaseTD_m->clone()));
    }
    if (rhs.amplitudeTD_m != nullptr) {
        setAmplitudeModel(std::shared_ptr<AbstractTimeDependence>(rhs.amplitudeTD_m->clone()));
    }
    if (rhs.frequencyTD_m != nullptr) {
        setFrequencyModel(std::shared_ptr<AbstractTimeDependence>(rhs.frequencyTD_m->clone()));
    }
    phaseName_m     = rhs.phaseName_m;
    amplitudeName_m = rhs.amplitudeName_m;
    frequencyName_m = rhs.frequencyName_m;
    halfWidth_m     = rhs.halfWidth_m;
    halfHeight_m    = rhs.halfHeight_m;
    setLength(rhs.length_m);
    return *this;
}

void VariableRFCavity::initNull() {
    length_m        = 0.;
    phaseName_m     = "";
    amplitudeName_m = "";
    frequencyName_m = "";
    halfHeight_m    = 0.;
    halfWidth_m     = 0;
    RefPartBunch_m  = nullptr;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getAmplitudeModel() const {
    return amplitudeTD_m;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getPhaseModel() const {
    return phaseTD_m;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getFrequencyModel() const {
    return frequencyTD_m;
}

void VariableRFCavity::setAmplitudeModel(
        const std::shared_ptr<AbstractTimeDependence> amplitude_td) {
    amplitudeTD_m = amplitude_td;
}

void VariableRFCavity::setPhaseModel(const std::shared_ptr<AbstractTimeDependence> phase_td) {
    phaseTD_m = phase_td;
}

void VariableRFCavity::setFrequencyModel(
        const std::shared_ptr<AbstractTimeDependence> frequency_td) {
    frequencyTD_m = frequency_td;
}

StraightGeometry& VariableRFCavity::getGeometry() { return geometry; }

const StraightGeometry& VariableRFCavity::getGeometry() const { return geometry; }

EMField& VariableRFCavity::getField() {
    throw GeneralClassicException("VariableRFCavity", "No field defined for VariableRFCavity");
}

const EMField& VariableRFCavity::getField() const {
    throw GeneralClassicException(
            "VariableRFCavity::getField", "No field defined for VariableRFCavity");
}

bool VariableRFCavity::apply(
        const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    const auto pc = RefPartBunch_m->getParticleContainer();
    Vector_t<double, 3> R{};
    Kokkos::deep_copy(
            Kokkos::View<Vector_t<double, 3>, Kokkos::HostSpace>(&R),
            Kokkos::subview(pc->R.getView(), i));
    Kokkos::fence();
    const double E0        = amplitudeTD_m->getValue(t) * Units::MVpm2Vpm;
    const double integralF = frequencyTD_m->getIntegral(t) * Units::MHz2Hz;
    const double phi       = phaseTD_m->getValue(t);
    return computeField(R, E, E0, integralF, phi, halfWidth_m, halfHeight_m);
}

bool VariableRFCavity::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    const double E0        = amplitudeTD_m->getValue(t) * Units::MVpm2Vpm;
    const double integralF = frequencyTD_m->getIntegral(t) * Units::MHz2Hz;
    const double phi       = phaseTD_m->getValue(t);
    return computeField(R, E, E0, integralF, phi, halfWidth_m, halfHeight_m);
}

bool VariableRFCavity::apply() {
    const auto pc          = RefPartBunch_m->getParticleContainer();
    const auto R           = pc->R.getView();
    const auto E           = pc->E.getView();
    const auto t           = RefPartBunch_m->getT();
    const double E0        = amplitudeTD_m->getValue(t) * Units::MVpm2Vpm;
    const double integralF = frequencyTD_m->getIntegral(t) * Units::MHz2Hz;
    const double phi       = phaseTD_m->getValue(t);
    const auto count       = pc->getLocalNum();
    const auto halfWidth   = halfWidth_m;
    const auto halfHeight  = halfHeight_m;
    // Kernel launch over all particles
    Kokkos::parallel_for(
            "VariableRFCavity::computeField()", count, KOKKOS_LAMBDA(const size_t i) {
                computeField(R(i), E(i), E0, integralF, phi, halfWidth, halfHeight);
            });

    return false;
}

bool VariableRFCavity::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    return apply(R, P, t, E, B);
}

void VariableRFCavity::initialise(
        PartBunch_t* bunch, double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
}

void VariableRFCavity::finalise() { RefPartBunch_m = nullptr; }

ElementBase* VariableRFCavity::clone() const { return new VariableRFCavity(*this); }

void VariableRFCavity::accept(BeamlineVisitor& visitor) const {
    initialiseTimeDependencies();
    visitor.visitVariableRFCavity(*this);
}

void VariableRFCavity::initialiseTimeDependencies() const {
    const std::shared_ptr<AbstractTimeDependence> phaseTD =
            AbstractTimeDependence::getTimeDependence(phaseName_m);
    phaseTD_m = std::shared_ptr<AbstractTimeDependence>(phaseTD->clone());
    const std::shared_ptr<AbstractTimeDependence> frequencyTD =
            AbstractTimeDependence::getTimeDependence(frequencyName_m);
    frequencyTD_m = std::shared_ptr<AbstractTimeDependence>(frequencyTD->clone());
    const std::shared_ptr<AbstractTimeDependence> amplitudeTD =
            AbstractTimeDependence::getTimeDependence(amplitudeName_m);
    amplitudeTD_m = std::shared_ptr<AbstractTimeDependence>(amplitudeTD->clone());

    if (halfHeight_m < 1e-9 || halfWidth_m < 1e-9)
        throw GeneralClassicException(
                "VariableRFCavity::initialise", "Height or width was not set on VariableRFCavity");
}

void VariableRFCavity::setLength(const double length) {
    length_m = length;
    geometry.setElementLength(length_m);
}
