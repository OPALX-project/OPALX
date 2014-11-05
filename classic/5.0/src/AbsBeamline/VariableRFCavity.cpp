/*
 *  Copyright (c) 2014, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "AbsBeamline/VariableRFCavity.h"

#include "Physics/Physics.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/GeneralClassicException.h"

VariableRFCavity::VariableRFCavity(const std::string &name) : Component(name) {
    initNull();  // initialise everything to NULL
}

VariableRFCavity::VariableRFCavity() : Component() {
    initNull();  // initialise everything to NULL
}

VariableRFCavity::VariableRFCavity(const VariableRFCavity& var) : Component() {
    initNull();  // initialise everything to NULL
    *this = var;
}

VariableRFCavity& VariableRFCavity::operator=(const VariableRFCavity& rhs) {
    if (&rhs == this) {
        return *this;
    }
    setName(rhs.getName());
    setPhaseModel(NULL);
    setAmplitudeModel(NULL);
    setFrequencyModel(NULL);
    if (rhs._phase_td != NULL)
        setPhaseModel(std::shared_ptr<AbstractTimeDependence>(rhs._phase_td->clone()));
    if (rhs._amplitude_td != NULL) {
        setAmplitudeModel(std::shared_ptr<AbstractTimeDependence>(rhs._amplitude_td->clone()));
    }
    if (rhs._frequency_td != NULL)
        setFrequencyModel(std::shared_ptr<AbstractTimeDependence>(rhs._frequency_td->clone()));
    phaseName_m = rhs.phaseName_m;
    amplitudeName_m = rhs.amplitudeName_m;
    frequencyName_m = rhs.frequencyName_m;
    halfWidth_m = rhs.halfWidth_m;
    halfHeight_m = rhs.halfHeight_m;
    setLength(rhs._length);
    return *this;
}

VariableRFCavity::~VariableRFCavity() {
    // if (_phase_td != NULL)
    //     delete _phase_td;
    // if (_amplitude_td != NULL)
    //     delete _amplitude_td;
    // if (_frequency_td != NULL)
    //     delete _frequency_td;
}

void VariableRFCavity::initNull() {
  _length = 0.;
  // _phase_td = NULL;
  // _amplitude_td = NULL;
  // _frequency_td = NULL;
  phaseName_m = "";
  amplitudeName_m = "";
  frequencyName_m = "";
  halfHeight_m = 0.;
  halfWidth_m = 0;
  RefPartBunch_m = NULL;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getAmplitudeModel() const {
    return _amplitude_td;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getPhaseModel() const {
    return _phase_td;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getFrequencyModel() const {
    return _frequency_td;
}

void VariableRFCavity::setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitude_td) {
    // if (_amplitude_td != NULL && amplitude_td != _amplitude_td)
    //     delete _amplitude_td;
    _amplitude_td = amplitude_td;
}

void VariableRFCavity::setPhaseModel(std::shared_ptr<AbstractTimeDependence> phase_td) {
    // if (_phase_td != NULL && phase_td != _phase_td)
    //     delete _phase_td;
    _phase_td = phase_td;
}

void VariableRFCavity::setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequency_td) {
    // if (_frequency_td != NULL && frequency_td != _frequency_td)
    //     delete _frequency_td;
    _frequency_td = frequency_td;
}

StraightGeometry &VariableRFCavity::getGeometry() {
    return geometry;
}

const StraightGeometry &VariableRFCavity::getGeometry() const {
    return geometry;
}

EMField &VariableRFCavity::getField() {
  throw GeneralClassicException("VariableRFCavity",
                      "No field defined for VariableRFCavity");
}

const EMField &VariableRFCavity::getField() const {
  throw GeneralClassicException("VariableRFCavity",
                      "No field defined for VariableRFCavity");
}


bool VariableRFCavity::apply(const size_t &i, const double &t,
                    Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->get_centroid(), t,
                 E, B);
}

bool VariableRFCavity::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);

    if (apply(i, t, Ev, Bv))
        return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

// If this is too slow: a quicker implementation would be to use templates not
// inheritance (vtable lookup is removed). This is in the inner
// tracking loop, so low level optimisation is possibly worthwhile.
//
// Do I need bound checking here? I have no "radius" parameter, but I do have a
// "length".
bool VariableRFCavity::apply(const Vector_t &R, const Vector_t &centroid,
           const double &t, Vector_t &E, Vector_t &B) {
    // std::cerr << "VariableRFCavity::apply " << R[0] << " " << R[1] << " " << R[2] << " * " << halfWidth_m << " " << halfHeight_m << std::endl;
    if (R[2] < 0. || R[2] > _length ||
        fabs(R[0]) > halfWidth_m || fabs(R[1]) > halfHeight_m)
        return true;
    double E0 = _amplitude_td->getValue(t);
    double f = _frequency_td->getValue(t);
    double phi = _phase_td->getValue(t);
    E = Vector_t(0., 0., E0*sin(Physics::two_pi*f*t+phi));
    // std::cerr << "                        t: " << t << " f: " << f << " phi: " << phi << " E0: " << E0 << " E[2]: " << E[2] << std::endl;
    return false;
}

void VariableRFCavity::initialise(PartBunch *bunch, double &startField, double &endField,
                const double &scaleFactor) {
    RefPartBunch_m = bunch;
}

void VariableRFCavity::finalise() {
    RefPartBunch_m = NULL;
}

ElementBase* VariableRFCavity::clone() const {
    return new VariableRFCavity(*this);
}

void VariableRFCavity::accept(BeamlineVisitor& visitor) const {
    VariableRFCavity* cavity = const_cast<VariableRFCavity*>(this);
    std::shared_ptr<AbstractTimeDependence> phaseTD =
                    AbstractTimeDependence::getTimeDependence(phaseName_m);
    cavity->setPhaseModel(std::shared_ptr<AbstractTimeDependence>(phaseTD->clone()));
    std::shared_ptr<AbstractTimeDependence> frequencyTD =
                    AbstractTimeDependence::getTimeDependence(frequencyName_m);
    cavity->setFrequencyModel(std::shared_ptr<AbstractTimeDependence>(frequencyTD->clone()));
        std::shared_ptr<AbstractTimeDependence> amplitudeTD =
                    AbstractTimeDependence::getTimeDependence(amplitudeName_m);
    cavity->setAmplitudeModel(std::shared_ptr<AbstractTimeDependence>(amplitudeTD->clone()));
    visitor.visitVariableRFCavity(*this);

    if (halfHeight_m < 1e-9 || halfWidth_m < 1e-9)
        throw GeneralClassicException("VariableRFCavity::accept",
                            "Height or width was not set on VariableRFCavity");
}

void VariableRFCavity::setLength(double length) {
    _length = length;
    geometry.setElementLength(length);
}
