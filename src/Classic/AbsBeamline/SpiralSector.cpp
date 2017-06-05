/*
 *  Copyright (c) 2017, Chris Rogers
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

#include <math.h>

#include "AbsBeamline/SpiralSector.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"

SpiralSector::SpiralSector(const std::string &name)
        : Component(name),
         planarArcGeometry_m(1., 1.), dummy() {
    setElType(isDrift);
}

SpiralSector::SpiralSector(const SpiralSector &right)
        : Component(right),
          planarArcGeometry_m(right.planarArcGeometry_m), dummy() {
    RefPartBunch_m = right.RefPartBunch_m;
    setElType(isDrift);
}

SpiralSector::~SpiralSector() {
}

ElementBase* SpiralSector::clone() const {
    return new SpiralSector(*this);
}


EMField &SpiralSector::getField() {
    return dummy;
}

const EMField &SpiralSector::getField() const {
    return dummy;
}

bool SpiralSector::apply(const size_t &i, const double &t,
                    Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

void SpiralSector::initialise(PartBunch *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
}

void SpiralSector::finalise() {
    RefPartBunch_m = NULL;
}

bool SpiralSector::bends() const {
    return true;
}

BGeometryBase& SpiralSector::getGeometry() {
    return planarArcGeometry_m;
}

const BGeometryBase& SpiralSector::getGeometry() const {
    return planarArcGeometry_m;
}

void SpiralSector::accept(BeamlineVisitor& visitor) const {
    visitor.visitSpiralSector(*this);
}


bool SpiralSector::apply(const Vector_t &R, const Vector_t &P,
                    const double &t, Vector_t &E, Vector_t &B) {
    Vector_t pos = R - centre_m;
    double radius = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    double phi = atan2(pos[2], pos[0]);
    double fringe = getFringe(phi);
    B[1] = fringe*Bz_m*pow(radius/r0_m, k_m);

    return true;
}

