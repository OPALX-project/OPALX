/* 
 *  Copyright (c) 2012, Chris Rogers
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

#include "Utilities/OpalRingSection.h"

OpalRingSection::OpalRingSection() 
  : component_m(NULL),
    componentPosition_m(0.), componentOrientation_m(0.),
    startPosition_m(0.), startOrientation_m(0.),
    endPosition_m(0.), endOrientation_m(0.) {
}

bool OpalRingSection::isOnOrPastStartPlane(const Vector_t& pos) {
    Vector_t pos_transformed = pos-startPosition_m;
    double dot_prod = pos_transformed(0)*startOrientation_m(0)+
                      pos_transformed(1)*startOrientation_m(1)+
                      pos_transformed(2)*startOrientation_m(2);
    return dot_prod >= 0.;
}

bool OpalRingSection::isPastEndPlane(const Vector_t& pos) {
    Vector_t pos_transformed = pos-endPosition_m;
    double dot_prod = pos_transformed(0)*endOrientation_m(0)+
                      pos_transformed(1)*endOrientation_m(1)+
                      pos_transformed(2)*endOrientation_m(2);
    return dot_prod > 0.;
}

bool OpalRingSection::getFieldValue(const Vector_t& pos,
                                    const Vector_t& centroid, const double& t,
                                    Vector_t& E, Vector_t& B) {
    // transform position into local coordinate system
    Vector_t pos_local = pos-componentPosition_m;
    rotate(pos_local);

    bool outOfBounds = component_m->apply(pos_local, centroid, t, E, B);
    if (outOfBounds) {
        return true;
    }
    // rotate fields back to global coordinate system
    rotate_back(E);
    rotate_back(B);
    B *= 10.;
    return false;
}

void OpalRingSection::updateComponentOrientation() {
    sin0_m = sin(componentOrientation_m(0));
    cos0_m = cos(componentOrientation_m(0));
    sin1_m = sin(componentOrientation_m(1));
    cos1_m = cos(componentOrientation_m(1));
    sin2_m = sin(componentOrientation_m(2));
    cos2_m = cos(componentOrientation_m(2));
}

