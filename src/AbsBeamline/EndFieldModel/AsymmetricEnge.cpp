/*
 *  Copyright (c) 2017-2026, Chris Rogers
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

#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"

namespace endfieldmodel {

    AsymmetricEnge::AsymmetricEnge(
            const std::vector<double> aStart, double x0Start, double lambdaStart,
            const std::vector<double> aEnd, double x0End, double lambdaEnd) {
        if (aStart.size()) {
            config_m.engeStart_m.a_m = Enge::makeView(aStart, "AsymmetricEngeStart");
        }
        config_m.engeStart_m.x0_m     = x0Start;
        config_m.engeStart_m.lambda_m = lambdaStart;
        // x0 is held in this
        if (aEnd.size()) {
            config_m.engeEnd_m.a_m = Enge::makeView(aEnd, "AsymmetricEngeEnd");
        }
        config_m.engeEnd_m.x0_m     = x0End;
        config_m.engeEnd_m.lambda_m = lambdaEnd;
    }

    void AsymmetricEnge::rescale(double scaleFactor) {
        config_m.engeStart_m.x0_m *= scaleFactor;
        config_m.engeStart_m.lambda_m *= scaleFactor;
        config_m.engeEnd_m.x0_m *= scaleFactor;
        config_m.engeEnd_m.lambda_m *= scaleFactor;
    }

    std::ostream& AsymmetricEnge::print(std::ostream& out) const {
        out << "AsymmetricEnge start ";
        auto conf = config_m.engeStart_m;
        out << "AsymmetricEnge start function l=" << conf.lambda_m << " x0=" << conf.x0_m << " c=";
        for (size_t ai = 0; ai < conf.a_m.extent(0); ++ai) {
            out << conf.a_m(ai) << " ";
        }
        conf = config_m.engeStart_m;
        out << "               end   function l=" << conf.lambda_m << " x0=" << conf.x0_m << " c=";
        for (size_t ai = 0; ai < conf.a_m.extent(0); ++ai) {
            out << conf.a_m(ai) << " ";
        }
        return out;
    }

    void AsymmetricEnge::initialise() {
        Enge::setEngeDiffIndices(EngeConfig::max_derivative, config_m.engeStart_m);
        Enge::setEngeDiffIndices(EngeConfig::max_derivative, config_m.engeEnd_m);
    }

    void AsymmetricEnge::finalise() {
        Enge::releaseEngeConfig(config_m.engeStart_m);
        Enge::releaseEngeConfig(config_m.engeEnd_m);
    }

}  // namespace endfieldmodel
