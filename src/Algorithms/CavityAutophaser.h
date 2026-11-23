//
// Class CavityAutophaser
//
// This class determines the phase of an RF cavity for which the reference particle
// is accelerated to the highest energy.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#ifndef CAVITYAUTOPHASER
#define CAVITYAUTOPHASER

#include <string>
#include <utility>
#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/PartData.h"

/// Name and phase (rad) of each cavity that autophase has set, in the order they were set.
using CavityPhases = std::vector<std::pair<std::string, double>>;

class CavityAutophaser {
public:
    /**
     * @param numRefinements Halvings of the phase step after the coarse search, the value of
     *        OPTION AUTOPHASE.
     * @param phases If given, getPhaseAtMaxEnergy() adds the phase it sets on the cavity; a
     *        checkpoint stores this list so that a restart keeps the phases.
     */
    CavityAutophaser(
            const PartData& ref, std::shared_ptr<ElementBase> cavity, int numRefinements,
            CavityPhases* phases = nullptr);

    ~CavityAutophaser();

    double getPhaseAtMaxEnergy(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, double t, double dt);

private:
    double guessCavityPhase(double t);
    std::pair<double, double> optimizeCavityPhase(double initialGuess, double t, double dt);

    double track(double t, const double dt, const double phase, std::ofstream* out = nullptr) const;

    const PartData& itsReference_m;
    std::shared_ptr<ElementBase> itsCavity_m;
    int numRefinements_m;
    CavityPhases* phases_m;

    Vector_t<double, 3> initialR_m;
    Vector_t<double, 3> initialP_m;
};

#endif
