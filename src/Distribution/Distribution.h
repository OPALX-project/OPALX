//
// Class Distribution
//   This class defines the initial beam that is injected or emitted into the simulation.
//
// Copyright (c) 2008 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

#include <string>

#include "AbstractObjects/Definition.h"
#include "Distribution/DistributionSettings.h"

/**
 * @brief The DISTRIBUTION command: the attributes of one particle distribution.
 *
 * A RUN reads the attribute values from the run description and works out the sampler values
 * with DistributionSettings.
 */
class Distribution : public Definition {
public:
    Distribution();

    virtual ~Distribution();

    virtual bool canReplaceBy(Object* object);

    /// Check the attribute values (DistributionSettings) when the statement is read.
    virtual void execute();

    static Distribution* find(const std::string& name);

private:
    // Not implemented.
    Distribution(const Distribution&)   = delete;
    void operator=(const Distribution&) = delete;
};

#endif  // OPAL_Distribution_HH
