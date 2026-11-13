//
// Class OpalVariableRFCavity
//   The class provides the user interface for the VARIABLE_RF_CAVITY object.
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
#ifndef OPAL_OPALVARIABLERFCAVITY_H
#define OPAL_OPALVARIABLERFCAVITY_H

#include <memory>

#include "Elements/OpalElement.h"

class OpalVariableRFCavity : public OpalElement {
public:
    /** enum maps string to integer value for UI definitions */
    enum {
        PHASE_MODEL = COMMON,
        AMPLITUDE_MODEL,
        FREQUENCY_MODEL,
        WIDTH,
        HEIGHT,
        SIZE  // size of the enum
    };

    /** Copy constructor **/
    OpalVariableRFCavity(const std::string& name, OpalVariableRFCavity* parent);

    /** Default constructor **/
    OpalVariableRFCavity();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
     **/
    OpalVariableRFCavity* clone();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
     */
    OpalVariableRFCavity* clone(const std::string& name) override;

    /** Destructor does nothing */
    ~OpalVariableRFCavity() override = default;

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

    // Not implemented.
    OpalVariableRFCavity(const OpalVariableRFCavity&) = delete;
    void operator=(const OpalVariableRFCavity&)       = delete;

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    static const std::string doc_string;
};

#endif  // OPAL_OPALVARIABLERFCAVITY_H
