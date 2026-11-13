//
// Class OpalDrift
//   The class of OPAL drift spaces.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalDrift_HH
#define OPAL_OpalDrift_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalDrift : public OpalElement {
public:
    enum {
        GEOMETRY = COMMON,  // geometry of boundary, one more enum member besides the common ones in
                            // OpalElement.
        NSLICES,            // The number of slices / steps per element for map tracking
        SIZE
    };

    /// Exemplar constructor.
    OpalDrift();

    virtual ~OpalDrift();

    /// Make clone.
    virtual OpalDrift* clone(const std::string& name);

    /// Test for drift.
    //  Return true.
    virtual bool isDrift() const;

    /// Update the embedded OPALX drift.
    virtual void update();

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    // Not implemented.
    OpalDrift(const OpalDrift&);
    void operator=(const OpalDrift&);

    // Clone constructor.
    OpalDrift(const std::string& name, OpalDrift* parent);
};

#endif  // OPAL_OpalDrift_HH
