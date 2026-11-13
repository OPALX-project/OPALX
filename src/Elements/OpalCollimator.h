//
// Class OpalCollimator
//   The COLLIMATOR element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalCollimator_HH
#define OPAL_OpalCollimator_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalCollimator : public OpalElement {
public:
    /// Exemplar constructor.
    OpalCollimator();

    virtual ~OpalCollimator();

    /// Make clone.
    virtual OpalCollimator* clone(const std::string& name);

    /**
     * @brief Build the runtime collimator of one occurrence from its attribute values.
     *
     * @throw OpalException if the input gives no APERTURE or a length L <= 0.
     *        Scraping only happens for z in [0, L) outside the aperture, so a
     *        collimator missing either would silently transmit everything.
     */
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    // Not implemented.
    OpalCollimator(const OpalCollimator&);
    void operator=(const OpalCollimator&);

    // Clone constructor.
    OpalCollimator(const std::string& name, OpalCollimator* parent);
};

#endif  // OPAL_OpalCollimator_HH
