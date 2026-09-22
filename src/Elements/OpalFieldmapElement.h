//
// Class OpalFieldmapElement
//   The FIELDMAP element.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_OpalFieldmapElement_HH
#define OPALX_OpalFieldmapElement_HH

#include "Elements/OpalElement.h"

/**
 * @class OpalFieldmapElement
 * @brief The FIELDMAP element: a straight element whose field is a tabulated map.
 *
 * The element takes its box from the map's own extent, so neither L nor ELEMEDGE may be
 * given: ELEMEDGE positions an element by path length along the reference orbit, which needs
 * a length along that orbit, and a map's extent is not one. Placement is the absolute 6D lab
 * pose (X, Y, Z, THETA, PHI, PSI).
 */
class OpalFieldmapElement : public OpalElement {
public:
    /// The attributes of class OpalFieldmapElement.
    enum {
        FMAPFN = COMMON,  // The field map filename.
        SCALE,            // Plain multiplier on the tabulated magnetic field.
        ESCALE,           // Plain multiplier on the tabulated electric field.
        ZREVERSE,         // Read the field map back to front.
        SIZE
    };

    /// Exemplar constructor.
    OpalFieldmapElement();

    virtual ~OpalFieldmapElement();

    /// Make clone.
    virtual OpalFieldmapElement* clone(const std::string& name) override;

    /// Update the embedded OPALX element.
    virtual void update() override;

private:
    // Not implemented.
    OpalFieldmapElement(const OpalFieldmapElement&);
    void operator=(const OpalFieldmapElement&);

    // Clone constructor.
    OpalFieldmapElement(const std::string& name, OpalFieldmapElement* parent);
};

#endif  // OPALX_OpalFieldmapElement_HH
