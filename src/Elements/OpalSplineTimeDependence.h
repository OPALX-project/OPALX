//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_OPALSPLINETIMEDEPENDENCE_HH
#define OPAL_OPALSPLINETIMEDEPENDENCE_HH

#include "Elements/OpalElement.h"

/** OpalSplineTimeDependence provides UI wrapper for the
 *  SplineTimeDependence
 */

class OpalSplineTimeDependence : public OpalElement {
public:
    /** Enumeration maps to UI parameters */
    enum {
        ORDER = COMMON,
        TIMES,
        VALUES,
        SIZE  // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalSplineTimeDependence();

    /** No memory allocated so does nothing */
    ~OpalSplineTimeDependence() override = default;

    /** Inherited copy constructor */
    OpalSplineTimeDependence* clone(const std::string& name) override;

    /** Receive parameters from the parser and hand them off to the
     *  SplineTimeDependence
     */
    void update() override;

    /** Calls print on the OpalElement */
    void print(std::ostream&) const override;

    // Not implemented.
    OpalSplineTimeDependence(const OpalSplineTimeDependence&) = delete;
    void operator=(const OpalSplineTimeDependence&)           = delete;

private:
    // Clone constructor.
    OpalSplineTimeDependence(const std::string& name, OpalSplineTimeDependence* parent);

    static const std::string doc_string;
};

#endif  // OPAL_OPALSPLINETIMEDEPENDENCE_HH
