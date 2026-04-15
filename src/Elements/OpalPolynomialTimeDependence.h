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

#ifndef OPAL_OPALPOLYNOMIALTIMEDEPENDENCE_HH
#define OPAL_OPALPOLYNOMIALTIMEDEPENDENCE_HH

#include "Elements/OpalElement.h"

/** OpalPolynomialTimeDependence provides UI wrapper for the
 *  PolynomialTimeDependence
 */

class OpalPolynomialTimeDependence : public OpalElement {
public:
    /** Enumeration maps to UI parameters */
    enum {
        P0 = COMMON,
        P1,
        P2,
        P3,
        COEFFICIENTS,
        SIZE  // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalPolynomialTimeDependence();

    /** No memory allocated so does nothing */
    ~OpalPolynomialTimeDependence() override = default;

    /** Inherited copy constructor */
    OpalPolynomialTimeDependence* clone(const std::string& name) override;

    /** Receive parameters from the parser and hand them off to the
     *  PolynomialTimeDependence
     */
    void update() override;

    /** Calls print on the OpalElement */
    void print(std::ostream&) const override;

    // Not implemented.
    OpalPolynomialTimeDependence(const OpalPolynomialTimeDependence&) = delete;
    void operator=(const OpalPolynomialTimeDependence&)               = delete;

private:
    // Clone constructor.
    OpalPolynomialTimeDependence(const std::string& name, OpalPolynomialTimeDependence* parent);

    static const std::string doc_string;
};

#endif  // OPAL_OPALPOLYNOMIALTIMEDEPENDENCE_HH
