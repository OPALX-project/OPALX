// Copyright (c) 2023, Chris Rogers
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

#ifndef OPAL_OpalOutputPlane_H
#define OPAL_OpalOutputPlane_H

#include "Elements/OpalElement.h"

// Class OpalOutputPlane
// ------------------------------------------------------------------------
/// Interface for output plane element.


class OpalOutputPlane final : public OpalElement {

public:

    // The attributes of class OpalOutputPlane.
    enum {
        CENTRE = COMMON, // plane centre position (3-vector)
        NORMAL,          // plane normal (3-vector)
        RADIUS,          // maximum allowed radial displacement from centre

        SIZE
    };

    // Exemplar constructor.
    OpalOutputPlane();
    ~OpalOutputPlane() override;

    // Make clone.
    OpalOutputPlane *clone(const std::string &name) override;

    /** Update with new parameters from UI parser */
    void update() override;

    // Not implemented.
    OpalOutputPlane(const OpalOutputPlane &) = delete;
    void operator=(const OpalOutputPlane &) = delete;

private:
    // Clone constructor.
    OpalOutputPlane(const std::string &name, OpalOutputPlane *parent);
};

#endif // OPAL_OpalOutputPlane_H
