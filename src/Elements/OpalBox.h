//
// Class OpalBox
//   The BOX element.
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
#ifndef OPAL_OpalBox_HH
#define OPAL_OpalBox_HH

#include "Elements/OpalElement.h"

class OpalBox : public OpalElement {
public:
    /// The attributes of class OpalBox.
    enum {
        WIDTH = COMMON,  // Full width along local x [m].
        HEIGHT,          // Full height along local y [m].
        SIZE
    };

    /// Exemplar constructor.
    OpalBox();

    virtual ~OpalBox();

    /// Make clone.
    virtual OpalBox* clone(const std::string& name);

    /**
     * @brief Update the embedded OPALX box.
     *
     * @throw OpalException if WIDTH, HEIGHT or L is not > 0, or if APERTURE is
     *        given. A box has no aperture: it removes particles inside its body.
     */
    virtual void update();

private:
    // Not implemented.
    OpalBox(const OpalBox&);
    void operator=(const OpalBox&);

    // Clone constructor.
    OpalBox(const std::string& name, OpalBox* parent);
};

#endif  // OPAL_OpalBox_HH
