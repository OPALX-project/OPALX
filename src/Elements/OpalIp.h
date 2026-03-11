//
// Class OpalIp
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
#ifndef OPAL_OpalIp_HH
#define OPAL_OpalIp_HH

#include "Elements/OpalElement.h"

class OpalIp: public OpalElement {

public:

    enum {
         GEOMETRY = COMMON,  // geometry of boundary, one more enum member besides the common ones in OpalElement.
         COLWINLEN,          // The lenght of the collision window
         SIZE
    };

    /// Exemplar constructor.
    OpalIp();

    virtual ~OpalIp();

    /// Make clone.
    virtual OpalIp* clone(const std::string& name);

    /// Test for drift.
    //  Return true.
    virtual bool isIp() const;

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalIp(const OpalIp&);
    void operator=(const OpalIp&);

    // Clone constructor.
    OpalIp(const std::string& name, OpalIp* parent);

};

#endif // OPAL_OpalIp_HH
