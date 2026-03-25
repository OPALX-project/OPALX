//
// Class OpalBeamBeam
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
#ifndef OPAL_OpalBeamBeam_HH
#define OPAL_OpalBeamBeam_HH

#include "Elements/OpalElement.h"

class OpalBeamBeam: public OpalElement {

public:

    enum {
         GEOMETRY = COMMON,  // geometry of boundary, one more enum member besides the common ones in OpalElement.
         COPY,               // Enable mirrored-bunch copy model
         VISUALIZE,          // Enable ASCII beam-beam-window visualization
         SIZE
    };

    /// Exemplar constructor.
    OpalBeamBeam();

    virtual ~OpalBeamBeam();

    /// Make clone.
    virtual OpalBeamBeam* clone(const std::string& name);

    /// Test for drift.
    //  Return true.
    virtual bool isBeamBeam() const;

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalBeamBeam(const OpalBeamBeam&);
    void operator=(const OpalBeamBeam&);

    // Clone constructor.
    OpalBeamBeam(const std::string& name, OpalBeamBeam* parent);

};

#endif // OPAL_OpalBeamBeam_HH
