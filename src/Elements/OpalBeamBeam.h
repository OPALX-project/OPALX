//
// Class OpalBeamBeam
//   Input-language definition of collective beam-beam elements.
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

/**
 * @brief Parse BEAMBEAM attributes and update the placed BeamBeamRep.
 *
 * COPY_TIME is a simulation time in seconds; zero disables the virtual mirrored
 * primary. BBRIGID suppresses only the physical primary's collective kick, not
 * external fields or witness gathering. WITNESS_CONTAINERS selects passive
 * secondary containers. VISUALIZE enables optional rank-zero ASCII output.
 * See BeamBeam and BEAMBEAM::Config for scope and numerical-window semantics.
 */
class OpalBeamBeam : public OpalElement {
public:
    enum {
        GEOMETRY = COMMON,  // geometry of boundary, one more enum member besides the common ones in
                            // OpalElement.
        COPY_TIME,          // Start mirrored-bunch copy model at or after this simulation time [s]
        BBRIGID,            // Suppress the BeamBeam collective kick on the physical source
        VISUALIZE,          // Enable ASCII beam-beam-window visualization
        WITNESS_CONTAINERS,  // Passive containers that gather source BeamBeam fields
        SIZE
    };

    /// Exemplar constructor.
    OpalBeamBeam();

    virtual ~OpalBeamBeam();

    /// Make clone.
    /// @param name Name for the cloned input definition.
    /// @return Caller-owned definition with its own BeamBeamRep.
    OpalBeamBeam* clone(const std::string& name) override;

    /// @return True: this input definition represents a BEAMBEAM element.
    virtual bool isBeamBeam() const;

    /// Validate input attributes and update the embedded BeamBeamRep.
    void update() override;

private:
    // Not implemented.
    OpalBeamBeam(const OpalBeamBeam&);
    void operator=(const OpalBeamBeam&);

    // Clone constructor.
    OpalBeamBeam(const std::string& name, OpalBeamBeam* parent);
};

#endif  // OPAL_OpalBeamBeam_HH
