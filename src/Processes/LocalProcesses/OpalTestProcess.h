//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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
#ifndef OPAL_OPALTESTPROCESS_H
#define OPAL_OPALTESTPROCESS_H

#include "Processes/LocalProcesses/OpalLocalProcess.h"

/**
 * @brief The "TESTPROCESS" definition: a physics-free local process for testing.
 *
 * Creates a TestCounterProcess that counts particles inside the element body each step
 * without modifying particle state. Intended for validating the local-process plumbing.
 */
class OpalTestProcess : public OpalLocalProcess {
public:
    /// The attributes of class OpalTestProcess (none beyond the base class).
    enum { SIZE };

    /// Exemplar constructor.
    OpalTestProcess();

    virtual OpalTestProcess* clone(const std::string& name);

private:
    // Not implemented.
    OpalTestProcess(const OpalTestProcess&);
    void operator=(const OpalTestProcess&);

    /// Clone constructor.
    OpalTestProcess(const std::string& name, OpalTestProcess* parent);

    virtual std::shared_ptr<LocalProcess> makeProcess();
};

#endif
