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
#include "Processes/LocalProcesses/OpalTestProcess.h"

#include "Processes/LocalProcesses/TestCounterProcess.h"

OpalTestProcess::OpalTestProcess()
    : OpalLocalProcess(
              SIZE, "TESTPROCESS",
              "The \"TESTPROCESS\" statement defines a physics-free local process that "
              "counts particles inside the element it is attached to.") {
    registerOwnership(AttributeHandler::STATEMENT);
}

OpalTestProcess::OpalTestProcess(const std::string& name, OpalTestProcess* parent)
    : OpalLocalProcess(name, parent) {}

OpalTestProcess* OpalTestProcess::clone(const std::string& name) {
    return new OpalTestProcess(name, this);
}

std::shared_ptr<LocalProcess> OpalTestProcess::makeProcess() {
    return std::make_shared<TestCounterProcess>(getOpalName());
}
