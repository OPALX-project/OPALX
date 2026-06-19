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
#include "Processes/LocalProcesses/OpalLocalProcess.h"

#include "AbstractObjects/OpalData.h"
#include "Processes/LocalProcesses/LocalProcess.h"
#include "Utilities/OpalException.h"

OpalLocalProcess::OpalLocalProcess(int size, const char* name, const char* help)
    : Definition(size, name, help) {}

OpalLocalProcess::OpalLocalProcess(const std::string& name, OpalLocalProcess* parent)
    : Definition(name, parent) {}

OpalLocalProcess::~OpalLocalProcess() {}

OpalLocalProcess* OpalLocalProcess::find(const std::string& name) {
    OpalLocalProcess* process =
            dynamic_cast<OpalLocalProcess*>(OpalData::getInstance()->find(name));

    if (process == nullptr) {
        throw OpalException(
                "OpalLocalProcess::find()",
                "Local process \"" + name
                        + "\" not found. The process definition must precede the "
                          "element that references it in the input file.");
    }
    return process;
}

std::shared_ptr<LocalProcess> OpalLocalProcess::getProcess() {
    if (!process_m) {
        process_m = makeProcess();
    }
    return process_m;
}

void OpalLocalProcess::execute() {}
