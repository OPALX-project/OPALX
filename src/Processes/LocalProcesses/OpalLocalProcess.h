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
#ifndef OPAL_OPALLOCALPROCESS_H
#define OPAL_OPALLOCALPROCESS_H

#include "AbstractObjects/Definition.h"

#include <memory>
#include <string>

class LocalProcess;

/**
 * @brief Abstract base for input-file definitions of local physics processes.
 *
 * Each concrete process type provides its own definition command (e.g. "TESTPROCESS";
 * later "SYNCHROTRONRADIATION", ...). A named instance is defined once in the input file
 * and attached to beamline elements by name via the PROCESSES element attribute:
 *
 * @code
 * TP1: TESTPROCESS;
 * B1:  SBEND, L=1.5, ANGLE=0.2, ELEMEDGE=0.5, PROCESSES={"TP1"};
 * @endcode
 *
 * @note The definition must precede the elements that reference it in the input file.
 *   Redefining an existing process name is an error.
 */
class OpalLocalProcess : public Definition {
public:
    virtual ~OpalLocalProcess();

    /**
     * @brief Find a named local-process definition in OpalData.
     *
     * @throw OpalException if no local-process definition with this name exists.
     */
    static OpalLocalProcess* find(const std::string& name);

    /**
     * @brief Return the runtime process instance.
     * 
     * @note The object is created the first time getProcess() is called. Repeat
     * calls return the same object.  
     *
     */
    std::shared_ptr<LocalProcess> getProcess();

    virtual void execute();

protected:
    /// Exemplar constructor.
    OpalLocalProcess(int size, const char* name, const char* help);

    /// Clone constructor.
    OpalLocalProcess(const std::string& name, OpalLocalProcess* parent);

    /// Factory hook implemented by concrete process definitions.
    virtual std::shared_ptr<LocalProcess> makeProcess() = 0;

private:
    // Not implemented.
    OpalLocalProcess();
    OpalLocalProcess(const OpalLocalProcess&);
    void operator=(const OpalLocalProcess&);

    /// Runtime process instance
    std::shared_ptr<LocalProcess> process_m;
};

#endif
