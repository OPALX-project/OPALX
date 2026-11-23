//
// Class DumpEMFields
//   The DUMPEMFIELDS command: a request to write the fields of the lattice on a
//   user-defined Cartesian grid.
//
// Copyright (c) 2017, Chris Rogers
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
#ifndef OPAL_BASICACTIONS_DUMPEMFIELDS_HH
#define OPAL_BASICACTIONS_DUMPEMFIELDS_HH

#include <string>
#include <vector>
#include "AbstractObjects/Action.h"
#include "Track/RunDescription.h"

/**
 * @brief The DUMPEMFIELDS command: a request to write the external fields of the tracked line
 *        on a 4D Cartesian grid in space and time.
 *
 * execute() checks the grid and registers the request. freeze() copies the registered requests
 * into the description of each later RUN, and the tracker writes them with writeEMFieldDumps()
 * when tracking starts.
 */
class DumpEMFields : public Action {
public:
    /// The common attributes of DumpEMFields.
    enum {
        FILE_NAME,
        COORDINATE_SYSTEM,
        X_START,
        DX,
        X_STEPS,
        Y_START,
        DY,
        Y_STEPS,
        Z_START,
        DZ,
        Z_STEPS,
        T_START,
        DT,
        T_STEPS,
        SIZE
    };

    DumpEMFields();

    ~DumpEMFields() override = default;

    /**
     * @brief Register the request with the attribute values it has now.
     * @throws OpalException if X_STEPS, Y_STEPS, Z_STEPS or T_STEPS is not a positive integer.
     */
    void execute() override;

    /// The requests of the DUMPEMFIELDS commands executed so far, in input order.
    static const std::vector<opalx::run::EMFieldDumpDescription>& requests();

    /** Print the attributes of DumpEMFields to standard out */
    void print(std::ostream& os) const override;

    /* For use only by unit tests, clear the registered requests */
    static void clearDumps() { requests_m.clear(); }

private:
    opalx::run::EMFieldDumpDescription describe() const;
    static void checkInt(double value, const std::string& name, double tolerance = 1e-9);

    static std::vector<opalx::run::EMFieldDumpDescription> requests_m;

    DumpEMFields(const DumpEMFields& dump)            = delete;
    DumpEMFields& operator=(const DumpEMFields& dump) = delete;
};

inline std::ostream& operator<<(std::ostream& os, const DumpEMFields& b) {
    b.print(os);
    return os;
}

#endif  // ifdef OPAL_BASICACTIONS_DUMPEMFIELDS_HH
