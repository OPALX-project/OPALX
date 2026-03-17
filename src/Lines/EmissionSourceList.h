//
// Class EmissionSourceList
//   The EMISSIONSOURCELIST definition. A named list of EmissionSource
//   elements, buildable like a Line: name: EMISSIONSOURCELIST = (ES1, ES2, ...).
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_EMISSIONSOURCELIST_HH
#define OPAL_EMISSIONSOURCELIST_HH

#include "AbstractObjects/Definition.h"
#include "Structure/EmissionSource.h"

#include <vector>

class Statement;

class EmissionSourceList : public Definition {
public:
    /// Exemplar constructor.
    EmissionSourceList();

    virtual ~EmissionSourceList();

    virtual bool canReplaceBy(Object* object);
    virtual EmissionSourceList* clone(const std::string& name);
    virtual void execute();

    virtual void parse(Statement& stat);

    static EmissionSourceList* find(const std::string& name);

    /// Return the list of EmissionSource pointers. Valid after parse/execute.
    const std::vector<EmissionSource*>& fetchSources() const { return sources_m; }

private:
    EmissionSourceList(const EmissionSourceList&);
    void operator=(const EmissionSourceList&);

    EmissionSourceList(const std::string& name, EmissionSourceList* parent);

    void parseList(Statement& stat);

    std::vector<EmissionSource*> sources_m;
};

#endif  // OPAL_EMISSIONSOURCELIST_HH
