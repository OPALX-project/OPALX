//
// Class SDDSColumnSet
//   This class writes rows of SDDS files.
//
// Copyright (c) 2019, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef SDDSWRITERCOLUMNSET_H
#define SDDSWRITERCOLUMNSET_H

#ifdef __CUDACC__
#pragma push_macro("__cpp_consteval")
#pragma push_macro("_NODISCARD")
#pragma push_macro("__builtin_LINE")

#define __cpp_consteval 201811L

#ifdef _NODISCARD
    #undef _NODISCARD
    #define _NODISCARD
#endif

#define consteval constexpr

#include <source_location>

#undef consteval
#pragma pop_macro("__cpp_consteval")
#pragma pop_macro("_NODISCARD")
#else
#include <source_location>
#endif

#include "Structure/SDDSColumn.h"
#include "Utilities/OpalException.h"

#include <string>
#include <ostream>
#include <vector>
#include <map>

class SDDSColumnSet {
public:
    SDDSColumnSet();

    void addColumn(const std::string& name,
                   const std::string& type,
                   const std::string& unit,
                   const std::string& desc,
                   std::ios_base::fmtflags flags = std::ios_base::scientific,
                   unsigned short precision = 15);

    template<typename T>
    void addColumnValue(const std::string& name,
                        const T& val);

    void writeHeader(std::ostream& os,
                     const std::string& indent) const;

    void writeRow(std::ostream& os) const;

    bool hasColumns() const;

private:
    std::vector<SDDSColumn> columns_m;
    std::map<std::string, size_t> name2idx_m;
};


inline
SDDSColumnSet::SDDSColumnSet()
{ }


template<typename T>
void SDDSColumnSet::addColumnValue(const std::string& name,
                                   const T& val) {

    auto it = name2idx_m.find(name);
    if (it == name2idx_m.end()) {
        throw OpalException("SDDSColumnSet::addColumnValue",
                            "column name '" + name + "' doesn't exists");
    }

    auto & col = columns_m[it->second];
    col.addValue(val);
}


inline
bool SDDSColumnSet::hasColumns() const {
    return !name2idx_m.empty();
}


#endif
