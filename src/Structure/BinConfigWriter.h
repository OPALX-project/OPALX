//
// Class BinConfigWriter
//   Writes binning configuration snapshots to a JSON file.
//
// Copyright (c) 2026, Paul Scherrer Institut
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

#ifndef OPAL_BIN_CONFIG_WRITER_H
#define OPAL_BIN_CONFIG_WRITER_H

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Helper to write binning configuration snapshots as a JSON array.
 *
 * The writer owns a single JSON file. On construction it opens the file in
 * truncate mode and emits the opening '['. Each call to writeEntry appends a
 * JSON object (with comma handling). The destructor closes the array and file.
 */
class BinConfigWriter {
public:
    explicit BinConfigWriter(const std::string& fileName);

    ~BinConfigWriter();

    void writeEntry(
            long long step, double time, bool preMerge, const std::vector<std::size_t>& binCounts,
            const std::vector<double>& binWidths, double xMin);

private:
    std::ofstream os_;
    bool first_;
};

#endif  // OPAL_BIN_CONFIG_WRITER_H
