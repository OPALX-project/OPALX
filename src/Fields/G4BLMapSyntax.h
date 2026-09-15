//
// G4BLMapSyntax
//   Line-level helpers shared by the G4beamline field map readers.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPALX_G4BLMapSyntax_HH
#define OPALX_G4BLMapSyntax_HH

#include <istream>
#include <map>
#include <sstream>
#include <string>

/**
 * @namespace G4BLMapSyntax
 * @brief Reading the line syntax common to every G4beamline field map file.
 *
 * Both the `cylinder` and the `grid` formats open with an optional `param` line followed by
 * a section keyword carrying `key=value` pairs, and both allow '#' comments and blank lines
 * anywhere. These three functions are all that is shared; the grid layout and the
 * interpolation differ completely between the two.
 */
namespace G4BLMapSyntax {

    /**
     * @brief Read the next line that carries content.
     *
     * Strips '#' comments and surrounding whitespace, and skips lines left empty.
     *
     * @note This deliberately does not use Fieldmap::getLine(). That reads through a 256
     *       character buffer (READ_BUFFER_LENGTH), while a cylinder data row is nR values
     *       wide -- about 660 characters for the muE4 WSX map -- and would be silently
     *       truncated.
     *
     * @return false at end of file.
     */
    inline bool nextLine(std::istream& in, std::string& line) {
        while (std::getline(in, line)) {
            const size_t comment = line.find('#');
            if (comment != std::string::npos) {
                line.erase(comment);
            }
            const size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) {
                continue;
            }
            const size_t last = line.find_last_not_of(" \t\r\n");
            line              = line.substr(first, last - first + 1);
            return true;
        }
        return false;
    }

    /// @brief First whitespace separated token of a line, the section keyword.
    inline std::string firstToken(const std::string& line) {
        std::istringstream is(line);
        std::string token;
        is >> token;
        return token;
    }

    /// @brief The `key=value` pairs of a line, skipping the leading keyword.
    inline std::map<std::string, std::string> keyValues(const std::string& line) {
        std::istringstream is(line);
        std::string token;
        is >> token;  // discard the section keyword

        std::map<std::string, std::string> pairs;
        while (is >> token) {
            const size_t equals = token.find('=');
            if (equals != std::string::npos) {
                pairs[token.substr(0, equals)] = token.substr(equals + 1);
            }
        }
        return pairs;
    }

}  // namespace G4BLMapSyntax

#endif  // OPALX_G4BLMapSyntax_HH
