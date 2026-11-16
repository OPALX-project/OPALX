#include "Utilities/InputVersion.h"

#include "BuildInfo.h"
#include "Utilities/OpalException.h"

#include <array>
#include <cctype>

namespace InputVersion {

    namespace {
        constexpr Version buildVersion{buildinfo::version_major, buildinfo::version_minor};

        // Append new entries at the end and raise the array size. No OPALX release has
        // changed the meaning of an existing input yet.
        constexpr std::array<Change, 0> changeTable{};

        static_assert(
                isValidTable(changeTable, buildVersion),
                "Input change table must be ordered by release and must not run ahead of the "
                "project version in CMakeLists.txt.");

        /// Reads a decimal number of at most four digits starting at pos and advances pos.
        bool readNumber(const std::string& text, std::size_t& pos, int& number) {
            const std::size_t begin = pos;
            number                  = 0;
            while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) {
                number = 10 * number + (text[pos] - '0');
                ++pos;
            }
            return pos > begin && pos - begin <= 4;
        }
    }  // namespace

    Version currentVersion() { return buildVersion; }

    std::span<const Change> registeredChanges() { return changeTable; }

    Version parse(const std::string& text) {
        Version version;
        int patch       = 0;
        std::size_t pos = 0;
        bool valid = readNumber(text, pos, version.major) && pos < text.size() && text[pos++] == '.'
                     && readNumber(text, pos, version.minor);
        if (valid && pos < text.size()) {
            valid = text[pos++] == '.' && readNumber(text, pos, patch) && pos == text.size();
        }
        if (!valid) {
            throw OpalException(
                    "InputVersion::parse",
                    "VERSION = \"" + text + "\" is not an OPALX release. Write it as "
                            + "\"major.minor\", for example OPTION, VERSION = \""
                            + toString(buildVersion) + "\";");
        }
        return version;
    }

    std::string toString(Version version) {
        return std::to_string(version.major) + "." + std::to_string(version.minor);
    }

    Result evaluate(
            const std::optional<Version>& input, Version build, std::span<const Change> table) {
        if (!input) {
            return {Status::NOTGIVEN, {}};
        }
        if (*input == build) {
            return {Status::CURRENT, {}};
        }
        if (*input > build) {
            return {Status::NEWER, {}};
        }

        Result result{Status::OLDERUNAFFECTED, {}};
        for (const Change& change : table) {
            if (change.since > *input) {
                result.changes.push_back(change);
            }
        }
        if (!result.changes.empty()) {
            result.status = Status::OLDERAFFECTED;
        }
        return result;
    }

}  // namespace InputVersion
