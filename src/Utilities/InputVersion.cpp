#include "Utilities/InputVersion.h"

#include "BuildInfo.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <array>
#include <cctype>

namespace InputVersion {

    namespace {
        constexpr Version buildVersion{buildinfo::version_major, buildinfo::version_minor};

        // Append new entries at the end and raise the array size. An entry names the first
        // release with the new meaning and says how to adapt an older input. A hypothetical
        // change in OPALX 1.1 would read:
        //
        //     constexpr std::array<Change, 1> changeTable{
        //             Change{{1, 1},
        //                    "LAG of RFCAVITY is given in degrees instead of radians. "
        //                    "Multiply the values of OPALX 1.0 inputs by 180/PI."}};
        //
        // Inputs with VERSION = "1.0" would then stop at TRACK or COF and print this text.
        // No OPALX release has changed the meaning of an existing input yet.
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

    void enforce() {
        static std::optional<std::string> checked;
        if (checked == Options::version) {
            return;
        }
        checked = Options::version;

        std::optional<Version> input;
        if (!Options::version.empty()) {
            input = parse(Options::version);
        }
        const std::string current = toString(buildVersion);
        const Result result       = evaluate(input, buildVersion, changeTable);

        switch (result.status) {
            case Status::NOTGIVEN: {
                Inform warning("Warning");
                warning << "The input does not state the OPALX release it was written for and "
                        << "is read as an OPALX " << current << " input.\n"
                        << "Add OPTION, VERSION = \"" << current << "\"; so that input changes "
                        << "of later releases are reported." << endl;
                break;
            }
            case Status::NEWER: {
                Inform warning("Warning");
                warning << "The input was written for OPALX " << toString(*input)
                        << ", but this is OPALX " << current << ".\n"
                        << "Input changes made after OPALX " << current
                        << " are unknown to this executable." << endl;
                break;
            }
            case Status::OLDERAFFECTED: {
                Inform error("Error");
                error << "\nThe input was written for OPALX " << toString(*input)
                      << ". Later releases changed the meaning of existing input:\n"
                      << endl;
                for (const Change& change : result.changes) {
                    error << "* OPALX " << toString(change.since) << ": " << change.text << "\n"
                          << endl;
                }
                error << "Adapt the input to these changes, then set\n"
                      << "OPTION, VERSION = \"" << current << "\";" << endl;
                throw OpalException(
                        "InputVersion::enforce",
                        "The input was written for OPALX " + toString(*input)
                                + " and is affected by later input changes.");
            }
            case Status::CURRENT:
            case Status::OLDERUNAFFECTED:
                break;
        }
    }

}  // namespace InputVersion
