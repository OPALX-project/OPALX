#ifndef OPALX_INPUT_VERSION_H
#define OPALX_INPUT_VERSION_H

#include <compare>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Compatibility check between an input file and the running OPALX release.
 *
 * An input states the release it was written for with `OPTION, VERSION = "major.minor";`.
 * The release number alone never stops a run. A run stops only if the change table holds an
 * entry newer than the stated release, that is, if a later release gave an existing input a
 * different meaning.
 *
 * | Input version                        | Outcome                         |
 * |--------------------------------------|---------------------------------|
 * | not given                            | warning, read as current release |
 * | current release                      | run                             |
 * | older, no newer table entry          | run                             |
 * | older, newer table entries           | abort and list the entries      |
 * | newer than the executable            | warning, run                    |
 *
 * Patch releases must not change the meaning of an input, so the patch number is ignored.
 */
namespace InputVersion {

    struct Version {
        int major = 0;
        int minor = 0;

        constexpr auto operator<=>(const Version&) const = default;
    };

    /**
     * @brief One release change after which an existing input still parses but means
     * something else (unit, default value, sign or frame convention).
     *
     * Renamed or removed keywords need no entry, because the parser rejects them.
     */
    struct Change {
        Version since;          ///< First release with the new meaning.
        std::string_view text;  ///< What changed and how to adapt an input.
    };

    enum class Status { NOTGIVEN, CURRENT, OLDERUNAFFECTED, OLDERAFFECTED, NEWER };

    struct Result {
        Status status;
        std::vector<Change> changes;  ///< Entries newer than the input; filled for OLDERAFFECTED.
    };

    /**
     * @brief Checks that a change table is ordered by release and holds no entry newer than
     * the build.
     *
     * An entry ahead of the build would abort every input, because no accepted `VERSION`
     * could reach it. The change that adds an entry for the next release therefore also
     * raises the project version.
     */
    constexpr bool isValidTable(std::span<const Change> table, Version build) {
        for (std::size_t i = 0; i < table.size(); ++i) {
            if (table[i].since > build) {
                return false;
            }
            if (i > 0 && table[i].since < table[i - 1].since) {
                return false;
            }
        }
        return true;
    }

    /// Release of this executable, without the patch number.
    Version currentVersion();

    /// Input changes of all OPALX releases, ordered by release.
    std::span<const Change> registeredChanges();

    /**
     * @brief Parses `"major.minor"`; an optional `.patch` is accepted and ignored.
     * @throws OpalException for any other text, including the packed integers of OPAL
     * such as `10900`.
     */
    Version parse(const std::string& text);

    std::string toString(Version version);

    /// Applies the table in the namespace description to one input version.
    Result evaluate(
            const std::optional<Version>& input, Version build, std::span<const Change> table);

}  // namespace InputVersion

#endif
