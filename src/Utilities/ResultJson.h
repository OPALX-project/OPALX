#ifndef OPAL_RESULT_JSON_H
#define OPAL_RESULT_JSON_H
#include <unistd.h>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "Utilities/OpalException.h"

/// Host optics JSON and atomic replacement shared by COF and MATCH.
namespace opalx::report {
    // Small streaming writer for the optics exports: JSON escaping, finite numbers,
    // classic-locale round-trip precision and two-space indentation.
    class ResultJsonWriter {
        std::ostringstream out;
        std::vector<bool> first;

        void quoted(const std::string& text) {
            constexpr char hex[] = "0123456789abcdef";
            out << '"';
            for (unsigned char c : text) {
                if (c == '"' || c == '\\')
                    out << '\\' << c;
                else if (c < 0x20)
                    out << "\\u00" << hex[c >> 4] << hex[c & 15];
                else
                    out << c;
            }
            out << '"';
        }
        void prefix(const char* key) {
            if (!first.empty()) {
                if (!first.back()) out << ',';
                out << '\n' << std::string(2 * first.size(), ' ');
                first.back() = false;
            }
            if (key) {
                quoted(key);
                out << ": ";
            }
        }
        void finite(double value) {
            if (!std::isfinite(value))
                throw OpalException("ResultJsonWriter", "Cannot serialize a nonfinite result.");
            out << value;
        }
        template <class F>
        void container(const char* key, char open, char close, F contents) {
            prefix(key);
            out << open;
            first.push_back(true);
            contents();
            const bool empty = first.back();
            first.pop_back();
            if (!empty) out << '\n' << std::string(2 * first.size(), ' ');
            out << close;
        }

    public:
        ResultJsonWriter() {
            out.imbue(std::locale::classic());
            out << std::setprecision(std::numeric_limits<double>::max_digits10);
        }
        template <class F>
        void object(const char* key, F contents) {
            container(key, '{', '}', contents);
        }
        template <class F>
        void array(const char* key, F contents) {
            container(key, '[', ']', contents);
        }
        void text(const char* key, const std::string& value) {
            prefix(key);
            quoted(value);
        }
        void number(const char* key, double value) {
            prefix(key);
            finite(value);
        }
        void boolean(const char* key, bool value) {
            prefix(key);
            out << (value ? "true" : "false");
        }
        void null(const char* key) {
            prefix(key);
            out << "null";
        }
        template <class Values>
        void numbers(const char* key, const Values& values) {
            prefix(key);
            out << '[';
            bool comma = false;
            for (double value : values) {
                if (comma) out << ", ";
                finite(value);
                comma = true;
            }
            out << ']';
        }
        std::string str() const { return out.str() + '\n'; }
    };
    // Keep the previous result intact until a complete replacement is closed.
    // A unique sibling file also avoids collisions between independent runs.
    inline void writeResult(
            const std::string& command, const std::string& output, const std::string& json) {
        const auto require = [&](bool ok, const std::string& message) {
            if (!ok) throw OpalException(command, message);
        };
        std::string temporary = output + ".tmp.XXXXXX";
        const int fd          = ::mkstemp(temporary.data());
        if (fd < 0)
            throw OpalException(
                    command,
                    "Cannot create temporary OUTPUT for " + output + ": " + std::strerror(errno));
        struct Cleanup {
            const std::string& path;
            ~Cleanup() {
                std::error_code ignored;
                std::filesystem::remove(path, ignored);
            }
        } cleanup{temporary};
        std::unique_ptr<FILE, decltype(&std::fclose)> file(::fdopen(fd, "w"), &std::fclose);
        if (!file) {
            const std::string message = std::strerror(errno);
            ::close(fd);
            throw OpalException(
                    command, "Cannot open temporary OUTPUT for " + output + ": " + message);
        }
        require(std::fwrite(json.data(), 1, json.size(), file.get()) == json.size(),
                "Error writing OUTPUT: " + output);
        require(std::fclose(file.release()) == 0, "Error closing OUTPUT: " + output);
        std::error_code error;
        std::filesystem::rename(temporary, output, error);
        require(!error, "Cannot replace OUTPUT " + output + ": " + error.message());
    }
}  // namespace opalx::report
#endif
