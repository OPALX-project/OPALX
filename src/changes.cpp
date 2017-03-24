#include <map>
#include <string>

namespace Versions {
    std::map<unsigned int, std::string> changes;

    void fillChanges() {
        if (changes.size() > 0) return;

        changes.insert({105,
                    "* The normalization of the 2-dimensional field maps has changed.\n"
                    "  Instead of normalizing with the overall maximum value of longitudinal\n"
                    "  component Opal now uses the maximum value on axis.\n"
                    "\n"
                    " The parser has been modified to check the type of all variables. All real \n"
                    " variables have to be prefixed with the keyword REAL.\n"
                    });
    }
}