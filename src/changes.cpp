#include <map>
#include <string>

#include "changes.h"

#include "Utility/Inform.h"
Inform* gmsg;

namespace Versions {
    std::map<unsigned int, std::string> changes;

    // No OPALX release has changed the meaning of an existing input yet.
    void fillChanges() {}
}  // namespace Versions