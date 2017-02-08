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
                    "* The design energy of dipoles is now expected to be in MeV instead of eV.\n"
                    });

        changes.insert({109,
                    "* Beamlines are now placed in 3-dimensional space. Make sure that you\n"
                    "  use apertures to limit the range of the elements. Default aperture \n"
                    "  has circular shape with diameter 1 meter.\n"
                    "\n"
                    "* Beamlines containing a cathode have to have a SOURCE element to indicate\n"
                    "  this fact\n"
                    });
    }
}