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
                    "* OPAL-T: The design energy of dipoles is now expected to be in MeV instead\n"
                    "  of eV.\n"
                    });

        changes.insert({109,
                    "* OPAL-T: Beamlines are now placed in 3-dimensional space. Make sure that\n"
                    "  you use apertures to limit the range of the elements. Default aperture \n"
                    "  has circular shape with diameter 1 meter.\n"
                    "\n"
                    "* OPAL-T: Beamlines containing a cathode have to have a SOURCE element to\n"
                    "  indicate this fact\n"
                    });
    }
}