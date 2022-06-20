#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyBeam.h"

#include "Structure/Beam.h"

namespace PyOpal {
namespace PyBeamNS {

std::string track_run_docstring = std::string();


const char* module_docstring = "build a tracking object";

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Beam>::attributes = {
    {"PARTICLE", "particle", "", PyOpalObjectNS::STRING},
    {"MASS", "mass", "", PyOpalObjectNS::DOUBLE},
    {"CHARGE", "charge", "", PyOpalObjectNS::DOUBLE},
    {"ENERGY", "energy", "", PyOpalObjectNS::DOUBLE},
    {"PC", "momentum", "", PyOpalObjectNS::DOUBLE},
    {"GAMMA", "gamma", "", PyOpalObjectNS::DOUBLE},
    {"BCURRENT", "beam_current", "", PyOpalObjectNS::DOUBLE},
    {"EX", "ex", "", PyOpalObjectNS::DOUBLE},
    {"EY", "ey", "", PyOpalObjectNS::DOUBLE},
    {"ET", "et", "", PyOpalObjectNS::DOUBLE},
    {"BFREQ", "beam_frequency", "", PyOpalObjectNS::DOUBLE},
    {"NPART", "number_of_particles", "", PyOpalObjectNS::DOUBLE},
    {"NSLICE", "number_of_slices", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<Beam>::classDocstring = "";

BOOST_PYTHON_MODULE(beam) {
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<Beam> aBeam;
    auto beamClass = aBeam.make_class("Beam");
    aBeam.addRegister(beamClass);
}

} // PyBeamNS
} // PyOpal

