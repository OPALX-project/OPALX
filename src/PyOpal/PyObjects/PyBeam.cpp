#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Structure/Beam.h"

namespace PyOpal {
namespace PyBeamNS {

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Beam>::attributes = {
    {"PARTICLE", "particle", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"MASS", "mass", "", PyOpalObjectNS::DOUBLE},
    {"CHARGE", "charge", "", PyOpalObjectNS::DOUBLE},
    {"ENERGY", "energy", "", PyOpalObjectNS::DOUBLE},
    {"PC", "momentum", "", PyOpalObjectNS::DOUBLE},
    {"GAMMA", "gamma", "", PyOpalObjectNS::DOUBLE},
    {"BCURRENT", "beam_current", "", PyOpalObjectNS::DOUBLE},
    {"BFREQ", "beam_frequency", "", PyOpalObjectNS::DOUBLE},
    {"NPART", "number_of_particles", "", PyOpalObjectNS::DOUBLE},
};

BOOST_PYTHON_MODULE(beam) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    PyOpalObjectNS::PyOpalObject<Beam> aBeam;
    auto beamClass = aBeam.make_class("Beam");
    aBeam.addRegister(beamClass);
}

} // PyBeamNS
} // PyOpal

