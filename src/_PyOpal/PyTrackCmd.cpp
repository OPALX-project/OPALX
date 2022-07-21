#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyTrackRun.h"

#include "Track/TrackCmd.h"

namespace PyOpal {
namespace PyTrackRunNS {

const char* module_docstring = "build a tracking object";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TrackCmd>::attributes = {
    {"LINE", "line", "", PyOpalObjectNS::STRING},
    {"BEAM", "beam", "", PyOpalObjectNS::STRING},
    {"DT", "time_steps", "", PyOpalObjectNS::DOUBLE}, // array
    {"DTSCINIT", "dt_space_charge", "", PyOpalObjectNS::DOUBLE},
    {"DTAU", "dtau", "", PyOpalObjectNS::DOUBLE},
    {"T0", "t0", "", PyOpalObjectNS::DOUBLE},
    {"MAXSTEPS", "max_steps", "", PyOpalObjectNS::DOUBLE},
    {"STEPSPERTURN", "steps_per_turn", "", PyOpalObjectNS::DOUBLE},
    {"ZSTART", "z_start", "", PyOpalObjectNS::DOUBLE},
    {"ZSTOP", "z_stop", "", PyOpalObjectNS::DOUBLE},
    {"TIMEINTEGRATOR", "time_integrator", "", PyOpalObjectNS::STRING},
    {"MAP_ORDER", "map_order", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<TrackCmd>::classDocstring = "";

void executeWrapper(PyOpalObjectNS::PyOpalObject<TrackCmd>& cmd) {
    std::shared_ptr<TrackCmd> objectPtr = cmd.getOpalShared();
    objectPtr->setIsParseable(false);
    objectPtr->execute();
}

BOOST_PYTHON_MODULE(track) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<TrackCmd> trackCmd;
    auto trackClass = trackCmd.make_class("Track");
    trackClass.def("execute", &executeWrapper);
}

} // PyTrackRun
} // PyOpal

