#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Track/TrackCmd.h"

namespace PyOpal {
namespace PyTrackCmdNS {

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TrackCmd>::attributes = {
    {"LINE", "line", "", PyOpalObjectNS::STRING},
    {"BEAM", "beam", "", PyOpalObjectNS::STRING},
    {"DT", "time_steps", "", PyOpalObjectNS::DOUBLE}, // array
    {"DTSCINIT", "dt_space_charge", "", PyOpalObjectNS::DOUBLE},
    {"DTAU", "dtau", "", PyOpalObjectNS::DOUBLE},
    {"T0", "t0", "", PyOpalObjectNS::DOUBLE},
    {"MAXSTEPS", "max_steps", "", PyOpalObjectNS::FLOAT_LIST},
    {"STEPSPERTURN", "steps_per_turn", "", PyOpalObjectNS::DOUBLE},
    {"ZSTART", "z_start", "", PyOpalObjectNS::DOUBLE},
    {"ZSTOP", "z_stop", "", PyOpalObjectNS::FLOAT_LIST},
    {"TIMEINTEGRATOR", "time_integrator", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"MAP_ORDER", "map_order", "", PyOpalObjectNS::DOUBLE},
};

// Can't use the default PyObject execute function because we need to call
// setIsParseable to false (otherwise OPAL will try to parse it as an OPAL file)
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

} // PyTrackCmd
} // PyOpal

