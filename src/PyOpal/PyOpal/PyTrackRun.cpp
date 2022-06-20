#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyTrackRun.h"
#include "PyOpal/Globals.h"

#include "Track/TrackRun.h"


extern Inform *gmsg;

namespace PyOpal {
namespace PyTrackRunNS {

std::string track_run_docstring = std::string();

const char* module_docstring = "build a tracking object";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TrackRun>::attributes = {
    {"METHOD", "method", "", PyOpalObjectNS::STRING},
    {"TURNS", "turns", "", PyOpalObjectNS::DOUBLE},
    {"MBMODE", "multibunch_mode", "", PyOpalObjectNS::STRING},
    {"PARAMB", "multibunch_control", "", PyOpalObjectNS::DOUBLE},
    {"MB_ETA", "multibunch_scale", "", PyOpalObjectNS::DOUBLE},
    {"MB_BINNING", "multibunch_binning", "", PyOpalObjectNS::DOUBLE},
    {"BEAM", "beam_name", "", PyOpalObjectNS::STRING},
    {"FIELDSOLVER", "field_solver", "", PyOpalObjectNS::STRING},
    {"BOUNDARYGEOMETRY", "boundary_geometry", "", PyOpalObjectNS::STRING},
    {"DISTRIBUTION", "distribution", "", PyOpalObjectNS::STRING},
    {"MULTIPACTING", "multipacting", "", PyOpalObjectNS::BOOL},
    {"KEEPALIVE", "keep_alive", "", PyOpalObjectNS::BOOL},
    {"OBJECTIVES", "objectives", "", PyOpalObjectNS::STRING},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<TrackRun>::classDocstring = "";

BOOST_PYTHON_MODULE(track_run) {
    std::cerr << "TRACK RUN MODULE 1 gmsg " << gmsg << std::endl;
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    std::cerr << "TRACK RUN MODULE 2 gmsg " << gmsg << std::endl;
    PyOpalObjectNS::PyOpalObject<TrackRun> trackRun;
    auto trackClass = trackRun.make_class("TrackRun");
    trackRun.addExecute(trackClass);
    std::cerr << "TRACK RUN MODULE DONE" << std::endl;
}

} // PyTrackRun
} // PyOpal

