#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyOpalObject.h"

#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyOpalLocalCartesianOffset {

using OpalOffset::OpalLocalCartesianOffset;

std::string track_run_docstring = std::string();


const char* module_docstring = "build a local cartesian offset";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalLocalCartesianOffset>::attributes = {
    {"END_POSITION_X", "end_position_x", "", PyOpalObjectNS::DOUBLE},
    {"END_POSITION_Y", "end_position_y", "", PyOpalObjectNS::DOUBLE},
    {"END_POSITION_Z", "end_position_z", "", PyOpalObjectNS::DOUBLE},
    {"END_NORMAL_X", "end_normal_x", "", PyOpalObjectNS::DOUBLE},
    {"END_NORMAL_Y", "end_normal_y", "", PyOpalObjectNS::DOUBLE},
    {"END_ROTATION_X", "end_rotation_x", "", PyOpalObjectNS::DOUBLE},
    {"END_ROTATION_Y", "end_rotation_y", "", PyOpalObjectNS::DOUBLE},
    {"END_ROTATION_Z", "end_rotation_z", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalLocalCartesianOffset>::classDocstring = "";

BOOST_PYTHON_MODULE(local_cartesian_offset) {
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalLocalCartesianOffset> element;
    auto elementClass = element.make_class("LocalCartesianOffset");
    element.addGetOpalElement(elementClass);
}

}
}
