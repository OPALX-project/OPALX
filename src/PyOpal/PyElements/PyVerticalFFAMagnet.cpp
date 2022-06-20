#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Elements/OpalVerticalFFAMagnet.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyVerticalFFAMagnet {

std::string track_run_docstring = std::string();


const char* module_docstring = "build a vertical_ffa_magnet";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::attributes = {
    {"B0", "b0", "", PyOpalObjectNS::DOUBLE},
    {"FIELD_INDEX", "field_index", "", PyOpalObjectNS::DOUBLE},
    {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
    {"MAX_HORIZONTAL_POWER", "max_horizontal_power", "", PyOpalObjectNS::INT},
    {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
    {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
    {"BB_LENGTH", "bb_length", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT_NEG_EXTENT", "height_neg_extent", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT_POS_EXTENT", "height_pos_extent", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::classDocstring = "";

BOOST_PYTHON_MODULE(vertical_ffa_magnet) {
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet> element;
    auto elementClass = element.make_class("VerticalFFAMagnet");
    element.addGetOpalElement(elementClass);
    elementClass.def("get_field_value", &PyOpalObjectNS::getFieldValue<OpalVerticalFFAMagnet>);
}

}
}
