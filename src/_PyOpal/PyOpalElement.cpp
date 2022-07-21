#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyOpalObject.h"

#include "Elements/OpalElement.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyOpalElement {

std::string track_run_docstring = std::string();


const char* module_docstring = "opal element base class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalElement>::attributes = {
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalElement>::classDocstring = "";


template <>
PyOpalObjectNS::PyOpalObject<OpalElement>::PyOpalObject() : object_m(NULL) {}

BOOST_PYTHON_MODULE(opal_element) {
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalElement> element;
    element.make_class("OpalElement");
}

}
}