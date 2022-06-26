#include "Distribution/Distribution.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"

namespace PyOpal {
namespace PyDistributionNS {

std::string distribution_docstring = std::string();

const char* module_docstring = "build a distribution object";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Distribution>::attributes = {
    {"TYPE", "type", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"FNAME", "fname", "", PyOpalObjectNS::STRING},
    {"INPUTMOUNITS", "momentum_units", "", PyOpalObjectNS::PREDEFINED_STRING},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<Distribution>::classDocstring = "";

void registerDistribution(PyOpalObjectNS::PyOpalObject<Distribution>& dist) {
    Object* obj = &(*dist.getOpalShared());
    obj->update();
    OpalData::getInstance()->define(obj);
}

BOOST_PYTHON_MODULE(distribution) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<Distribution> distributionObject;
    auto distributionClass = distributionObject.make_class("Distribution");
    distributionObject.addExecute(distributionClass);
    distributionClass.def("register", &registerDistribution);

}

} // PyTrackRun
} // PyOpal

