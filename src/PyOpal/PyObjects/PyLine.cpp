#include "PyOpal/ExceptionTranslation.h"
#include "Lines/Line.h"
#include "PyOpal/PyLine.h"
#include "PyOpal/Globals.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyLineNS {

const char* module_docstring =
"The line module handles building of a line of elements in OPAL";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TBeamline<Element> >::attributes = {
    {"L", "length", "", PyOpalObjectNS::DOUBLE},
    {"ORIGIN", "origin", "", PyOpalObjectNS::STRING},
    {"ORIENTATION", "orientation", "", PyOpalObjectNS::STRING},
    {"X", "x", "", PyOpalObjectNS::DOUBLE},
    {"Y", "y", "", PyOpalObjectNS::DOUBLE},
    {"Z", "z", "", PyOpalObjectNS::DOUBLE},
    {"THETA", "theta", "", PyOpalObjectNS::DOUBLE},
    {"PHI", "phi", "", PyOpalObjectNS::DOUBLE},
    {"PSI", "psi", "", PyOpalObjectNS::DOUBLE}
};

template <>
std::string PyOpalObjectNS::PyOpalObject<TBeamline<Element> >::classDocstring = "";

BOOST_PYTHON_MODULE(line) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    PyLine element;
    auto lineClass = element.make_class("Line");
    // https://docs.python.org/3/library/collections.abc.html
    // I tried to pull everything from:
    // mutable sequence, sequence, reversible, iterator, iterable, container
    // so should look like a python list
    lineClass
        .def("__len__", &PyLine::getLength)
        .def("__getitem__", &PyLine::getElement)
        .def("__setitem__", &PyLine::setElement)
        .def("append", &PyLine::append)
        /*
        .def("__iter__", &PyElement<Line>::dummyGet<double>)
        .def("__next__", &PyElement<Line>::dummyGet<double>)
        .def("__delitem__", &PyElement<Line>::dummyGet<double>)
        .def("__contains__", &PyElement<Line>::dummySet<double>)
        .def("__reversed__", &PyElement<Line>::dummySet<double>)
        .def("index", &PyElement<Line>::dummySet<double>)
        .def("count", &PyElement<Line>::dummySet<double>)
        .def("reverse", &PyElement<Line>::dummySet<double>)
        .def("extend", &PyElement<Line>::dummySet<double>)
        .def("pop", &PyElement<Line>::dummySet<double>)
        .def("remove", &PyElement<Line>::dummySet<double>)
        .def("append", &PyElement<Line>::dummySet<double>)
        .def("__iadd__", &PyElement<Line>::dummySet<double>)*/
        ;
    element.addGetOpalElement(lineClass);
    lineClass.def("register", &PyLine::registerObject);

}

}
}
