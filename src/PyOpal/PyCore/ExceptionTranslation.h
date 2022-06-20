#ifndef PYOPAL_PYCORE_EXCEPTIONTRANSLATION_H
#define PYOPAL_PYCORE_EXCEPTIONTRANSLATION_H 1

#include <Python.h>
#include <structmember.h>

#include <string>
#include <exception>

#include <boost/python.hpp>

namespace PyOpal {
namespace py = boost::python;

/** Exception translation uses boost::python hooks to wrap C++ exceptions */
namespace ExceptionTranslation {

/** Register exception translations with boost */
void registerExceptions();

/** Translates the std::exceptions into a RuntimeError */
template <class T>
void translateException(T const& exception) {
    PyErr_SetString(PyExc_RuntimeError, exception.what());
}

/** Translates the OpalExceptions into a RuntimeError */
template <class T>
void translateOpalException(T const& exception) {
    std::string msg = exception.what()+" in C++ method "+exception.where();
    PyErr_SetString(PyExc_RuntimeError, msg.c_str());
} 


}
}

#endif