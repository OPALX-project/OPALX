#include <Python.h>
#include <structmember.h>

#include <exception>
#include <iostream>
#include <boost/python.hpp>

#include "Classic/Fields/Interpolation/NDGrid.h"
#include "Classic/Fields/Interpolation/PPSolveFactory.h"
#include "Classic/Fields/Interpolation/PolynomialPatch.h"

#include "PyOpal/ExceptionTranslation.h"

namespace PyOpal {

namespace PyPolynomialPatch {

using namespace interpolation;
namespace py = boost::python;

PolynomialPatch* initialiseFromSolveFactory(
                  NDGrid* points,
                  boost::python::list values,
                  int polyPatchOrder,
                  int smoothingOrder) {
    int gLength = py::len(values);
    std::vector<std::vector<double> > valuesVec(gLength);
    for (int i = 0; i < gLength; ++i) {
        int lineLength = py::len(values[i]);
        valuesVec[i] = std::vector<double>(lineLength);
        for (int j = 0; j < lineLength; ++j) {
            valuesVec[i][j] = py::extract<double>(values[i][j]);
        }
    }
    // points is owned by Python; clone to get a clean copy...
    Mesh* pointsClone = points->clone();
    PolynomialPatch* patch = PPSolveFactory(pointsClone,
                                            valuesVec,
                                            polyPatchOrder,
                                            polyPatchOrder).solve();
    return patch;
}

py::list function(PolynomialPatch* patch, py::list point) {
    int pointDim = patch->getPointDimension();
    int valueDim = patch->getValueDimension();
    if (py::len(point) != pointDim) {
        //error
    }
    std::vector<double> pointVec(pointDim);
    for (int i = 0; i < pointDim; ++i) {
        pointVec[i] = py::extract<double>(point[i]);
    }
    std::vector<double> valueVec(valueDim);
    patch->function(&pointVec[0], &valueVec[0]);
    py::list value = py::list();
    for (int i = 0; i < valueDim; ++i) {
        value.append(valueVec[i]);
    }
    return value;
}


const char* module_docstring = "polynomial_patch module returns the field";

BOOST_PYTHON_MODULE(polynomial_patch) {
    ExceptionTranslation::registerExceptions();
    py::class_<PolynomialPatch, boost::noncopyable>("PolynomialPatch")
        .def("initialise_from_solve_factory", &initialiseFromSolveFactory,
         py::return_value_policy<py::manage_new_object>())
        .staticmethod("initialise_from_solve_factory")
        .def("function", &function)
    ;
}

} // namespace PyPolynomialPatch
} // namespace PyOpal
