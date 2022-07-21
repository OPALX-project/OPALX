#include <boost/python.hpp>

#include "Utilities/OpalException.h"
#include "AbsBeamline/Ring.h"
#include "Track/TrackRun.h"
#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelCyclotronTracker.h"

#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"

namespace PyOpal {
namespace Field {

std::string field_docstring = 
  "field module enables user to get the field at a point";

std::string get_field_value_docstring =
  "Get the field value at a point in the field map.\n"
  "\n"
  "The field lookup is performed against the last RINGDEFINITION that was\n"
  "instantiated. This should be instantiated by calling\n"
  "pyopal.parser.initialise_from_opal_file\n"
  "\n"
  "Parameters\n"
  "----------\n"
  "x : float\n"
  "    x position [m]\n"
  "y : float\n"
  "    y position [m]\n"
  "z : float\n"
  "    z position [m]\n"
  "t: float\n"
  "    time [ns]\n"
  "\n"
  "Returns\n"
  "-------\n"
  "The function returns a tuple containing 7 values:\n"
  "out of bounds : int\n"
  "    1 if the event was out of the field map boundary, else 0.\n"
  "Bx : float\n"
  "    x magnetic field [T]\n"
  "By : float\n"
  "    y magnetic field [T]\n"
  "Bz : float\n"
  "    z magnetic field [T]\n"
  "Ex : float\n"
  "    x electric field\n"
  "Ey : float\n"
  "    y electric field\n"
  "Ez : float\n"
  "    z electric field\n";

py::object get_field_value_cyclotron(double x,
                                     double y,
                                     double z,
                                     double t,
                                     ParallelCyclotronTracker* tracker) {
    if (tracker == NULL) {
        throw(OpalException("PyField::get_field_value_cyclotron",
                            "ParallelCyclotronTracker was NULL"));
    }
    Vector_t R(x, y, z);
    Vector_t P, B, E;
    int outOfBounds = tracker->computeExternalFields_m(R, P, t, E, B);
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0], B[1], B[2], E[0], E[1], E[2]);
    return value;

}

py::object get_field_value(double x, double y, double z, double t) {
    /*
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring != NULL) {
        return get_field_value_ring(x, y, z, t, ring);
    }*/
    std::shared_ptr<Tracker> tracker = TrackRun::getTracker();
    ParallelCyclotronTracker* trackerCycl = 
                        dynamic_cast<ParallelCyclotronTracker*>(tracker.get());
    if (trackerCycl != nullptr) {
        return get_field_value_cyclotron(x, y, z, t, trackerCycl);
    }
    throw(OpalException("PyField::get_field_value",
                        "Could not find a ParallelCyclotronTracker"));
}

BOOST_PYTHON_MODULE(field) {
    ExceptionTranslation::registerExceptions();
    py::scope().attr("__doc__") = field_docstring.c_str();
    py::def("get_field_value",
            get_field_value,
            py::args("x", "y", "z", "t"),
            get_field_value_docstring.c_str()
    );
}

}
}

