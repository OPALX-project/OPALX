#include <boost/python.hpp>

#include "Utilities/OpalException.h"
#include "AbsBeamline/Ring.h"
#include "Track/TrackRun.h"
#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelCyclotronTracker.h"

#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyObjects/PyField.h"


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
  "The function returns a tuple containing 6 values:\n"
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

py::object get_field_value_parallelt(double x,
                                     double y,
                                     double z,
                                     double t,
                                     ParallelTTracker* tracker) {
    throw OpalException("PyField::get_field_value_parallelt",
                        "Not implemented");
    if (tracker == NULL) {
        throw(OpalException("PyField::get_field_value_parallelt",
                            "ParallelTTracker was NULL"));
    }
    int outOfBounds = 0;
    Vector_t R(x, y, z);
    Vector_t P(0, 0, 0);
    Vector_t E, B;
    // outOfBounds = tracker->getFieldValue(R, t, B, E);
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0], B[1], B[2],
                                          E[0], E[1], E[2]);
    return value;
}

py::object get_field_value_cyclotron(double x,
                                     double y,
                                     double z,
                                     double t,
                                     ParallelCyclotronTracker* tracker) {
    throw OpalException("PyField::get_field_value_parallelcyclotron",
                        "Not implemented");
    if (tracker == NULL) {
        throw(OpalException("PyField::get_field_value_cyclotron",
                            "ParallelCyclotronTracker was NULL"));
    }
    boost::python::tuple value = boost::python::make_tuple(1,
                                          -1, -1, -1, -1, -1, -1);
    return value;

}

py::object get_field_value_ring(double x,
                                     double y,
                                     double z,
                                     double t,
                                     Ring* ring) {
    if (ring == NULL) {
        throw(OpalException("PyField::get_field_value_ring",
                            "Ring was NULL"));
    }
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe a "
           "RingDefinition was not defined or KeepAlive was False";
        throw(OpalException("PyField::get_field_value", err));
    }
    Vector_t R(x, y, z);
    Vector_t P(0, 0, 0);
    Vector_t E, B;
    int outOfBounds = ring->apply(R, P, t, E, B);
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0]/10., B[1]/10., B[2]/10.,
                                          E[0], E[1], E[2]);
    return value;
}

py::object get_field_value(double x, double y, double z, double t) {
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring != NULL) {
        return get_field_value_ring(x, y, z, t, ring);
    }
    /*
    Tracker* tracker = TrackRun::getTracker();
    ParallelTTracker* trackerT = dynamic_cast<ParallelTTracker*>(tracker);
    if (trackerT != NULL) {
        return get_field_value_parallelt(x, y, z, t, trackerT);       
    }
    ParallelCyclotronTracker* trackerCycl = dynamic_cast<ParallelCyclotronTracker*>(tracker);
    if (trackerCycl != NULL) {
        return get_field_value_cyclotron(x, y, z, t, trackerCycl);
    }
    */
    throw(OpalException("PyField::get_field_value",
                        "Could not find a Ring, ParallelTTracker or ParallelCyclotronTracker"));
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

