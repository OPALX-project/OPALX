
#include "Structure/FieldSolver.h"

#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

namespace PyOpal {
namespace PyFieldSolverNS {

std::string track_run_docstring = std::string();


const char* module_docstring = "build a tracking object";

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<FieldSolver>::attributes = {
    {"FSTYPE", "field_solver_type", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"MX", "mesh_size_x", "", PyOpalObjectNS::DOUBLE},
    {"MY", "mesh_size_y", "", PyOpalObjectNS::DOUBLE},
    {"MT", "mesh_size_t", "", PyOpalObjectNS::DOUBLE},
    {"PARFFTX", "parallelise_fft_x", "", PyOpalObjectNS::BOOL},
    {"PARFFTY", "parallelise_fft_y", "", PyOpalObjectNS::BOOL},
    {"PARFFTT", "parallelise_fft_t", "", PyOpalObjectNS::BOOL},
    {"BCFFTX", "fft_boundary_x", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BCFFTY", "fft_boundary_y", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BCFFTZ", "fft_boundary_z", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"GREENSF", "greens_function", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BBOXINCR", "bounding_box_increase", "", PyOpalObjectNS::DOUBLE},
    {"GEOMETRY", "geometry", "", PyOpalObjectNS::UPPER_CASE_STRING},
    {"ITSOLVER", "iterative_solver", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"INTERPL", "interpolation", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"TOL", "tolerance", "", PyOpalObjectNS::DOUBLE},
    {"MAXITERS", "max_iterations", "", PyOpalObjectNS::DOUBLE},
    {"PRECMODE", "preconditioner_mode", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"RC", "cutoff_radius", "", PyOpalObjectNS::DOUBLE},
    {"ALPHA", "alpha", "", PyOpalObjectNS::DOUBLE},
    {"EPSILON", "epsilon", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<FieldSolver>::classDocstring = "";

void registerFieldSolver(PyOpalObjectNS::PyOpalObject<FieldSolver>& fs) {
    Object* obj = &(*fs.getOpalShared());
    OpalData::getInstance()->define(obj);
}

BOOST_PYTHON_MODULE(field_solver) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<FieldSolver> fs;
    auto fsClass = fs.make_class("FieldSolver");
    fs.addRegister(fsClass);
}

} // PyFieldSolverNS
} // PyOpal

