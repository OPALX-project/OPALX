#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyFieldSolver.h"

#include "Structure/FieldSolver.h"

namespace PyOpal {
namespace PyFieldSolverNS {

std::string track_run_docstring = std::string();


const char* module_docstring = "build a tracking object";

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<FieldSolver>::attributes = {
    {"FSTYPE", "field_solver_type", "", PyOpalObjectNS::STRING},
    {"MX", "mesh_size_x", "", PyOpalObjectNS::DOUBLE},
    {"MY", "mesh_size_y", "", PyOpalObjectNS::DOUBLE},
    {"MT", "mesh_size_t", "", PyOpalObjectNS::DOUBLE},
    {"PARFFTX", "parallelise_fft_x", "", PyOpalObjectNS::BOOL},
    {"PARFFTY", "parallelise_fft_y", "", PyOpalObjectNS::BOOL},
    {"PARFFTT", "parallelise_fft_t", "", PyOpalObjectNS::BOOL},
    {"BCFFTX", "fft_boundary_x", "", PyOpalObjectNS::STRING},
    {"BCFFTY", "fft_boundary_y", "", PyOpalObjectNS::STRING},
    {"BCFFTZ", "fft_boundary_z", "", PyOpalObjectNS::STRING},
    {"GREENSF", "greens_function", "", PyOpalObjectNS::STRING},
    {"BBOXINCR", "bounding_box_increase", "", PyOpalObjectNS::DOUBLE},
    {"GEOMETRY", "geometry", "", PyOpalObjectNS::STRING},
    {"ITSOLVER", "iterative_solver", "", PyOpalObjectNS::STRING},
    {"INTERPL", "interpolation", "", PyOpalObjectNS::STRING},
    {"TOL", "tolerance", "", PyOpalObjectNS::DOUBLE},
    {"MAXITERS", "max_iterations", "", PyOpalObjectNS::DOUBLE},
    {"PRECMODE", "preconditioner_mode", "", PyOpalObjectNS::STRING},
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
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<FieldSolver> fs;
    auto fsClass = fs.make_class("FieldSolver");
    fs.addRegister(fsClass);
}

} // PyFieldSolverNS
} // PyOpal

