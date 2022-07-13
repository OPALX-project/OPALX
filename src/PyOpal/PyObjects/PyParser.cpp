#include <sstream>
#include <string>
#include <boost/python.hpp>

#include "AbstractObjects/OpalData.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "Main.cpp"

namespace PyOpal {
namespace PyParser {

std::string initialise_from_opal_file_docstring = 
"Initialise from opal file\n"
"- file_name: string corresponding to the file name of the OPAL\n"
"  file.\n"
"Note that if file_name is not valid, OPAL may terminate the python script\n"
"execution abnormally (without the usual python exit semantics).\n"
"\n"
"Returns an integer; 0 for successful execution or non-zero if an error\n"
"occurred.\n";

int initialise_from_opal_file(std::string file_name) {
    std::string exe("parser");
    char* argvr[3];
    // argv must be NULL terminated array
    argvr[0] = exe.data();
    argvr[1] = file_name.data();
    argvr[2] = nullptr;
    int error_code = opalMain(2, argvr);
    return error_code;
}

std::string list_objects_docstring = 
"List the objects that are known by the Opal parser, either through calls to\n"
"initialise_from_opal_file or through calls directly to the python API\n"
"\n"
"Returns a list of strings, each one corresponding to the name of a\n"
"particular object\n";

boost::python::list list_objects() {
    std::vector<std::string> names = OpalData::getInstance()->getAllNames();
    boost::python::list pynames;
    for (size_t i = 0; i < names.size(); ++i) {
        pynames.append(names[i]);
    } 
    return pynames;
}


std::string module_docstring =
"The parser module is used to load an OPAL input file from within python";
 
BOOST_PYTHON_MODULE(parser) {
    PyOpal::ExceptionTranslation::registerExceptions();
    boost::python::scope().attr("__doc__") = module_docstring.c_str();
    boost::python::def("initialise_from_opal_file",
            initialise_from_opal_file,
            boost::python::args("file_name"),
            initialise_from_opal_file_docstring.c_str()
    );
    boost::python::def("list_objects",
            list_objects,
            list_objects_docstring.c_str()
    );
}

} // namespace PyParser
} // namespace PyOpal