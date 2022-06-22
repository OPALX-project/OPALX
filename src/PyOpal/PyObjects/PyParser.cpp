#include <string>
#include <boost/python.hpp>

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
    /*
    argvr[0] = new char[exe.length()+2]();
    strcpy(argvr[0], exe.c_str());

    argvr[1] = new char[file_name.length()+2]();
    strcpy(argvr[0], file_name.c_str());
    */
    argvr[2] = nullptr;
    int error_code = opalMain(2, argvr);
    return error_code;
}

std::string module_docstring =
"The parser module is used to load an OPAL input file from within python";
 
BOOST_PYTHON_MODULE(parser) { // parser is a python internal library
    PyOpal::ExceptionTranslation::registerExceptions();
    boost::python::scope().attr("__doc__") = module_docstring.c_str();
    boost::python::def("initialise_from_opal_file",
            initialise_from_opal_file,
            boost::python::args("file_name"),
            initialise_from_opal_file_docstring.c_str()
    );
}

} // namespace PyParser
} // namespace PyOpal