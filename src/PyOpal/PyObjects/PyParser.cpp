#include <string>
#include <boost/python.hpp>

#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "Main.cpp"

namespace PyOpal {
namespace PyParser {

std::string initialise_from_opal_file_docstring = 
std::string("Initialise from opal file\n")+
std::string("- file_name: string corresponding to the file name of the OPAL\n")+
std::string("  file.\n")+
std::string("Returns None.\n");

void initialise_from_opal_file(std::string file_name) {
    std::string exe("parser");
    char* argvr[3];
    // argv must be NULL terminated array (a week of my life figuring that one)
    argvr[0] = new char[exe.length()+1]();
    strcpy(argvr[0], exe.c_str());
    argvr[1] = new char[file_name.length()+1]();
    strcpy(argvr[0], file_name.c_str());
    argvr[2] = nullptr;
    opalMain(2, argvr);
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