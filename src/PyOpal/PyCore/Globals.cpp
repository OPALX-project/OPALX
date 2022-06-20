#include <gsl/gsl_errno.h>
#include "Utilities/OpalException.h"

// Note the gymnastics here - we only want to define gmsg and ippl once
#define PYOPAL_GLOBALS_C
#include "PyOpal/PyCore/Globals.h"

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int line,
                         int gsl_errno) {
        throw OpalException(file, reason);
        if (line || gsl_errno) {;} // disable gcc warning; does nothing
    }
}

namespace PyOpal {
namespace Globals {
void Initialise() {
    if (gmsg == nullptr) {
        gmsg = new Inform("OPAL");
    }
    if (ippl == nullptr) {
        int argc = 3;

        char* argvr[] = {
            (char*)("pyopal"), 
            (char*)("--processes"),
            (char*)("3"),
            nullptr
        };
        char** argv = argvr;
        // Ippl is a typedef of IpplInfo in ippl/Utilities
        ippl = new Ippl(argc, argv);
    }
    gsl_set_error_handler(&errorHandlerGSL);
}
}
}
