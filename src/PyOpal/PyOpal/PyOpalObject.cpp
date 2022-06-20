#include "PyOpal/PyOpalObject.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyOpalObjectNS {

std::map<AttributeType, std::string> attributeName = std::map<AttributeType, std::string>({
    {DOUBLE, "float"},
    {STRING, "string"},
    {BOOL, "bool"},
    {INT, "int"},
    {FLOATLIST, "list of floats"}
});

} // PyElementNS
} // PyOpal
