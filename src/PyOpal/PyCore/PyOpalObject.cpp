#include "PyOpal/PyCore/PyOpalObject.h"

namespace PyOpal {
namespace PyOpalObjectNS {

std::map<AttributeType, std::string> attributeName = std::map<AttributeType, std::string>({
    {DOUBLE, "float"},
    {STRING, "string"},
    {BOOL, "bool"},
    {INT, "int"},
    {FLOATLIST, "list of floats"}
});

} // PyOpalObjectNS
} // PyOpal
