#include "PyOpal/PyCore/PyOpalObject.h"

namespace PyOpal {
namespace PyOpalObjectNS {

std::map<AttributeType, std::string> attributeName = std::map<AttributeType, std::string>({
    {DOUBLE, "float"},
    {STRING, "string"},
    {PREDEFINED_STRING, "predefined string"},
    {UPPER_CASE_STRING, "upper case string"},
    {BOOL, "bool"},
    {INT, "int"},
    {FLOAT_LIST, "list of floats"}
});
} // PyOpalObjectNS
} // PyOpal
