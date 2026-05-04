#include "Utilities/GeneralOpalException.h"

GeneralOpalException::GeneralOpalException(const std::string& meth, const std::string& msg)
    : OpalException(meth, msg) {}

GeneralOpalException::GeneralOpalException(const GeneralOpalException& rhs) : OpalException(rhs) {}

GeneralOpalException::~GeneralOpalException() {}
