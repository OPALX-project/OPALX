#include "Utilities/EarlyLeaveException.h"


EarlyLeaveException::EarlyLeaveException(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}


EarlyLeaveException::EarlyLeaveException(const EarlyLeaveException &rhs):
    ClassicException(rhs)
{}


EarlyLeaveException::~EarlyLeaveException()
{}
