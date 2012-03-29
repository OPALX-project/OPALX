#include "SwitcherError.h"

SwitcherError::SwitcherError(const string &meth, const string &msg):
    ClassicException(meth, msg)
{}

SwitcherError::SwitcherError(const SwitcherError &rhs):
    ClassicException(rhs)
{}

SwitcherError::~SwitcherError()
{}
