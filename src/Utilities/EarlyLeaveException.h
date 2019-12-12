#ifndef OPAL_EARLY_LEAVE_EXCEPTION_HH
#define OPAL_EARLY_LEAVE_EXCEPTION_HH

#include "Utilities/ClassicException.h"


class EarlyLeaveException: public ClassicException {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    explicit EarlyLeaveException(const std::string &meth, const std::string &msg);

    EarlyLeaveException(const EarlyLeaveException &);
    virtual ~EarlyLeaveException();

    /// Return the message string for the exception.
    using ClassicException::what;

    /// Return the name of the method or function which detected the exception.
    using ClassicException::where;

private:

    // Not implemented.
    EarlyLeaveException();
};

#endif
