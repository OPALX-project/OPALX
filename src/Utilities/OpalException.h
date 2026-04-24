#ifndef OPAL_OpalException_HH
#define OPAL_OpalException_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalException.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalException
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <string>


// Class OpalException
// ------------------------------------------------------------------------
/// The base class for all OPAL exceptions.
///
/// Stores the detecting method and diagnostic message used by the top-level
/// OPAL error handler. Concrete subclasses should keep domain-specific type
/// information and derive from this class directly or indirectly.
class OpalException {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    explicit OpalException(const std::string &meth, const std::string &msg);

    OpalException(const OpalException &);
    virtual ~OpalException();

    /// Return the message string for the exception.
    virtual const std::string &what() const;

    /// Return the name of the method or function which detected the exception.
    virtual const std::string &where() const;

private:

    // Not implemented.
    OpalException();

    // The method detecting the exception and the message.
    const std::string message;
    const std::string method;
};

#endif // OPAL_OpalException_HH
