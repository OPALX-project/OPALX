#ifndef GENERALOPALEXCEPTION_H
#define GENERALOPALEXCEPTION_H

#include "Utilities/OpalException.h"

/// Generic OPAL exception for code paths that do not yet expose a more
/// domain-specific exception type.
///
/// Prefer dedicated subclasses such as parse, format, or fieldmap errors when
/// adding new code. This type preserves existing generic error semantics while
/// keeping the hierarchy rooted in OpalException.
class GeneralOpalException : public OpalException {
public:
    GeneralOpalException(const std::string& meth, const std::string& msg);

    GeneralOpalException(const GeneralOpalException&);
    ~GeneralOpalException() override;

private:
    // Not implemented.
    GeneralOpalException();
};

#endif
