#ifndef OPAL_ErrorParser_HH
#define OPAL_ErrorParser_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorParser
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "AbstractObjects/Directory.h"


// Class ErrorParser
// ------------------------------------------------------------------------
/// The parser used in the error module.
//  As long as control remains in this class, OPAL recognises only the
//  commands allowed in error mode.  Thus this parser has its own command
//  directory with a find() method which is used to find commands.

class ErrorParser: public OpalParser {

public:

    ErrorParser();
    virtual ~ErrorParser();

protected:

    /// Find object by name in the error command directory.
    virtual Object *find(const std::string &) const;

private:

    // Not implemented.
    ErrorParser(const ErrorParser &);
    void operator=(const ErrorParser &);

    // The sub-command directory.
    Directory ErrorDirectory;
};

#endif // OPAL_ErrorParser_HH
