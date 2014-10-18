#ifndef OPAL_MatchParser_HH
#define OPAL_MatchParser_HH

// ------------------------------------------------------------------------
// $RCSfile: MatchParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchParser
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "AbstractObjects/Directory.h"


// Class MatchParser
// ------------------------------------------------------------------------
/// The parser used in the OPAL match module.
//  As long as control remains in this class, OPAL recognises only the
//  commands allowed in matching mode.  Thus this parser has its own command
//  directory with a find() method which is used to find commands.

class MatchParser: public OpalParser {

public:

    MatchParser();
    virtual ~MatchParser();

protected:

    /// Find object by name in the match command directory.
    virtual Object *find(const std::string &) const;

private:

    // Not implemented.
    MatchParser(const MatchParser &);
    void operator=(const MatchParser &);

    // The sub-command directory.
    Directory MatchDirectory;
};

#endif // OPAL_MatchParser_HH
