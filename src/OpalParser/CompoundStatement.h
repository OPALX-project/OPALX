#ifndef CLASSIC_CompoundStatement_HH
#define CLASSIC_CompoundStatement_HH 1

// ------------------------------------------------------------------------
// $RCSfile: CompoundStatement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CompoundStatement
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/Statement.h"
#include "OpalParser/MacroStream.h"
#include <memory>
#include <iosfwd>

class TokenStream;
class Parser;


// class CompoundStatement
// ------------------------------------------------------------------------
/// Compound statement.
//  A statement of the form "{ <statement>; <statement>; }"
//  (may be used for FOR, IF, or WHILE blocks).
//  The compound statement is stored as a MacroStream which is sent to
//  the parser for execution.

class CompoundStatement: public Statement {

public:

    /// Constructor.
    //  Parse the statement on the given token stream.
    CompoundStatement(TokenStream &);

    virtual ~CompoundStatement();

    /// Execute.
    //  Use the given parser to execute the contained statements.
    virtual void execute(const Parser &);

private:

    // Not implemented.
    CompoundStatement();
    CompoundStatement(const CompoundStatement &);
    void operator=(const CompoundStatement &);

    // Token list contained in this compound.
    std::shared_ptr<MacroStream> tokens;
};

#endif // CLASSIC_CompoundStatement_HH
