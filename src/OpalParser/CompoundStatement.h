#ifndef OPALX_CompoundStatement_HH
#define OPALX_CompoundStatement_HH 1

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

#include <iosfwd>
#include <vector>
#include "OpalParser/Statement.h"

class TokenStream;
class Parser;

// class CompoundStatement
// ------------------------------------------------------------------------
/// Compound statement.
//  A statement of the form "{ <statement>; <statement>; }"
//  (may be used for FOR, IF, or WHILE blocks).
//  Stored tokens are copied into a parser-owned stream for each execution.

class CompoundStatement : public Statement {
public:
    /// Constructor.
    //  Parse the statement on the given token stream.
    CompoundStatement(TokenStream&);

    virtual ~CompoundStatement();

    /// Execute.
    //  Use the given parser to execute the contained statements.
    virtual void execute(const Parser&);

private:
    // Not implemented.
    CompoundStatement();
    CompoundStatement(const CompoundStatement&);
    void operator=(const CompoundStatement&);

    // Token list contained in this compound.
    std::vector<Token> tokens;
};

#endif  // OPALX_CompoundStatement_HH
