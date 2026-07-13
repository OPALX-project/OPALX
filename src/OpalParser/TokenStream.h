#ifndef OPALX_TokenStream_HH
#define OPALX_TokenStream_HH

// ------------------------------------------------------------------------
// $RCSfile: TokenStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TokenStream
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <string>
#include "OpalParser/Token.h"

// Class TokenStream
// ------------------------------------------------------------------------
/// Abstract interface for a stream of input tokens.

class TokenStream {
public:
    /// Constructor.
    //  Store the stream name.
    TokenStream(const std::string& name);

    virtual ~TokenStream();

    /// Put token back to stream.
    //  This allows to reparse the token.
    void putBack(const Token& token);

    /// Read single token from stream.
    virtual Token readToken() = 0;

    /// Return stream name.
    const std::string& getName() const;

    /// Return line number.
    int getLine() const;

protected:
    // Current stream name.
    std::string stream_name;

    // Current input line number.
    int curr_line;

    // Last token put back, if any.
    bool put_back_flag;
    Token put_back;

private:
    // Not implemented.
    TokenStream();
    TokenStream(const TokenStream&);
    void operator=(const TokenStream&);
};

#endif  // OPALX_TokenStream_HH
