#ifndef CLASSIC_TerminalStream_HH
#define CLASSIC_TerminalStream_HH

// ------------------------------------------------------------------------
// $RCSfile: TerminalStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TerminalStream
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/AbsFileStream.h"


// Class TerminalStream
// ------------------------------------------------------------------------
/// A stream of input tokens.
//  The source of tokens is the terminal.

class TerminalStream: public AbsFileStream {

public:

  /// Constructor.
  //  The C-style string program name may be used in the readline package.
  TerminalStream(const char program[]);

  virtual ~TerminalStream();

  /// Read next input line.
  virtual bool fillLine();

private:

  // Not implemented.
  TerminalStream(const TerminalStream &);
  void operator=(const TerminalStream &);
};

#endif // CLASSIC_TerminalStream_HH
