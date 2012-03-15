#ifndef CLASSIC_MacroCmd_HH
#define CLASSIC_MacroCmd_HH

// ------------------------------------------------------------------------
// $RCSfile: MacroCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MacroCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "MadParser/Macro.h"
#include "MadParser/MacroStream.h"
#include "MemoryManagement/Pointer.h"
#include <iosfwd>


// Class MacroCmd
// ------------------------------------------------------------------------
//: Encapsulate the buffer for the ``archetypes'' of all macros.
//  The macro is stored as a MacroStream.  For execution, first the
//  parameters are replaced, then the resulting stream is sent to the parser.

class MacroCmd: public Macro {

public:

  MacroCmd();
  MacroCmd(const string &name, MacroCmd *parent);
  virtual ~MacroCmd();

  //: Execute the macro command.
  virtual void execute();

  //: Make a macro instance.
  //  Expects parse pointer in the statement to be set on the first argument.
  //  The parser is used to determine the parse mode
  //  (normal, error, match, edit, track).
  virtual Object *makeInstance
  (const string &name, Statement &, const Parser *);

  //: Make a macro template.
  //  Expects parse pointer in the statement to be set on the first argument.
  virtual Object *makeTemplate(const string &, TokenStream &, Statement &);

private:

  // Not implemented.
  MacroCmd(const MacroCmd &);
  void operator=(const MacroCmd &);

  // The stream of tokens representing the macro command.
  Pointer<MacroStream> body;

  // Pointer to the parser to be used in execution.
  const Parser *itsParser;
};

#endif // CLASSIC_MacroCmd_HH
