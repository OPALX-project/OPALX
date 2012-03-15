#ifndef OPAL_SequenceTemplate_HH
#define OPAL_SequenceTemplate_HH

// ------------------------------------------------------------------------
// $RCSfile: SequenceTemplate.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SequenceTemplate
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:41:42 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "OpalParser/Macro.h"
#include "Parser/SimpleStatement.h"
#include "OpalParser/MacroStream.h"
#include <list>

class Sequence;
class Statement;
class TokenStream;


// Class SequenceTemplate
// ------------------------------------------------------------------------
/// An ``archetype'' for a SEQUENCE with arguments.
//  The model is stored in form of a MacroStream.  A call to the macro
//  sequence is expanded by first replacing the arguments, and then parsing
//  the resulting stream as a SEQUENCE definition.

class SequenceTemplate: public Macro {

  friend class Sequence;

public:

  SequenceTemplate();
  virtual ~SequenceTemplate();

  /// Make clone.
  //  Throw OpalException, since the template cannot be cloned.
  virtual SequenceTemplate *clone(const string &name);

  /// Make line instance.
  //  The instance gets the name [b]name[/b], and its actual arguments
  //  are read from [b]stat[/b].  The parser is ignored.
  virtual Object *makeInstance
  (const string &name, Statement &, const Parser *);

  /// Make a sequence template.
  //  Return NULL, since one cannot make a template from a template.
  virtual Object *makeTemplate(const string &name, TokenStream &, Statement &);

  /// Parse the sequence template.
  void parseTemplate(TokenStream &, Statement &);

private:

  // Not implemented.
  SequenceTemplate(const SequenceTemplate &);
  void operator=(const SequenceTemplate &);

  // Clone constructor.
  SequenceTemplate(const string &name, Object *parent);

  // The contained beam sequence element list.
  MacroStream body;
};

#endif // OPAL_SequenceTemplate_HH
