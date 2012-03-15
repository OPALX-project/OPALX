#ifndef OPAL_EditParser_HH
#define OPAL_EditParser_HH

// ------------------------------------------------------------------------
// $RCSfile: EditParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditParser
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "AbstractObjects/Directory.h"

class Statement;


// Class EditParser
// ------------------------------------------------------------------------
/// The parser for the OPAL sequence editor.
//  As long as control remains in this class, OPAL recognises only the
//  commands allowed during sequence editing.  Thus this parser has its own
//  command directory with a find() method which is used to find commands.

class EditParser: public OpalParser {

public:

  EditParser();
  virtual ~EditParser();

protected:

  /// Find object by name in the sequence editor command directory.
  virtual Object *find(const string &) const;

  /// Parse and execute current statement.
  virtual void parse(Statement &) const;

  /// Parse definition.
  //  Special version for INSTALL command.
  virtual void parseInstall(Statement &) const;

private:

  // Not implemented.
  EditParser(const EditParser &);
  void operator=(const EditParser &);

  // The sub-command directory.
  Directory editDirectory;
};

#endif // OPAL_EditParser_HH
