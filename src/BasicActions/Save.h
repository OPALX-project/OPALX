#ifndef OPAL_Save_HH
#define OPAL_Save_HH

// ------------------------------------------------------------------------
// $RCSfile: Save.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Save
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Save
// ------------------------------------------------------------------------
/// The SAVE command.

class Save: public Action {

public:

  /// Exemplar constructor.
  Save();

  virtual ~Save();

  /// Make clone.
  virtual Save *clone(const string &name);

  /// Execute the command.
  virtual void execute();

  /// Parse command (special for one-attribute command).
 void parse(Statement &);

private:

  // Not implemented.
  Save(const Save &);
  void operator=(const Save &);

  // Clone constructor.
  Save(const string &name, Save *parent);
};

#endif // OPAL_Save_H
