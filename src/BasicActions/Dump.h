#ifndef OPAL_Dump_HH
#define OPAL_Dump_HH

// ------------------------------------------------------------------------
// $RCSfile: Dump.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Dump
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Dump
// ------------------------------------------------------------------------
/// The DUMP command.

class Dump: public Action {

public:

  /// Exemplar constructor.
  Dump();

  virtual ~Dump();

  /// Make clone.
  virtual Dump *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  Dump(const Dump &);
  void operator=(const Dump &);

  // Clone constructor.
  Dump(const string &name, Dump *parent);
};

#endif // OPAL_Dump_H
