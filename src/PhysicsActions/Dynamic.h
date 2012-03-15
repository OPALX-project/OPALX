#ifndef OPAL_Dynamic_HH
#define OPAL_Dynamic_HH

// ------------------------------------------------------------------------
// $RCSfile: Dynamic.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Dynamic
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:19:44 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Dynamic.
// ------------------------------------------------------------------------
/// The DYNAMIC command.

class Dynamic: public Action {

public:

  /// Exemplar constructor.
  Dynamic();

  virtual ~Dynamic();

  /// Make clone.
  virtual Dynamic *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  Dynamic(const Dynamic &);
  void operator=(const Dynamic &);

  // Clone constructor.
  Dynamic(const string &name, Dynamic *parent);
};

#endif // OPAL_Dynamic_HH
