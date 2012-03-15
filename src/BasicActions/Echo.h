#ifndef OPAL_Echo_HH
#define OPAL_Echo_HH

// ------------------------------------------------------------------------
// $RCSfile: Echo.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Echo
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Echo
// ------------------------------------------------------------------------
/// The ECHO command.

class Echo: public Action {

public:

  /// Exemplar constructor.
  Echo();

  virtual ~Echo();

  /// Make clone.
  virtual Echo *clone(const string &name);

  /// Execute the command.
  virtual void execute();

  /// Parse command (special for one-attribute command).
  virtual void parse(Statement &);

private:

  // Not implemented.
  Echo(const Echo &);
  void operator=(const Echo &);

  // Clone constructor.
  Echo(const string &name, Echo *parent);
};

#endif // OPAL_Echo_HH
