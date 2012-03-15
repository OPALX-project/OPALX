#ifndef OPAL_System_HH
#define OPAL_System_HH

// ------------------------------------------------------------------------
// $RCSfile: System.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: System
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class System
// ------------------------------------------------------------------------
/// The SYSTEM command.

class System: public Action {

public:

  /// Exemplar constructor.
  System();

  virtual ~System();

  /// Make clone.
  virtual System *clone(const string &name);

  /// Execute the command.
  virtual void execute();

  /// Parse command (special for one-attribute command).
  virtual void parse(Statement &);

private:

  // Not implemented.
  System(const System &);
  void operator=(const System &);

  // Clone constructor.
  System(const string &name, System *parent);
};

#endif // OPAL_System_HH
