#ifndef OPAL_ErrorCmd_HH
#define OPAL_ErrorCmd_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ErrorCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class ErrorCmd
// ------------------------------------------------------------------------
/// The ERROR command.

class ErrorCmd: public Action {

public:

  /// Exemplar constructor.
  ErrorCmd();

  virtual ~ErrorCmd();

  /// Make clone.
  virtual ErrorCmd *clone(const string &name);
  
  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  ErrorCmd(const ErrorCmd &);
  void operator=(const ErrorCmd &);

  // Clone constructor.
  ErrorCmd(const string &name, ErrorCmd *parent);
};

#endif // OPAL_ErrorCmd_HH
