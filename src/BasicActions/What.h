#ifndef OPAL_What_HH
#define OPAL_What_HH

// ------------------------------------------------------------------------
// $RCSfile: What.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: What
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class What
// ------------------------------------------------------------------------
/// The WHAT command.

class What: public Action {

public:

  /// Exemplar constructor.
  What();

  virtual ~What();

  /// Make clone.
  virtual What *clone(const string &name);

  /// Execute the command.
  virtual void execute();

  /// Parse command (special for one-attribute command).
  virtual void parse(Statement &);

private:

  // Not implemented.
  What(const What &);
  void operator=(const What &);

  // Clone constructor.
  What(const string &name, What *parent);
};

#endif // OPAL_What_HH
