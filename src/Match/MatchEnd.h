#ifndef OPAL_MatchEnd_HH
#define OPAL_MatchEnd_HH

// ------------------------------------------------------------------------
// $RCSfile: MatchEnd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchEnd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

class MatchParser;


// Class MatchEnd
// ------------------------------------------------------------------------
/// The ENDMATCH command.

class MatchEnd: public Action {

public:

  /// Exemplar constructor.
  MatchEnd();

  virtual ~MatchEnd();

  /// Make clone.
  virtual MatchEnd *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  MatchEnd(const MatchEnd &);
  void operator=(const MatchEnd &);

  // Clone constructor.
  MatchEnd(const string &name, MatchEnd *parent);
};

#endif // OPAL_MatchEnd_HH
