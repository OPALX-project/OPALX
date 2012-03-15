#ifndef OPAL_MatchOption_HH
#define OPAL_MatchOption_HH

// ------------------------------------------------------------------------
// $RCSfile: MatchOption.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchOption
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

class Sequence;
class MatchParser;


// Class MatchOption
// ------------------------------------------------------------------------
/// The OPTION command for matching.

class MatchOption: public Action {

public:

  /// Exemplar constructor.
  MatchOption();

  virtual ~MatchOption();

  /// Make clone.
  virtual MatchOption *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  MatchOption(const MatchOption &);
  void operator=(const MatchOption &);

  // Clone constructor.
  MatchOption(const string &name, MatchOption *parent);
};

#endif // OPAL_MatchOption_HH
