#ifndef OPAL_Match_HH
#define OPAL_Match_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Match.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Match
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/MatchParser.h"
#include "Match/MatchState.h"
#include <string>
#include <list>

using std::string;

class AbstractFun;
class AbstractVar;
class MatchParser;

template <class T> class Vector;


// Class Match
// ------------------------------------------------------------------------
/// Matching block.
//  This class encapsulates all data for matching which do not depend
//  on the matching method.
//  It acts as a communication area between the various matching commands.

class Match {

public:
  
  /// Constructor.
  Match();

  ~Match();

  /// Add a matching variable.
  void addVariable(AbstractVar *);

  /// Delete a matching variable.
  //  Identified by name.
  void deleteVariable(const string &name);

  /// Find a matching variable.
  //  Identified by name.
  AbstractVar *findVariable(const string &name);

  /// Get values of matching variables.
  void getVariables(Vector<double> &x) const;

  /// Set values of matching variables.
  void setVariables(const Vector<double> &x);

  /// Get total number of variables.
  int countVariables() const;

  /// Add a set of matching function(s).
  void addFunction(AbstractFun *);

  /// Evaluate the matching functions.
  //  Set the matching variables to  [b]x[/b], cache the function values,
  //  and return them in  [b]f[/b].  The boolean return value indicates
  //  success (true) or failure (false).
  bool evaluate(const Vector<double> &x, Vector<double> &f);

  /// Get cached values of matching functions.
  void getFunctions(Vector<double> &f) const;

  /// Return total number of functions.
  int countFunctions() const;

  /// Print the results of minimisation.
  void print(const char *method, MatchState state);

  /// Get the flag for printing.
  int getPrintLevel() const;

  /// Set the flag for printing.
  void setPrintLevel(int);

  /// Return count of function evaluations.
  int getCallCount() const;

  /// The parser used during for matching.
  MatchParser parser;

  /// The block of match data.
  static Match *block;

private:

  // Not implemented.
  Match(const Match &);
  void operator=(const Match &);

  // The list of access functions to variables.
  typedef std::list<AbstractVar *> VarList;
  VarList theVariables;
  
  // The list of matching functions.
  typedef std::list<AbstractFun *> FunList;
  FunList theFunctions;

  // Counter for matched values.
  int constraintCount;

  // Counter for functions calls.
  int callCount;

  // The flag for printing.
  int printLevel;
};

#endif // OPAL_Match_HH
