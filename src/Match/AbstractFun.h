#ifndef OPAL_AbstractFun_HH
#define OPAL_AbstractFun_HH 1

// ------------------------------------------------------------------------
// $RCSfile: AbstractFun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractFun
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algebra/Vector.h"
#include <iosfwd>
#include <string>

using std::string;


// Class AbstractFun
// ------------------------------------------------------------------------
/// Abstract base for matching constraints.
//  The interface allows a constraint to return several constrained 
//  functions.

class AbstractFun {

public:

  AbstractFun();
  virtual ~AbstractFun();

  /// Get number of constrained values.
  virtual int countConstraints() const = 0;

  /// Evaluate the matching function(s).
  //  Increment [b]n[/b] for each constrained value and store the value
  //  in vector [b]f[/b].
  virtual void evaluate(Vector<double> &f, int &n) const = 0;

  /// Print the function name and value(s).
  virtual void print(std::ostream &) const = 0;

private:

  // Not implemented.
  AbstractFun(const AbstractFun &);
  void operator=(const AbstractFun &);
};

#endif // OPAL_AbstractFun_HH

