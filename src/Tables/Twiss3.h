#ifndef OPAL_Twiss3_HH
#define OPAL_Twiss3_HH

// ------------------------------------------------------------------------
// $RCSfile: Twiss3.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Twiss3
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include <iosfwd>

class Twiss;


// Class Twiss3
// ------------------------------------------------------------------------
/// The TWISS3 command.

class Twiss3: public Action {

public:

  /// Exemplar constructor.
  Twiss3();

  virtual ~Twiss3();

  /// Make clone.
  virtual Twiss3 *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  Twiss3(const Twiss3 &);
  void operator=(const Twiss3 &);

  // Clone constructor.
  Twiss3(const string &name, Twiss3 *parent);

  // Do the listing.
  void format(std::ostream &, const Twiss *);

  /// Print Twiss table in Mais-Ripken representation.
  void formatPrint(std::ostream &, const Twiss *) const;

  /// Print Twiss table in Mais-Ripken representation.
  void formatTFS(std::ostream &, const Twiss *) const;
};

#endif // OPAL_Twiss3_HH
