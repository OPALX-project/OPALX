#ifndef OPAL_EditCycle_HH
#define OPAL_EditCycle_HH

// ------------------------------------------------------------------------
// $RCSfile: EditCycle.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditCycle
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditCycle
// ------------------------------------------------------------------------
/// The sequence editor CYCLE command.

class EditCycle: public Editor {

public:

  /// Exemplar constructor.
  EditCycle();

  virtual ~EditCycle();

  /// Make clone.
  virtual EditCycle *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  EditCycle(const EditCycle &);
  void operator=(const EditCycle &);

  // Clone constructor.
  EditCycle(const string &name, EditCycle *parent);
};

#endif // OPAL_EditCycle_HH
