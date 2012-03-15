#ifndef OPAL_EditRemove_HH
#define OPAL_EditRemove_HH

// ------------------------------------------------------------------------
// $RCSfile: EditRemove.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditRemove
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditRemove
// ------------------------------------------------------------------------
/// The sequence editor REMOVE command.

class EditRemove: public Editor {

public:

  /// Exemplar constructor.
  EditRemove();

  virtual ~EditRemove();

  /// Make clone.
  virtual EditRemove *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  EditRemove(const EditRemove &);
  void operator=(const EditRemove &);

  // Clone constructor.
  EditRemove(const string &name, EditRemove *parent);
};

#endif // OPAL_EditRemove_HH
