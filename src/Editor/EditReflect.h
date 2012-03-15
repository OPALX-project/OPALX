#ifndef OPAL_EditReflect_HH
#define OPAL_EditReflect_HH

// ------------------------------------------------------------------------
// $RCSfile: EditReflect.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditReflect
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditReflect
// ------------------------------------------------------------------------
/// The sequence editor REFLECT command.

class EditReflect: public Editor {

public:

  /// Exemplar constructor.
  EditReflect();

  virtual ~EditReflect();

  /// Make clone.
  virtual EditReflect *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  EditReflect(const EditReflect &);
  void operator=(const EditReflect &);

  // Clone constructor.
  EditReflect(const string &name, EditReflect *parent);
};

#endif // OPAL_EditReflect_HH
