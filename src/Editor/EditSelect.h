#ifndef OPAL_EditSelect_HH
#define OPAL_EditSelect_HH

// ------------------------------------------------------------------------
// $RCSfile: EditSelect.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditSelect
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditSelect
// ------------------------------------------------------------------------
/// The sequence editor SELECT command.

class EditSelect: public Editor {

public:

  /// Exemplar constructor.
  EditSelect();

  virtual ~EditSelect();

  /// Make clone.
  virtual EditSelect *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  EditSelect(const EditSelect &);
  void operator=(const EditSelect &);

  // Clone constructor.
  EditSelect(const string &name, EditSelect *parent);
};

#endif // OPAL_EditSelect_HH
